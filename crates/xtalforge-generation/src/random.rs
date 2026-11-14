//! AIRSS-style random search (Task 20.2.1).
//!
//! The search draws random symmetric crystals, relaxes them with a native
//! [`Calculator`], optionally removes duplicates with the Task 16.3.2
//! deduplicator and ranks the survivors by energy per atom. The whole run is
//! driven by one seed: the i-th structure uses `seed + i`, so a search can be
//! reproduced and resumed from its checkpoint (Task 20.2.3).
//!
//! The default calculator is UFF4MOF. Any calculator can be supplied through
//! [`SearchCalculator`], which also covers the ML bridge and engine-backed
//! calculators of Task 18.4 when the caller owns one.

use serde::{Deserialize, Serialize};
use xtalforge_analysis::dataset::{Dataset, DatasetRecord};
use xtalforge_analysis::dedup::{DedupOptions, deduplicate};
use xtalforge_core::{Error, Structure};
use xtalforge_ff::{
    Calculator, CalculatorOptions, Minimizer, OptimizeOptions, UffCalculator, minimize,
};
use xtalforge_symmetry::MatcherOptions;

use crate::wyckoff::{CrystalOptions, generate_crystal};

/// Which minimizer the search uses.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default, Serialize, Deserialize)]
#[serde(rename_all = "snake_case")]
pub enum SearchMinimizer {
    /// FIRE.
    #[default]
    Fire,
    /// L-BFGS.
    Lbfgs,
}

/// Options for [`random_search`].
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct RandomSearchOptions {
    /// Element symbols and counts in the conventional cell.
    pub composition: Vec<(String, u32)>,
    /// Space groups to sample, in round-robin order.
    pub space_groups: Vec<u32>,
    /// Number of structures to generate.
    pub n_structures: usize,
    /// Explicit seed.
    pub seed: u64,
    /// Volume per atom for the generator.
    pub volume_per_atom: f64,
    /// Minimum distance for the generator.
    pub min_distance: f64,
    /// Attempts per generated structure.
    pub max_attempts: usize,
    /// Whether to relax every structure with UFF.
    pub relax: bool,
    /// Whether the cell is relaxed too.
    pub relax_cell: bool,
    /// Maximum optimizer steps.
    pub max_steps: usize,
    /// Whether to remove duplicates before ranking.
    pub deduplicate: bool,
    /// Matcher tolerances for deduplication.
    pub matcher: MatcherOptions,
    /// Minimizer.
    pub minimizer: SearchMinimizer,
}

impl RandomSearchOptions {
    /// A UFF search over `space_groups` with one structure per group.
    #[must_use]
    pub fn new(
        composition: Vec<(String, u32)>,
        space_groups: Vec<u32>,
        n_structures: usize,
        seed: u64,
    ) -> Self {
        Self {
            composition,
            space_groups,
            n_structures,
            seed,
            volume_per_atom: 20.0,
            min_distance: 1.0,
            max_attempts: 32,
            relax: true,
            relax_cell: false,
            max_steps: 200,
            deduplicate: true,
            matcher: MatcherOptions::default(),
            minimizer: SearchMinimizer::Fire,
        }
    }
}

/// One ranked search result.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct SearchHit {
    /// Relaxed structure.
    pub structure: Structure,
    /// Potential energy in eV.
    pub energy_ev: f64,
    /// Energy per atom in eV.
    pub energy_per_atom: f64,
    /// Space group of the generated candidate.
    pub space_group: u32,
    /// Index of the candidate in generation order.
    pub candidate: usize,
    /// Seed used for this candidate.
    pub seed: u64,
    /// Whether the structure converged.
    pub converged: bool,
    /// Whether the structure is the representative of a duplicate group.
    pub duplicate: bool,
}

/// Result of a random search.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct RandomSearchReport {
    /// Ranked hits, lowest energy per atom first.
    pub hits: Vec<SearchHit>,
    /// Number of generated candidates.
    pub generated: usize,
    /// Number of candidates that relaxed successfully.
    pub relaxed: usize,
    /// Number of duplicates removed.
    pub duplicates_removed: usize,
    /// Calculator name.
    pub calculator: String,
    /// Seed.
    pub seed: u64,
}

impl RandomSearchReport {
    /// Serializes the report to JSON.
    ///
    /// # Errors
    ///
    /// Returns [`Error::InvalidInput`] when serialization fails.
    pub fn to_json(&self) -> Result<String, Error> {
        serde_json::to_string(self).map_err(|error| {
            Error::InvalidInput(format!("could not serialize search report: {error}"))
        })
    }

    /// The best hit, when the search produced one.
    #[must_use]
    pub fn best(&self) -> Option<&SearchHit> {
        self.hits.first()
    }
}

/// A calculator supplied by the caller.
pub struct SearchCalculator<'a> {
    /// Calculator instance.
    pub calculator: &'a dyn Calculator,
    /// Name reported in the result.
    pub name: String,
}

/// Runs an AIRSS-style random search with the default UFF4MOF calculator.
///
/// # Errors
///
/// Returns the generator, optimizer or deduplication errors.
pub fn random_search(options: &RandomSearchOptions) -> Result<RandomSearchReport, Error> {
    let calculator = UffCalculator::new(CalculatorOptions::default()).map_err(crate::ff_error)?;
    random_search_with(
        options,
        &SearchCalculator {
            calculator: &calculator,
            name: "uff".to_owned(),
        },
    )
}

/// Runs an AIRSS-style random search with an explicit calculator.
///
/// # Errors
///
/// Returns [`Error::InvalidInput`] for invalid options, and the generator,
/// optimizer or deduplication errors.
pub fn random_search_with(
    options: &RandomSearchOptions,
    calculator: &SearchCalculator<'_>,
) -> Result<RandomSearchReport, Error> {
    validate(options)?;
    let optimizer = OptimizeOptions {
        max_steps: options.max_steps,
        cell: options.relax_cell,
        ..OptimizeOptions::default()
    };
    let minimizer = match options.minimizer {
        SearchMinimizer::Fire => Minimizer::Fire,
        SearchMinimizer::Lbfgs => Minimizer::Lbfgs,
    };
    let mut candidates = Vec::with_capacity(options.n_structures);
    let mut relaxed = 0usize;
    for index in 0..options.n_structures {
        let group = options
            .space_groups
            .get(index % options.space_groups.len().max(1))
            .copied()
            .unwrap_or(1);
        let seed = options.seed.wrapping_add(index as u64);
        let crystal = generate_crystal(&CrystalOptions {
            space_group: group,
            composition: options.composition.clone(),
            volume_per_atom: options.volume_per_atom,
            min_distance: options.min_distance,
            max_attempts: options.max_attempts,
            seed,
        })?;
        let (structure, energy, converged) = if options.relax {
            let result = minimize(
                calculator.calculator,
                &crystal.structure,
                &optimizer,
                minimizer,
            )
            .map_err(crate::ff_error)?;
            relaxed += 1;
            (result.structure, result.energy, result.converged)
        } else {
            let energy = calculator
                .calculator
                .energy(&crystal.structure)
                .map_err(crate::ff_error)?;
            (crystal.structure, energy, false)
        };
        candidates.push(SearchHit {
            energy_per_atom: energy / atom_count(&structure)?.max(1) as f64,
            structure,
            energy_ev: energy,
            space_group: group,
            candidate: index,
            seed,
            converged,
            duplicate: false,
        });
    }
    let mut duplicates_removed = 0usize;
    if options.deduplicate && candidates.len() > 1 {
        let mut dataset = Dataset::new("random-search");
        for (index, hit) in candidates.iter().enumerate() {
            dataset.push(DatasetRecord::new(
                format!("candidate-{index}"),
                hit.structure.clone(),
            )?)?;
        }
        let report = deduplicate(
            &dataset,
            &DedupOptions {
                matcher: options.matcher,
                ..DedupOptions::default()
            },
        )?;
        let kept: Vec<usize> = report.report.kept.clone();
        duplicates_removed = report.report.n_removed();
        let mut deduplicated = Vec::with_capacity(kept.len());
        for index in kept {
            if let Some(mut hit) = candidates.get(index).cloned() {
                hit.duplicate = false;
                deduplicated.push(hit);
            }
        }
        // Mark the removed entries if the caller asks for them later; the
        // ranked report only keeps representatives.
        candidates = deduplicated;
    }
    candidates.sort_by(|left, right| {
        left.energy_per_atom
            .total_cmp(&right.energy_per_atom)
            .then_with(|| left.candidate.cmp(&right.candidate))
    });
    let generated = options.n_structures;
    Ok(RandomSearchReport {
        hits: candidates,
        generated,
        relaxed,
        duplicates_removed,
        calculator: calculator.name.clone(),
        seed: options.seed,
    })
}

fn atom_count(structure: &Structure) -> Result<usize, Error> {
    Ok(structure.p1()?.atoms.len())
}

fn validate(options: &RandomSearchOptions) -> Result<(), Error> {
    if options.n_structures == 0 {
        return Err(Error::InvalidInput(
            "n_structures must be at least one".to_owned(),
        ));
    }
    if options.space_groups.is_empty() {
        return Err(Error::InvalidInput(
            "space_groups must list at least one group".to_owned(),
        ));
    }
    if options.max_steps == 0 {
        return Err(Error::InvalidInput(
            "max_steps must be at least one".to_owned(),
        ));
    }
    if options.composition.is_empty() {
        return Err(Error::InvalidInput(
            "the composition needs at least one element".to_owned(),
        ));
    }
    for group in &options.space_groups {
        if !(1..=230).contains(group) {
            return Err(Error::InvalidInput(format!(
                "space-group number must be in 1..=230, got {group}"
            )));
        }
    }
    Ok(())
}
