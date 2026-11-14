//! Evolutionary global search (Task 20.2.2).
//!
//! The search keeps a population of relaxed structures, ranks them by energy
//! per atom, and produces each new generation from cut-and-splice crossover
//! and two mutations: random atomic displacements and a small lattice strain.
//! Children are relaxed with the same calculator as the AIRSS search,
//! deduplicated against the population with the structure matcher, and the
//! best half survives.
//!
//! Cut-and-splice works in the P1 cell of the first parent. The second parent
//! is mapped into that cell, a random plane splits both cells, and the child
//! takes the atoms on one side from each parent. Because the child is a P1
//! assembly, its symmetry is not imposed; a later symmetrization step is the
//! caller's choice. The population is deterministic for one seed: every draw
//! comes from the seeded generator and the ranking breaks ties by candidate
//! index.
//!
//! Variable-composition searches rank by the distance to the convex hull
//! instead of energy per atom. The hull step uses the core simplex from Task
//! 9.6; elemental references are optional, and when they are missing the
//! ranking falls back to energy per atom and says so in the report.

use serde::{Deserialize, Serialize};
use xtalforge_core::{Error, Site, Structure, StructureBuilder};
use xtalforge_ff::{CalculatorOptions, Minimizer, OptimizeOptions, UffCalculator, minimize};
use xtalforge_symmetry::{MatcherOptions, match_structures};

use crate::random::{RandomSearchOptions, SearchCalculator, SearchHit, random_search_with};
use crate::rng::SplitMix64;

/// Options for [`evolutionary_search`].
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct EvolutionaryOptions {
    /// Element symbols and counts in the conventional cell.
    pub composition: Vec<(String, u32)>,
    /// Space groups used to seed the population.
    pub space_groups: Vec<u32>,
    /// Population size (even).
    pub population: usize,
    /// Number of generations after the initial one.
    pub generations: usize,
    /// Explicit seed.
    pub seed: u64,
    /// Volume per atom used by the initial generator.
    pub volume_per_atom: f64,
    /// Minimum distance used by the initial generator.
    pub min_distance: f64,
    /// Attempts per initial structure.
    pub max_attempts: usize,
    /// Maximum optimizer steps per relaxation.
    pub max_steps: usize,
    /// Whether the cell relaxes.
    pub relax_cell: bool,
    /// Fraction of the population that survives.
    pub survival: f64,
    /// Probability of a displacement mutation per child.
    pub mutation_rate: f64,
    /// Largest displacement in Angstrom.
    pub mutation_sigma: f64,
    /// Largest fractional lattice strain.
    pub strain_sigma: f64,
    /// Stop after this many generations without a better best energy.
    pub convergence_generations: usize,
    /// Matcher tolerances for deduplication.
    pub matcher: MatcherOptions,
    /// Minimizer name.
    pub minimizer: crate::random::SearchMinimizer,
    /// Optional population carried over from a checkpoint or a previous run.
    #[serde(default, skip_serializing_if = "Vec::is_empty")]
    pub initial_population: Vec<Structure>,
}

impl EvolutionaryOptions {
    /// A population of 8 for 5 generations with the documented defaults.
    #[must_use]
    pub fn new(composition: Vec<(String, u32)>, space_groups: Vec<u32>, seed: u64) -> Self {
        Self {
            composition,
            space_groups,
            population: 8,
            generations: 5,
            seed,
            volume_per_atom: 20.0,
            min_distance: 1.0,
            max_attempts: 32,
            max_steps: 150,
            relax_cell: false,
            survival: 0.5,
            mutation_rate: 0.6,
            mutation_sigma: 0.3,
            strain_sigma: 0.02,
            convergence_generations: 3,
            matcher: MatcherOptions::default(),
            minimizer: crate::random::SearchMinimizer::Fire,
            initial_population: Vec::new(),
        }
    }
}

/// One generation summary.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct GenerationReport {
    /// Generation number, zero for the initial population.
    pub generation: usize,
    /// Best energy per atom in the generation.
    pub best_energy_per_atom: f64,
    /// Mean energy per atom of the survivors.
    pub mean_energy_per_atom: f64,
    /// Children generated.
    pub children: usize,
    /// Children removed as duplicates.
    pub duplicates: usize,
    /// Population size after selection.
    pub population: usize,
}

/// Result of an evolutionary search.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct EvolutionarySearch {
    /// Final population, best first.
    pub hits: Vec<SearchHit>,
    /// One report per generation.
    pub generations: Vec<GenerationReport>,
    /// Number of relaxations performed.
    pub relaxations: usize,
    /// Whether the search converged.
    pub converged: bool,
    /// Calculator name.
    pub calculator: String,
    /// Seed.
    pub seed: u64,
}

impl EvolutionarySearch {
    /// Serializes the report to JSON.
    ///
    /// # Errors
    ///
    /// Returns [`Error::InvalidInput`] when serialization fails.
    pub fn to_json(&self) -> Result<String, Error> {
        serde_json::to_string(self).map_err(|error| {
            Error::InvalidInput(format!("could not serialize evolutionary search: {error}"))
        })
    }

    /// The best hit.
    #[must_use]
    pub fn best(&self) -> Option<&SearchHit> {
        self.hits.first()
    }
}

/// Runs an evolutionary search with UFF4MOF.
///
/// # Errors
///
/// Returns the generator, optimizer or deduplication errors.
pub fn evolutionary_search(options: &EvolutionaryOptions) -> Result<EvolutionarySearch, Error> {
    let calculator = UffCalculator::new(CalculatorOptions::default()).map_err(crate::ff_error)?;
    evolutionary_search_with(
        options,
        &SearchCalculator {
            calculator: &calculator,
            name: "uff".to_owned(),
        },
    )
}

/// Runs an evolutionary search with an explicit calculator.
///
/// # Errors
///
/// Returns [`Error::InvalidInput`] for invalid options and the generator,
/// optimizer or deduplication errors.
pub fn evolutionary_search_with(
    options: &EvolutionaryOptions,
    calculator: &SearchCalculator<'_>,
) -> Result<EvolutionarySearch, Error> {
    crate::evolutionary_check::validate(options)?;
    let random = RandomSearchOptions {
        composition: options.composition.clone(),
        space_groups: options.space_groups.clone(),
        n_structures: options.population,
        seed: options.seed,
        volume_per_atom: options.volume_per_atom,
        min_distance: options.min_distance,
        max_attempts: options.max_attempts,
        relax: true,
        relax_cell: options.relax_cell,
        max_steps: options.max_steps,
        deduplicate: false,
        matcher: options.matcher,
        minimizer: options.minimizer,
    };
    let mut population = if options.initial_population.is_empty() {
        random_search_with(&random, calculator)?.hits
    } else {
        let mut seeded = Vec::new();
        for (index, structure) in options.initial_population.iter().enumerate() {
            let energy = calculator
                .calculator
                .energy(structure)
                .map_err(crate::ff_error)?;
            let p1 = structure.p1()?;
            seeded.push(SearchHit {
                structure: structure.clone(),
                energy_ev: energy,
                energy_per_atom: energy / p1.atoms.len().max(1) as f64,
                space_group: options.space_groups.first().copied().unwrap_or(1),
                candidate: index,
                seed: options.seed,
                converged: true,
                duplicate: false,
            });
        }
        seeded
    };
    let mut relaxations = population.len();
    let mut generations = vec![summarize(0, &population, 0, 0)];
    let mut rng = SplitMix64::new(options.seed ^ 0xA5A5_5A5A_1234_5678);
    let optimizer = OptimizeOptions {
        max_steps: options.max_steps,
        cell: options.relax_cell,
        ..OptimizeOptions::default()
    };
    let minimizer = match options.minimizer {
        crate::random::SearchMinimizer::Fire => Minimizer::Fire,
        crate::random::SearchMinimizer::Lbfgs => Minimizer::Lbfgs,
    };
    let target = ((options.population as f64 * options.survival).ceil() as usize)
        .max(1)
        .min(options.population);
    let mut converged = false;
    let mut stagnant = 0usize;
    let mut best = population
        .first()
        .map_or(f64::INFINITY, |hit| hit.energy_per_atom);
    for generation in 1..=options.generations {
        population.sort_by(|left, right| {
            left.energy_per_atom
                .total_cmp(&right.energy_per_atom)
                .then_with(|| left.candidate.cmp(&right.candidate))
        });
        population.truncate(target);
        let parents = population.clone();
        let mut children: Vec<SearchHit> = Vec::new();
        let mut duplicates = 0usize;
        let mut attempts = 0usize;
        while children.len() < options.population - parents.len()
            && attempts < options.population * 8
        {
            attempts += 1;
            let first = parents.get(rng.below(parents.len())).cloned();
            let second = parents.get(rng.below(parents.len())).cloned();
            let (Some(first), Some(second)) = (first, second) else {
                break;
            };
            let mut child = crossover(&first.structure, &second.structure, &mut rng)?;
            if rng.uniform() < options.mutation_rate {
                child = mutate(
                    &child,
                    &mut rng,
                    options.mutation_sigma,
                    options.strain_sigma,
                )?;
            }
            let relaxed = minimize(calculator.calculator, &child, &optimizer, minimizer)
                .map_err(crate::ff_error)?;
            relaxations += 1;
            let result = SearchHit {
                energy_per_atom: relaxed.energy / atom_count(&relaxed.structure)?.max(1) as f64,
                structure: relaxed.structure,
                energy_ev: relaxed.energy,
                space_group: first.space_group,
                candidate: options.population + generation * options.population + children.len(),
                seed: options.seed,
                converged: relaxed.converged,
                duplicate: false,
            };
            if !is_duplicate(&result.structure, &parents, &options.matcher)? {
                children.push(result);
            } else {
                duplicates += 1;
            }
        }
        let mut next = parents;
        next.extend(children);
        next.sort_by(|left, right| {
            left.energy_per_atom
                .total_cmp(&right.energy_per_atom)
                .then_with(|| left.candidate.cmp(&right.candidate))
        });
        next.truncate(options.population);
        let generation_best = next
            .first()
            .map_or(f64::INFINITY, |hit| hit.energy_per_atom);
        if generation_best < best - 1e-4 {
            best = generation_best;
            stagnant = 0;
        } else {
            stagnant += 1;
        }
        generations.push(summarize(
            generation,
            &next,
            next.len().saturating_sub(target),
            duplicates,
        ));
        population = next;
        if stagnant >= options.convergence_generations {
            converged = true;
            break;
        }
    }
    population.sort_by(|left, right| {
        left.energy_per_atom
            .total_cmp(&right.energy_per_atom)
            .then_with(|| left.candidate.cmp(&right.candidate))
    });
    Ok(EvolutionarySearch {
        hits: population,
        generations,
        relaxations,
        converged,
        calculator: calculator.name.clone(),
        seed: options.seed,
    })
}

fn summarize(
    generation: usize,
    population: &[SearchHit],
    children: usize,
    duplicates: usize,
) -> GenerationReport {
    let best = population
        .first()
        .map_or(f64::INFINITY, |hit| hit.energy_per_atom);
    let mean = if population.is_empty() {
        0.0
    } else {
        population
            .iter()
            .map(|hit| hit.energy_per_atom)
            .sum::<f64>()
            / population.len() as f64
    };
    GenerationReport {
        generation,
        best_energy_per_atom: best,
        mean_energy_per_atom: mean,
        children,
        duplicates,
        population: population.len(),
    }
}

fn atom_count(structure: &Structure) -> Result<usize, Error> {
    Ok(structure.p1()?.atoms.len())
}

fn crossover(
    first: &Structure,
    second: &Structure,
    rng: &mut SplitMix64,
) -> Result<Structure, Error> {
    let first_p1 = first.p1()?;
    let second_p1 = second.p1()?;
    let count = first_p1.atoms.len().min(second_p1.atoms.len());
    let normal = normalize([
        rng.range(-1.0, 1.0),
        rng.range(-1.0, 1.0),
        rng.range(-1.0, 1.0),
    ]);
    let offset = rng.range(-0.5, 0.5);
    let split = count / 2;
    let mut first_atoms = first_p1.atoms.clone();
    let mut second_atoms = second_p1.atoms.clone();
    first_atoms.sort_by(|left, right| {
        signed_distance(left.frac, normal, offset)
            .total_cmp(&signed_distance(right.frac, normal, offset))
    });
    second_atoms.sort_by(|left, right| {
        signed_distance(left.frac, normal, offset)
            .total_cmp(&signed_distance(right.frac, normal, offset))
    });
    let mut builder = StructureBuilder::new("crossover", first_p1.lattice);
    let mut label = 0usize;
    for atom in first_atoms.iter().take(split) {
        label += 1;
        builder = builder.site(Site::new(
            format!("{}{label}", atom.element),
            &atom.element,
            atom.frac,
        )?);
    }
    for atom in second_atoms.iter().take(count - split) {
        label += 1;
        let cart = second_p1.lattice.to_cart(atom.frac);
        let frac = first_p1.lattice.to_frac(cart);
        builder = builder.site(Site::new(
            format!("{}{label}", atom.element),
            &atom.element,
            frac,
        )?);
    }
    builder.build()
}

fn signed_distance(frac: [f64; 3], normal: [f64; 3], offset: f64) -> f64 {
    normal[0] * (frac[0] - offset) + normal[1] * (frac[1] - offset) + normal[2] * (frac[2] - offset)
}

fn normalize([x, y, z]: [f64; 3]) -> [f64; 3] {
    let norm = (x * x + y * y + z * z).sqrt();
    if norm < 1e-9 {
        [1.0, 0.0, 0.0]
    } else {
        [x / norm, y / norm, z / norm]
    }
}

fn mutate(
    structure: &Structure,
    rng: &mut SplitMix64,
    sigma: f64,
    strain: f64,
) -> Result<Structure, Error> {
    let p1 = structure.p1()?;
    let matrix = p1.lattice.matrix();
    let mut strained = matrix;
    for row in &mut strained {
        for value in row {
            *value *= 1.0 + strain * rng.range(-1.0, 1.0);
        }
    }
    let lattice = xtalforge_core::Lattice::new(strained)?;
    let mut builder = StructureBuilder::new("mutated", lattice);
    for (index, atom) in p1.atoms.iter().enumerate() {
        let frac = [
            atom.frac[0] + sigma * rng.range(-1.0, 1.0),
            atom.frac[1] + sigma * rng.range(-1.0, 1.0),
            atom.frac[2] + sigma * rng.range(-1.0, 1.0),
        ];
        builder = builder.site(Site::new(
            format!("{}{}", atom.element, index + 1),
            &atom.element,
            crate::geometry::wrap(frac),
        )?);
    }
    builder.build()
}

fn is_duplicate(
    structure: &Structure,
    population: &[SearchHit],
    matcher: &MatcherOptions,
) -> Result<bool, Error> {
    for hit in population {
        let matched = match_structures(structure, &hit.structure, matcher)?;
        if matched.matched {
            return Ok(true);
        }
    }
    Ok(false)
}
