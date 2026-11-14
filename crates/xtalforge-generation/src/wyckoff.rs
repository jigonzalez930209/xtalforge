#![allow(clippy::indexing_slicing)] // fixed three-component vectors and bounded orbit arrays

//! Wyckoff-constrained random crystals (Task 20.1.1).
//!
//! The generator chooses a combination of Wyckoff positions whose
//! multiplicities add up to the requested number of atoms per conventional
//! cell, evaluates each coordinate template at random free parameters, and
//! builds one representative site per orbit. The structure keeps the Hall
//! symbol and the explicit operations of the requested group, so the normal
//! P1 expansion produces the full cell.
//!
//! A candidate is accepted only when its periodic minimum distance is at
//! least `min_distance` and a symmetry search on the expanded cell returns the
//! requested ITA number. Special free parameters can collapse an orbit onto a
//! smaller one or add an accidental symmetry; those candidates are rejected
//! and the next attempt draws new parameters.
//!
//! Every attempt consumes the seeded SplitMix64 stream in a fixed order,
//! so the same seed, composition and group give the same crystal.

use std::collections::BTreeSet;

use serde::{Deserialize, Serialize};
use serde_json::{Value, json};
use xtalforge_core::{Error, Site, Structure, StructureBuilder};
use xtalforge_symmetry::{
    SpaceGroup, SymmetryOptions, WyckoffPosition, crystal_system, detect_symmetry,
};

use crate::geometry::{evaluate_template, lattice_for_system, minimum_distance, orbit, wrap};
use crate::rng::SplitMix64;

/// Options for [`generate_crystal`].
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct CrystalOptions {
    /// ITA space-group number (`1` to `230`).
    pub space_group: u32,
    /// Element symbols and their counts in the conventional cell.
    pub composition: Vec<(String, u32)>,
    /// Target volume per atom in cubic Angstrom.
    pub volume_per_atom: f64,
    /// Smallest allowed periodic distance in Angstrom.
    pub min_distance: f64,
    /// Number of attempts before the generator reports failure.
    pub max_attempts: usize,
    /// Explicit seed.
    pub seed: u64,
}

impl CrystalOptions {
    /// Options for one composition in a space group with the documented
    /// defaults: 20 A^3 per atom, 1.0 A minimum distance, 64 attempts.
    #[must_use]
    pub fn new(space_group: u32, composition: Vec<(String, u32)>, seed: u64) -> Self {
        Self {
            space_group,
            composition,
            volume_per_atom: 20.0,
            min_distance: 1.0,
            max_attempts: 64,
            seed,
        }
    }
}

/// Result of one generation.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct GeneratedCrystal {
    /// The generated structure with its space group and operations.
    pub structure: Structure,
    /// Requested ITA number.
    pub space_group: u32,
    /// Wyckoff letters used, one entry per orbit.
    pub wyckoffs: Vec<String>,
    /// Attempt that produced the crystal (1-based).
    pub attempt: usize,
    /// Seed that produced the crystal.
    pub seed: u64,
}

impl GeneratedCrystal {
    /// Serializes the result to JSON.
    ///
    /// # Errors
    ///
    /// Returns [`Error::InvalidInput`] when serialization fails.
    pub fn to_json(&self) -> Result<String, Error> {
        serde_json::to_string(self)
            .map_err(|error| Error::InvalidInput(format!("could not serialize crystal: {error}")))
    }

    /// One-line human summary.
    #[must_use]
    pub fn summary(&self) -> String {
        format!(
            "{} in {} (seed {}, attempt {}, Wyckoff {})",
            self.structure.name(),
            self.space_group,
            self.seed,
            self.attempt,
            self.wyckoffs.join(" ")
        )
    }
}

/// Generates one random crystal.
///
/// # Errors
///
/// Returns [`Error::InvalidInput`] for an invalid space group, composition,
/// volume or attempt count, and a final [`Error::InvalidInput`] when no
/// attempt satisfies the distance and symmetry constraints.
pub fn generate_crystal(options: &CrystalOptions) -> Result<GeneratedCrystal, Error> {
    validate(options)?;
    let group = SpaceGroup::from_number(options.space_group)?;
    let positions = WyckoffPosition::list(group.hall_number)?;
    let total: u32 = options.composition.iter().map(|(_, count)| *count).sum();
    let volume = options.volume_per_atom * f64::from(total);
    let mut rng = SplitMix64::new(options.seed);
    let mut last_reason = String::from("no candidate was tried");
    for attempt in 1..=options.max_attempts {
        let lattice = lattice_for_system(crystal_system(group.number), volume, &mut rng)?;
        let multiplicities: Vec<(usize, char, bool)> = positions
            .iter()
            .map(|position| {
                (
                    position.multiplicity,
                    position.letter,
                    position.coordinates.contains(['x', 'y', 'z']),
                )
            })
            .collect();
        let Some(selection) = choose_orbit_letters(&options.composition, &multiplicities, &mut rng)
        else {
            last_reason = "the composition cannot be expressed with the Wyckoff \
                multiplicities of this group"
                .to_owned();
            continue;
        };
        match build_candidate(options, &group, &positions, &selection, lattice, &mut rng) {
            Ok(Some((structure, wyckoffs))) => {
                return Ok(GeneratedCrystal {
                    structure,
                    space_group: group.number,
                    wyckoffs,
                    attempt,
                    seed: options.seed,
                });
            }
            Ok(None) => {
                last_reason =
                    "a candidate collapsed its orbit or violated the minimum distance".to_owned();
            }
            Err(error) => last_reason = error.message(),
        }
    }
    Err(Error::InvalidInput(format!(
        "could not generate a {} crystal after {} attempts: {last_reason}",
        options.space_group, options.max_attempts
    )))
}

/// Generates one crystal per space-group number in `groups`.
///
/// Failed groups are reported in place, so the output length always matches
/// the input length. The seeds advance deterministically by group index.
///
/// # Errors
///
/// Returns the validation errors of [`generate_crystal`] for every element.
pub fn generate_crystal_range(
    groups: &[u32],
    options: &CrystalOptions,
) -> Result<Vec<Result<GeneratedCrystal, Error>>, Error> {
    validate(options)?;
    let mut out = Vec::with_capacity(groups.len());
    for (offset, group) in groups.iter().enumerate() {
        let mut per_group = options.clone();
        per_group.space_group = *group;
        per_group.seed = options.seed.wrapping_add(offset as u64);
        out.push(generate_crystal(&per_group));
    }
    Ok(out)
}

/// Number of structurally distinct random crystals for one seed range.
///
/// Used by the statistics documented in the algorithm page; it counts the
/// crystals whose detected space group equals the request.
///
/// # Errors
///
/// Returns the first validation error.
pub fn generation_statistics(
    options: &CrystalOptions,
    seeds: u64,
) -> Result<GenerationStatistics, Error> {
    validate(options)?;
    let mut succeeded = 0u64;
    let mut failed = 0u64;
    let mut attempts = 0u64;
    for offset in 0..seeds {
        let mut per_seed = options.clone();
        per_seed.seed = options.seed.wrapping_add(offset);
        match generate_crystal(&per_seed) {
            Ok(crystal) => {
                succeeded += 1;
                attempts += crystal.attempt as u64;
            }
            Err(_) => failed += 1,
        }
    }
    Ok(GenerationStatistics {
        seeds,
        succeeded,
        failed,
        mean_attempts: if succeeded == 0 {
            0.0
        } else {
            attempts as f64 / succeeded as f64
        },
    })
}

/// Outcome of a seed sweep.
#[derive(Debug, Clone, Copy, PartialEq, Serialize, Deserialize)]
pub struct GenerationStatistics {
    /// Number of seeds tried.
    pub seeds: u64,
    /// Crystals that satisfied the acceptance checks.
    pub succeeded: u64,
    /// Seeds that exhausted their attempts.
    pub failed: u64,
    /// Mean attempts per successful crystal.
    pub mean_attempts: f64,
}

fn validate(options: &CrystalOptions) -> Result<(), Error> {
    if !(1..=230).contains(&options.space_group) {
        return Err(Error::InvalidInput(format!(
            "space-group number must be in 1..=230, got {}",
            options.space_group
        )));
    }
    if options.composition.is_empty() {
        return Err(Error::InvalidInput(
            "the composition needs at least one element".to_owned(),
        ));
    }
    let mut seen = BTreeSet::new();
    for (element, count) in &options.composition {
        let info = xtalforge_core::elements::element_by_symbol(element)?;
        if !seen.insert(info.symbol.clone()) {
            return Err(Error::InvalidInput(format!(
                "element `{element}` is listed twice"
            )));
        }
        if *count == 0 {
            return Err(Error::InvalidInput(format!(
                "element `{element}` has a zero count"
            )));
        }
    }
    let total: u32 = options.composition.iter().map(|(_, count)| *count).sum();
    if total > 400 {
        return Err(Error::InvalidInput(format!(
            "the conventional cell is capped at 400 atoms, got {total}"
        )));
    }
    if !options.volume_per_atom.is_finite() || options.volume_per_atom <= 0.0 {
        return Err(Error::InvalidInput(
            "volume per atom must be positive and finite".to_owned(),
        ));
    }
    if !options.min_distance.is_finite() || options.min_distance < 0.0 {
        return Err(Error::InvalidInput(
            "minimum distance must be finite and non-negative".to_owned(),
        ));
    }
    if options.max_attempts == 0 {
        return Err(Error::InvalidInput(
            "max_attempts must be at least one".to_owned(),
        ));
    }
    Ok(())
}

/// One chosen orbit: position index in the group list and element.
pub(crate) struct Selection {
    pub(crate) element: String,
    pub(crate) position: usize,
}

pub(crate) fn choose_orbit_letters(
    composition: &[(String, u32)],
    positions: &[(usize, char, bool)],
    rng: &mut SplitMix64,
) -> Option<Vec<Selection>> {
    // A composition may need a common scaling factor to match the orbit
    // multiplicities (for example one formula unit in a face-centred cell).
    for scale in 1..=6u32 {
        let targets: Vec<(String, usize)> = composition
            .iter()
            .map(|(element, count)| (element.clone(), (*count * scale) as usize))
            .collect();
        let total: usize = targets.iter().map(|(_, count)| *count).sum();
        if total > 400 {
            break;
        }
        let mut combinations: Vec<Vec<Vec<usize>>> = Vec::new();
        let mut possible = true;
        for (_, target) in &targets {
            let found = position_combinations(*target, positions);
            if found.is_empty() {
                possible = false;
                break;
            }
            combinations.push(found);
        }
        if !possible {
            continue;
        }
        let order: Vec<usize> = (0..targets.len()).collect();
        for _ in 0..96 {
            let mut used = BTreeSet::new();
            let mut selection = Vec::new();
            let mut ok = true;
            for element_index in &order {
                let Some((element, _)) = targets.get(*element_index) else {
                    ok = false;
                    break;
                };
                let Some(candidates) = combinations.get(*element_index) else {
                    ok = false;
                    break;
                };
                let candidate = candidates.get(rng.below(candidates.len()))?;
                for position in candidate {
                    let (_, letter, flexible) = *positions.get(*position)?;
                    if !used.insert(letter) && !flexible {
                        ok = false;
                        break;
                    }
                }
                if !ok {
                    break;
                }
                for position in candidate {
                    selection.push(Selection {
                        element: element.clone(),
                        position: *position,
                    });
                }
            }
            if ok && !selection.is_empty() {
                return Some(selection);
            }
        }
    }
    None
}

/// Combinations of at most four distinct positions whose multiplicities sum
/// to `target`.
fn position_combinations(target: usize, positions: &[(usize, char, bool)]) -> Vec<Vec<usize>> {
    let mut out = Vec::new();
    let mut current = Vec::new();
    collect(target, positions, 0, 4, &mut current, &mut out);
    // Deterministic order: fewer orbits first, then by position index.
    out.sort_by(|left, right| left.len().cmp(&right.len()).then_with(|| left.cmp(right)));
    out.truncate(256);
    out
}

fn collect(
    remaining: usize,
    positions: &[(usize, char, bool)],
    start: usize,
    depth: usize,
    current: &mut Vec<usize>,
    out: &mut Vec<Vec<usize>>,
) {
    if remaining == 0 {
        out.push(current.clone());
        return;
    }
    if depth == 0 {
        return;
    }
    for index in start..positions.len() {
        let Some(position) = positions.get(index) else {
            continue;
        };
        if position.0 > remaining {
            continue;
        }
        current.push(index);
        collect(
            remaining - position.0,
            positions,
            index + 1,
            depth - 1,
            current,
            out,
        );
        current.pop();
    }
}

fn build_candidate(
    options: &CrystalOptions,
    group: &SpaceGroup,
    positions: &[WyckoffPosition],
    selection: &[Selection],
    lattice: xtalforge_core::Lattice,
    rng: &mut SplitMix64,
) -> Result<Option<(Structure, Vec<String>)>, Error> {
    let mut builder = StructureBuilder::new(format!("random-{}", group.number), lattice)
        .hall(group.hall_symbol.clone())
        .metadata("generator", json!("wyckoff"))
        .metadata("seed", json!(options.seed))
        .metadata("space_group", json!(group.number));
    for operation in &group.operations {
        builder = builder.operation(*operation);
    }
    let mut letters = Vec::new();
    let mut label_index = 0usize;
    for chosen in selection {
        let Some(position) = positions.get(chosen.position) else {
            return Ok(None);
        };
        let free = [rng.uniform(), rng.uniform(), rng.uniform()];
        let point = wrap(evaluate_template(&position.coordinates, free)?);
        let images = orbit(&[point], &group.operations, &lattice);
        if images.len() != position.multiplicity {
            return Ok(None);
        }
        label_index += 1;
        let site = Site::new(
            format!("{}{label_index}", chosen.element),
            &chosen.element,
            point,
        )?;
        builder = builder.site(site);
        letters.push(position.letter.to_string());
    }
    let structure = builder.build()?;
    if minimum_distance(&structure)? < options.min_distance {
        return Ok(None);
    }
    let dataset = detect_symmetry(&structure, SymmetryOptions::default())?;
    if dataset.number != group.number {
        return Ok(None);
    }
    Ok(Some((structure, letters)))
}

/// JSON description of one generated crystal for the bindings.
///
/// # Errors
///
/// Returns [`Error::InvalidInput`] when serialization fails.
pub fn crystal_json(crystal: &GeneratedCrystal) -> Result<Value, Error> {
    let value = serde_json::to_value(crystal)
        .map_err(|error| Error::InvalidInput(format!("could not serialize crystal: {error}")))?;
    Ok(value)
}
