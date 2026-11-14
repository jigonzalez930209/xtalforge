#![allow(clippy::indexing_slicing, clippy::needless_range_loop)] // fixed-size numeric kernels

//! Special quasirandom structures (Task 9.5.3).
//!
//! A binary disordered sublattice is decorated with a fixed seed, and the
//! decoration is relaxed by swapping A/B sites so the pair and triplet
//! correlations approach the random targets. The report lists every shell.
//!
//! The multicomponent extension with pair, triplet and quadruplet targets
//! lives in [`crate::sqs_multicomponent`] and is re-exported here as
//! [`multicomponent_sqs`].

pub use crate::sqs_multicomponent::{
    MulticomponentSqsOptions, MulticomponentSqsResult, SqsCorrelation, multicomponent_sqs,
};

use std::collections::BTreeMap;

use serde::{Deserialize, Serialize};

use xtalforge_core::{Error, Site, Structure, StructureBuilder};

/// One pair shell of the SQS report.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct SqsShell {
    /// Mean distance of the shell in Angstrom.
    pub distance: f64,
    /// Number of pairs in the shell.
    pub pairs: usize,
    /// Achieved pair correlation.
    pub correlation: f64,
    /// Target pair correlation.
    pub target: f64,
}

/// Result of [`special_quasirandom`].
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct SqsResult {
    /// Decorated supercell.
    pub structure: Structure,
    /// Species per disordered position.
    pub assignment: Vec<String>,
    /// Pair correlations by shell.
    pub shells: Vec<SqsShell>,
    /// Achieved triplet correlation over nearest-neighbour triangles.
    pub triplet_correlation: f64,
    /// Weighted objective; zero is a perfect match.
    pub objective: f64,
    /// Number of accepted swaps.
    pub accepted: usize,
    /// Number of attempted swaps.
    pub attempted: usize,
}

/// Options for [`special_quasirandom`].
#[derive(Debug, Clone, PartialEq)]
pub struct SqsOptions {
    /// Supercell repeats.
    pub supercell: [u32; 3],
    /// Seed of the deterministic generator.
    pub seed: u64,
    /// Pair cutoff in Angstrom.
    pub cutoff: f64,
    /// Number of swap attempts.
    pub iterations: usize,
    /// Weight of the pair objective.
    pub pair_weight: f64,
    /// Weight of the triplet objective.
    pub triplet_weight: f64,
    /// Initial annealing temperature.
    pub temperature: f64,
}

impl Default for SqsOptions {
    fn default() -> Self {
        Self {
            supercell: [2, 2, 2],
            seed: 0,
            cutoff: 6.0,
            iterations: 20_000,
            pair_weight: 1.0,
            triplet_weight: 0.25,
            temperature: 0.5,
        }
    }
}

/// Generates a special quasirandom structure for a binary disorder.
///
/// # Errors
///
/// Returns [`Error::InvalidInput`] for a structure that does not have one
/// partial-occupancy group with two species, invalid options, or a builder
/// error.
pub fn special_quasirandom(
    structure: &Structure,
    options: &SqsOptions,
) -> Result<SqsResult, Error> {
    if options.supercell.contains(&0) {
        return Err(Error::InvalidInput(
            "SQS supercell repeats must be positive".to_owned(),
        ));
    }
    if options.iterations == 0 || !options.temperature.is_finite() {
        return Err(Error::InvalidInput(
            "SQS iterations and temperature must be valid".to_owned(),
        ));
    }
    let p1 = structure.p1()?;
    let disorder = binary_disorder(&p1)?;
    let [na, nb, nc] = options.supercell;
    let host = host_without_groups(&p1, &disorder.positions)?;
    let host = host.supercell(na, nb, nc)?;
    let mut positions: Vec<[f64; 3]> = Vec::new();
    for cell_a in 0..na {
        for cell_b in 0..nb {
            for cell_c in 0..nc {
                for position in &disorder.positions {
                    positions.push([
                        (position[0] + f64::from(cell_a)) / f64::from(na),
                        (position[1] + f64::from(cell_b)) / f64::from(nb),
                        (position[2] + f64::from(cell_c)) / f64::from(nc),
                    ]);
                }
            }
        }
    }
    if positions.is_empty() {
        return Err(Error::InvalidInput(
            "the SQS has no disordered position".to_owned(),
        ));
    }
    let count_a = ((positions.len() as f64) * disorder.fraction).round() as usize;
    if count_a == 0 || count_a >= positions.len() {
        return Err(Error::InvalidInput(
            "the supercell is too small for both species of the SQS".to_owned(),
        ));
    }
    let shells = build_shells(&positions, host.lattice().matrix(), options.cutoff)?;
    let triples = triangle_triples(&shells);
    let mut rng = SplitMix64::new(options.seed);
    let mut assignment = vec![false; positions.len()];
    for slot in assignment.iter_mut().take(count_a) {
        *slot = true;
    }
    rng.shuffle(&mut assignment);
    let mut objective = evaluate(
        &assignment,
        &shells,
        &triples,
        options.pair_weight,
        options.triplet_weight,
    );
    let mut accepted = 0;
    let mut attempted = 0;
    let mut temperature = options.temperature.max(1e-9);
    for _ in 0..options.iterations {
        let first = (rng.next_u64() % positions.len() as u64) as usize;
        let second = (rng.next_u64() % positions.len() as u64) as usize;
        if assignment.get(first).copied() == assignment.get(second).copied() {
            continue;
        }
        attempted += 1;
        assignment.swap(first, second);
        let candidate = evaluate(
            &assignment,
            &shells,
            &triples,
            options.pair_weight,
            options.triplet_weight,
        );
        let delta = candidate - objective;
        if delta <= 0.0 || rng.uniform() < (-delta / temperature.max(1e-9)).exp() {
            objective = candidate;
            accepted += 1;
        } else {
            assignment.swap(first, second);
        }
        temperature *= 0.9995;
    }
    let mut builder = StructureBuilder::new("sqs", *host.lattice());
    for site in host.sites() {
        builder = builder.site(site.clone());
    }
    let mut labels = Vec::with_capacity(positions.len());
    for (index, (position, &is_a)) in positions.iter().zip(&assignment).enumerate() {
        let species = if is_a {
            &disorder.species[0]
        } else {
            &disorder.species[1]
        };
        labels.push(species.clone());
        builder = builder.site(Site::new(
            format!("{species}{}", index + 1),
            species,
            *position,
        )?);
    }
    let mut reported = Vec::with_capacity(shells.len());
    for shell in &shells {
        reported.push(SqsShell {
            distance: shell.distance,
            pairs: shell.members.len(),
            correlation: pair_correlation(&assignment, &shell.members),
            target: 0.0,
        });
    }
    let (triplet_correlation, _) = triplet_value(&assignment, &triples);
    objective = evaluate(
        &assignment,
        &shells,
        &triples,
        options.pair_weight,
        options.triplet_weight,
    );
    Ok(SqsResult {
        structure: builder.build()?,
        assignment: labels,
        shells: reported,
        triplet_correlation,
        objective,
        accepted,
        attempted,
    })
}

struct BinaryDisorder {
    positions: Vec<[f64; 3]>,
    species: [String; 2],
    fraction: f64,
}

struct Shell {
    distance: f64,
    members: Vec<(usize, usize)>,
}

fn binary_disorder(p1: &xtalforge_core::P1Cell) -> Result<BinaryDisorder, Error> {
    let mut groups: BTreeMap<String, Vec<usize>> = BTreeMap::new();
    for (index, atom) in p1.atoms.iter().enumerate() {
        let key = atom.disorder_group.map_or_else(
            || {
                format!(
                    "{:.5}:{:.5}:{:.5}",
                    atom.frac[0], atom.frac[1], atom.frac[2]
                )
            },
            |group| format!("group:{group}"),
        );
        groups.entry(key).or_default().push(index);
    }
    for members in groups.into_values() {
        if members.len() < 2 {
            continue;
        }
        let mut positions: Vec<[f64; 3]> = Vec::new();
        let mut species: Vec<String> = Vec::new();
        let mut weights: BTreeMap<String, f64> = BTreeMap::new();
        for index in members {
            let Some(atom) = p1.atoms.get(index) else {
                continue;
            };
            if !positions.iter().any(|position| {
                (position[0] - atom.frac[0]).abs() < 1e-5
                    && (position[1] - atom.frac[1]).abs() < 1e-5
                    && (position[2] - atom.frac[2]).abs() < 1e-5
            }) {
                positions.push(atom.frac);
            }
            if !species.contains(&atom.element) {
                species.push(atom.element.clone());
            }
            *weights.entry(atom.element.clone()).or_default() += atom.occupancy;
        }
        if positions.len() != 1 || species.len() != 2 {
            return Err(Error::InvalidInput(
                "SQS generation needs one disordered position with two species".to_owned(),
            ));
        }
        species.sort();
        let first_weight = weights.get(&species[0]).copied().unwrap_or(0.0);
        let second_weight = weights.get(&species[1]).copied().unwrap_or(0.0);
        let total = first_weight + second_weight;
        if total <= 0.0 {
            return Err(Error::InvalidInput(
                "the disordered species have no occupancy".to_owned(),
            ));
        }
        return Ok(BinaryDisorder {
            positions,
            species: [species[0].clone(), species[1].clone()],
            fraction: first_weight / total,
        });
    }
    Err(Error::InvalidInput(
        "the structure has no two-species partial occupancy".to_owned(),
    ))
}

fn host_without_groups(
    p1: &xtalforge_core::P1Cell,
    positions: &[[f64; 3]],
) -> Result<Structure, Error> {
    let mut builder = StructureBuilder::new("host", p1.lattice);
    for atom in &p1.atoms {
        let disordered = positions.iter().any(|position| {
            (position[0] - atom.frac[0]).abs() < 1e-5
                && (position[1] - atom.frac[1]).abs() < 1e-5
                && (position[2] - atom.frac[2]).abs() < 1e-5
        });
        if !disordered {
            builder = builder.site(
                Site::new(atom.label.clone(), &atom.element, atom.frac)?
                    .with_occupancy(atom.occupancy)?,
            );
        }
    }
    builder.build()
}

fn build_shells(
    positions: &[[f64; 3]],
    matrix: [[f64; 3]; 3],
    cutoff: f64,
) -> Result<Vec<Shell>, Error> {
    if !cutoff.is_finite() || cutoff <= 0.0 {
        return Err(Error::InvalidInput(
            "the SQS cutoff must be positive and finite".to_owned(),
        ));
    }
    let mut shells: BTreeMap<u64, Shell> = BTreeMap::new();
    for first in 0..positions.len() {
        for second in (first + 1)..positions.len() {
            let (Some(a), Some(b)) = (positions.get(first), positions.get(second)) else {
                continue;
            };
            let delta = [b[0] - a[0], b[1] - a[1], b[2] - a[2]];
            let mut distance = f64::INFINITY;
            for i in -1..=1 {
                for j in -1..=1 {
                    for k in -1..=1 {
                        let shift = [
                            delta[0] + f64::from(i),
                            delta[1] + f64::from(j),
                            delta[2] + f64::from(k),
                        ];
                        let cart = [
                            shift[0] * matrix[0][0]
                                + shift[1] * matrix[1][0]
                                + shift[2] * matrix[2][0],
                            shift[0] * matrix[0][1]
                                + shift[1] * matrix[1][1]
                                + shift[2] * matrix[2][1],
                            shift[0] * matrix[0][2]
                                + shift[1] * matrix[1][2]
                                + shift[2] * matrix[2][2],
                        ];
                        distance = distance.min(
                            (cart[0] * cart[0] + cart[1] * cart[1] + cart[2] * cart[2]).sqrt(),
                        );
                    }
                }
            }
            if distance > cutoff || distance < 1e-9 {
                continue;
            }
            let bin = (distance * 100.0).round() as u64;
            let entry = shells.entry(bin).or_insert(Shell {
                distance: 0.0,
                members: Vec::new(),
            });
            entry.distance += distance;
            entry.members.push((first, second));
        }
    }
    let mut output = Vec::new();
    for (_, mut shell) in shells {
        shell.distance /= shell.members.len().max(1) as f64;
        output.push(shell);
    }
    Ok(output)
}

fn triangle_triples(shells: &[Shell]) -> Vec<(usize, usize, usize)> {
    let Some(first) = shells.first() else {
        return Vec::new();
    };
    let mut adjacency: BTreeMap<usize, Vec<usize>> = BTreeMap::new();
    for (a, b) in &first.members {
        adjacency.entry(*a).or_default().push(*b);
        adjacency.entry(*b).or_default().push(*a);
    }
    let mut triples = Vec::new();
    for (a, neighbours) in &adjacency {
        for (position, b) in neighbours.iter().enumerate() {
            for c in neighbours.iter().skip(position + 1) {
                if adjacency.get(b).is_some_and(|list| list.contains(c)) {
                    let sorted = [*a, *b, *c];
                    let mut triple = sorted;
                    triple.sort();
                    triples.push((triple[0], triple[1], triple[2]));
                }
            }
        }
    }
    triples.sort();
    triples.dedup();
    triples
}

fn pair_correlation(assignment: &[bool], members: &[(usize, usize)]) -> f64 {
    if members.is_empty() {
        return 0.0;
    }
    let mut sum = 0.0;
    for (first, second) in members {
        let (Some(&a), Some(&b)) = (assignment.get(*first), assignment.get(*second)) else {
            continue;
        };
        sum += if a == b { 1.0 } else { -1.0 };
    }
    sum / members.len() as f64
}

fn triplet_value(assignment: &[bool], triples: &[(usize, usize, usize)]) -> (f64, usize) {
    if triples.is_empty() {
        return (0.0, 0);
    }
    let mut sum = 0.0;
    for (first, second, third) in triples {
        let (Some(&a), Some(&b), Some(&c)) = (
            assignment.get(*first),
            assignment.get(*second),
            assignment.get(*third),
        ) else {
            continue;
        };
        let sign_a = if a { 1.0 } else { -1.0 };
        let sign_b = if b { 1.0 } else { -1.0 };
        let sign_c = if c { 1.0 } else { -1.0 };
        sum += sign_a * sign_b * sign_c;
    }
    (sum / triples.len() as f64, triples.len())
}

fn evaluate(
    assignment: &[bool],
    shells: &[Shell],
    triples: &[(usize, usize, usize)],
    pair_weight: f64,
    triplet_weight: f64,
) -> f64 {
    let mut pair = 0.0;
    for shell in shells {
        let correlation = pair_correlation(assignment, &shell.members);
        pair += correlation * correlation;
    }
    let (triplet, _) = triplet_value(assignment, triples);
    pair_weight * pair + triplet_weight * triplet * triplet
}

struct SplitMix64 {
    state: u64,
}

impl SplitMix64 {
    fn new(seed: u64) -> Self {
        Self { state: seed }
    }

    fn next_u64(&mut self) -> u64 {
        self.state = self.state.wrapping_add(0x9E37_79B9_7F4A_7C15);
        let mut z = self.state;
        z = (z ^ (z >> 30)).wrapping_mul(0xBF58_476D_1CE4_E5B9);
        z = (z ^ (z >> 27)).wrapping_mul(0x94D0_49BB_1331_11EB);
        z ^ (z >> 31)
    }

    fn uniform(&mut self) -> f64 {
        (self.next_u64() >> 11) as f64 / (1u64 << 53) as f64
    }

    fn shuffle<T>(&mut self, values: &mut [T]) {
        for index in (1..values.len()).rev() {
            let other = (self.next_u64() % (index as u64 + 1)) as usize;
            values.swap(index, other);
        }
    }
}
