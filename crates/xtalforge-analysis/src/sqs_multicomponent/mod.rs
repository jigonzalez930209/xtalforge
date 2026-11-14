//! Multicomponent special quasirandom structures (Task 19.2.3).
//!
//! The binary [`special_quasirandom`](crate::special_quasirandom) entry point
//! stays unchanged. This module adds a general entry point for two or more
//! species with explicit pair, triplet and quadruplet correlation targets.
//!
//! # Search
//!
//! The initial decoration is a random shuffle of the exact species counts,
//! seeded by [`MulticomponentSqsOptions::seed`]. A greedy local search then
//! swaps pairs of unlike species and accepts a swap only when the weighted
//! squared deviation from the targets decreases. The whole procedure is
//! deterministic for a fixed seed and fixed input.
//!
//! # Correlation definition
//!
//! Correlations use the centred occupation variables of
//! [`crate::cluster_expansion`]: the random target is zero for every shell and
//! a binary `1/2` alloy is proportional to the Ising correlation by `1/2`.

use serde::{Deserialize, Serialize};

use xtalforge_core::{Error, Site, Structure, StructureBuilder};

use crate::cluster_expansion::{
    ClusterCorrelations, ClusterExpansionOptions, ClusterLattice, ClusterShellData,
    build_cluster_model, correlations_for_assignment,
};
use crate::disorder::disorder_groups;

mod support;

use support::{SplitMix64, allocate_counts, host_without_position};

/// One reported shell of a multicomponent SQS.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct SqsCorrelation {
    /// Mean pair distance of the shell, in Angstrom.
    pub distance: f64,
    /// Number of clusters in the shell.
    pub clusters: usize,
    /// Achieved centred correlation.
    pub correlation: f64,
    /// Requested target correlation.
    pub target: f64,
}

/// Result of [`multicomponent_sqs`].
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct MulticomponentSqsResult {
    /// Ordered SQS supercell.
    pub structure: Structure,
    /// Species assigned to each disordered position, in cell order.
    pub assignment: Vec<String>,
    /// Species labels in index order.
    pub species: Vec<String>,
    /// Normalized species fractions.
    pub fractions: Vec<f64>,
    /// Pair shells with achieved and target correlations.
    pub pair_shells: Vec<SqsCorrelation>,
    /// Triplet shells with achieved and target correlations.
    pub triplet_shells: Vec<SqsCorrelation>,
    /// Quadruplet shells with achieved and target correlations.
    pub quadruplet_shells: Vec<SqsCorrelation>,
    /// Weighted squared deviation from every target; zero is a perfect match.
    pub objective: f64,
    /// Number of accepted swaps.
    pub accepted: usize,
    /// Number of attempted swaps between unlike species.
    pub attempted: usize,
}

/// Options for [`multicomponent_sqs`].
#[derive(Debug, Clone, PartialEq)]
pub struct MulticomponentSqsOptions {
    /// Diagonal supercell repeats.
    pub supercell: [u32; 3],
    /// Explicit seed of the deterministic shuffle and swap sequence.
    pub seed: u64,
    /// Largest pair distance inside a cluster, in Angstrom.
    pub cutoff: f64,
    /// Number of greedy swap attempts.
    pub iterations: usize,
    /// Highest cluster order, from 2 to 4.
    pub max_order: usize,
    /// Distance tolerance that groups clusters into one shell, in Angstrom.
    pub distance_tolerance: f64,
    /// Largest total number of clusters accepted by the model.
    pub max_clusters: usize,
    /// Occupancy below which a species is treated as absent.
    pub concentration_tolerance: f64,
    /// Weight of the pair objective.
    pub pair_weight: f64,
    /// Weight of the triplet objective.
    pub triplet_weight: f64,
    /// Weight of the quadruplet objective.
    pub quadruplet_weight: f64,
    /// Pair targets in shell order; empty means all zeros.
    pub pair_targets: Vec<f64>,
    /// Triplet targets in shell order; empty means all zeros.
    pub triplet_targets: Vec<f64>,
    /// Quadruplet targets in shell order; empty means all zeros.
    pub quadruplet_targets: Vec<f64>,
}

impl Default for MulticomponentSqsOptions {
    fn default() -> Self {
        Self {
            supercell: [2, 2, 2],
            seed: 0,
            cutoff: 6.0,
            iterations: 20_000,
            max_order: 4,
            distance_tolerance: 0.1,
            max_clusters: 200_000,
            concentration_tolerance: 1e-6,
            pair_weight: 1.0,
            triplet_weight: 0.25,
            quadruplet_weight: 0.1,
            pair_targets: Vec::new(),
            triplet_targets: Vec::new(),
            quadruplet_targets: Vec::new(),
        }
    }
}

/// Builds a multicomponent special quasirandom structure.
///
/// The structure must contain exactly one disorder group with one fractional
/// position and at least two species with positive occupancy. The group is
/// expanded with the diagonal `supercell`; the exact species counts follow the
/// largest-remainder rule.
///
/// # Errors
///
/// Returns [`Error::InvalidInput`] for bad options, a structure without exactly
/// one single-position disorder group, fewer than two species, a supercell too
/// small for every species, a target vector whose length does not match its
/// shell count, or the errors of the supercell construction and cluster model.
pub fn multicomponent_sqs(
    structure: &Structure,
    options: &MulticomponentSqsOptions,
) -> Result<MulticomponentSqsResult, Error> {
    validate_options(options)?;
    let groups = disorder_groups(structure)?;
    if groups.len() != 1 {
        return Err(Error::InvalidInput(
            "multicomponent SQS needs exactly one disorder group".to_owned(),
        ));
    }
    let Some(group) = groups.first() else {
        return Err(Error::InvalidInput(
            "the structure has no disorder group".to_owned(),
        ));
    };
    if group.positions.len() != 1 {
        return Err(Error::InvalidInput(
            "multicomponent SQS needs one position per disorder group".to_owned(),
        ));
    }
    let occupants: Vec<_> = group
        .aggregate_occupants()
        .into_iter()
        .filter(|occupant| occupant.occupancy > options.concentration_tolerance)
        .collect();
    if occupants.len() < 2 {
        return Err(Error::InvalidInput(
            "multicomponent SQS needs at least two species with occupancy".to_owned(),
        ));
    }
    let total: f64 = occupants.iter().map(|occupant| occupant.occupancy).sum();
    if total <= 0.0 {
        return Err(Error::InvalidInput(
            "the disorder group has no occupancy".to_owned(),
        ));
    }
    let species: Vec<String> = occupants
        .iter()
        .map(|occupant| occupant.element.clone())
        .collect();
    let fractions: Vec<f64> = occupants
        .iter()
        .map(|occupant| occupant.occupancy / total)
        .collect();
    let Some(position) = group.positions.first() else {
        return Err(Error::InvalidInput(
            "the disorder group has no position".to_owned(),
        ));
    };
    let p1 = structure.as_p1()?;
    let host = host_without_position(&p1, position.frac)?;
    let [na, nb, nc] = options.supercell;
    let host = host.supercell(na, nb, nc)?;
    let mut positions = Vec::with_capacity((na * nb * nc) as usize);
    for cell_a in 0..na {
        for cell_b in 0..nb {
            for cell_c in 0..nc {
                positions.push([
                    (position.frac[0] + f64::from(cell_a)) / f64::from(na),
                    (position.frac[1] + f64::from(cell_b)) / f64::from(nb),
                    (position.frac[2] + f64::from(cell_c)) / f64::from(nc),
                ]);
            }
        }
    }
    let cluster_lattice = ClusterLattice::new(
        *host.lattice(),
        positions.clone(),
        species.clone(),
        fractions.clone(),
    )?;
    let cluster_options = ClusterExpansionOptions {
        cutoff: options.cutoff,
        max_order: options.max_order,
        distance_tolerance: options.distance_tolerance,
        max_clusters: options.max_clusters,
    };
    let model = build_cluster_model(&cluster_lattice, &cluster_options)?;
    let Some(counts) = allocate_counts(&fractions, positions.len(), options) else {
        return Err(Error::InvalidInput(
            "the supercell is too small for every SQS species".to_owned(),
        ));
    };
    let pair_targets = targets(&options.pair_targets, model.pair.len(), "pair")?;
    let triplet_targets = targets(&options.triplet_targets, model.triplet.len(), "triplet")?;
    let quadruplet_targets = targets(
        &options.quadruplet_targets,
        model.quadruplet.len(),
        "quadruplet",
    )?;
    let mut assignment = Vec::with_capacity(positions.len());
    for (species_index, count) in counts.iter().enumerate() {
        for _ in 0..*count {
            assignment.push(species_index);
        }
    }
    let mut rng = SplitMix64::new(options.seed);
    rng.shuffle(&mut assignment);
    let mut correlations = correlations_for_assignment(&model, &assignment)?;
    let mut objective = evaluate(
        &correlations,
        &pair_targets,
        &triplet_targets,
        &quadruplet_targets,
        options,
    );
    let mut accepted = 0;
    let mut attempted = 0;
    for _ in 0..options.iterations {
        if assignment.len() < 2 {
            break;
        }
        let first = (rng.next_u64() % assignment.len() as u64) as usize;
        let second = (rng.next_u64() % assignment.len() as u64) as usize;
        if assignment.get(first).copied() == assignment.get(second).copied() {
            continue;
        }
        attempted += 1;
        assignment.swap(first, second);
        let candidate_correlations = correlations_for_assignment(&model, &assignment)?;
        let candidate = evaluate(
            &candidate_correlations,
            &pair_targets,
            &triplet_targets,
            &quadruplet_targets,
            options,
        );
        if candidate + 1e-12 < objective {
            correlations = candidate_correlations;
            objective = candidate;
            accepted += 1;
        } else {
            assignment.swap(first, second);
        }
    }
    let pair_shells = report_shells(&model.pair, &correlations.pair, &pair_targets)?;
    let triplet_shells = report_shells(&model.triplet, &correlations.triplet, &triplet_targets)?;
    let quadruplet_shells = report_shells(
        &model.quadruplet,
        &correlations.quadruplet,
        &quadruplet_targets,
    )?;
    let mut labels = Vec::with_capacity(assignment.len());
    for index in &assignment {
        let Some(label) = species.get(*index) else {
            return Err(Error::InvalidInput(
                "the SQS assignment references an unknown species".to_owned(),
            ));
        };
        labels.push(label.clone());
    }
    let mut builder = StructureBuilder::new("sqs", *host.lattice());
    for site in host.sites() {
        builder = builder.site(site.clone());
    }
    for (index, (frac, label)) in positions.iter().zip(&labels).enumerate() {
        builder = builder.site(Site::new(format!("{label}{}", index + 1), label, *frac)?);
    }
    Ok(MulticomponentSqsResult {
        structure: builder.build()?,
        assignment: labels,
        species,
        fractions,
        pair_shells,
        triplet_shells,
        quadruplet_shells,
        objective,
        accepted,
        attempted,
    })
}

fn validate_options(options: &MulticomponentSqsOptions) -> Result<(), Error> {
    if options.supercell.contains(&0) {
        return Err(Error::InvalidInput(
            "SQS supercell repeats must be positive".to_owned(),
        ));
    }
    if options.iterations == 0 {
        return Err(Error::InvalidInput(
            "SQS iterations must be positive".to_owned(),
        ));
    }
    if !options.cutoff.is_finite() || options.cutoff <= 0.0 {
        return Err(Error::InvalidInput(
            "SQS cutoff must be positive and finite".to_owned(),
        ));
    }
    if !(2..=4).contains(&options.max_order) {
        return Err(Error::InvalidInput(
            "SQS max_order must be between 2 and 4".to_owned(),
        ));
    }
    if !options.distance_tolerance.is_finite() || options.distance_tolerance <= 0.0 {
        return Err(Error::InvalidInput(
            "SQS distance tolerance must be positive and finite".to_owned(),
        ));
    }
    if options.max_clusters == 0 {
        return Err(Error::InvalidInput(
            "SQS max_clusters must be positive".to_owned(),
        ));
    }
    if !options.concentration_tolerance.is_finite() || options.concentration_tolerance < 0.0 {
        return Err(Error::InvalidInput(
            "SQS concentration tolerance must be non-negative and finite".to_owned(),
        ));
    }
    for weight in [
        options.pair_weight,
        options.triplet_weight,
        options.quadruplet_weight,
    ] {
        if !weight.is_finite() || weight < 0.0 {
            return Err(Error::InvalidInput(
                "SQS weights must be non-negative and finite".to_owned(),
            ));
        }
    }
    Ok(())
}

fn targets(values: &[f64], count: usize, name: &str) -> Result<Vec<f64>, Error> {
    if values.is_empty() {
        return Ok(vec![0.0; count]);
    }
    if values.len() != count {
        return Err(Error::InvalidInput(format!(
            "{name} target count {} does not match {count} shells",
            values.len()
        )));
    }
    if !values.iter().all(|value| value.is_finite()) {
        return Err(Error::InvalidInput(format!(
            "{name} targets must be finite"
        )));
    }
    Ok(values.to_vec())
}

fn evaluate(
    correlations: &ClusterCorrelations,
    pair_targets: &[f64],
    triplet_targets: &[f64],
    quadruplet_targets: &[f64],
    options: &MulticomponentSqsOptions,
) -> f64 {
    squared_deviation(&correlations.pair, pair_targets) * options.pair_weight
        + squared_deviation(&correlations.triplet, triplet_targets) * options.triplet_weight
        + squared_deviation(&correlations.quadruplet, quadruplet_targets)
            * options.quadruplet_weight
}

fn squared_deviation(values: &[f64], targets: &[f64]) -> f64 {
    values
        .iter()
        .zip(targets)
        .map(|(value, target)| {
            let difference = value - target;
            difference * difference
        })
        .sum()
}

fn report_shells(
    shells: &[ClusterShellData],
    correlations: &[f64],
    targets: &[f64],
) -> Result<Vec<SqsCorrelation>, Error> {
    if shells.len() != correlations.len() || shells.len() != targets.len() {
        return Err(Error::InvalidInput(
            "SQS shell and correlation counts disagree".to_owned(),
        ));
    }
    Ok(shells
        .iter()
        .zip(correlations)
        .zip(targets)
        .map(|((shell, correlation), target)| SqsCorrelation {
            distance: shell.distance,
            clusters: shell.clusters.len(),
            correlation: *correlation,
            target: *target,
        })
        .collect())
}
