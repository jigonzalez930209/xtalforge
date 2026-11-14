//! Cluster-expansion features (Task 19.2.4).
//!
//! The module builds the pair, triplet and quadruplet shells of a decorated
//! lattice and evaluates one cluster correlation per shell. It is the feature
//! side of a cluster expansion: no fitting engine is included, the vectors are
//! exported for an external least-squares or machine-learning fit.
//!
//! # Correlation definition
//!
//! For species fractions `c_k`, the centred occupation of site `i` is
//! `y_k(i) = 1 - c_k` when the site carries species `k` and `y_k(i) = -c_k`
//! otherwise. For a cluster `C` of order `n` the correlation is
//!
//! ```text
//! Gamma(C) = sum_k prod_{i in C} y_k(i).
//! ```
//!
//! A random decoration has expectation zero, so all targets are zero. For a
//! binary alloy at `c = 1/2` the pair correlation is the usual Ising
//! correlation scaled by `1/2`; the scaling is stated here so values are never
//! compared with an uncentred implementation by accident.
//!
//! Shells group clusters with the same sorted set of pair distances within
//! `distance_tolerance`. The shell order is the lexicographic order of that
//! quantised key, which makes every exported array independent of the input
//! site order.

use std::collections::BTreeSet;

use serde::{Deserialize, Serialize};

use xtalforge_core::{Error, Lattice};

mod correlation;
mod shells;

pub(crate) use correlation::{
    ClusterCorrelations, correlations_for_assignment, feature_names, species_indices,
};
pub(crate) use shells::{ClusterModel, ClusterShellData, build_cluster_model};

/// Highest cluster order supported: pair, triplet and quadruplet.
pub const MAX_CLUSTER_ORDER: usize = 4;

/// A periodic lattice of decorated sites.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct ClusterLattice {
    /// Periodic cell.
    pub lattice: Lattice,
    /// Fractional site positions, one per decorated site.
    pub positions: Vec<[f64; 3]>,
    /// Species labels. The order defines every species index in the output.
    pub species: Vec<String>,
    /// Species fractions, normalized to sum to one when constructed.
    pub fractions: Vec<f64>,
}

impl ClusterLattice {
    /// Creates a validated cluster lattice and normalizes the fractions.
    ///
    /// # Errors
    ///
    /// Returns [`Error::InvalidInput`] for an empty or duplicate site list,
    /// non-finite positions, repeated or empty species labels, a fraction
    /// count that does not match the species count, negative or non-finite
    /// fractions, or fractions that sum to zero.
    pub fn new(
        lattice: Lattice,
        positions: Vec<[f64; 3]>,
        species: Vec<String>,
        fractions: Vec<f64>,
    ) -> Result<Self, Error> {
        if positions.is_empty() {
            return Err(Error::InvalidInput(
                "a cluster lattice needs at least one position".to_owned(),
            ));
        }
        if !positions.iter().flatten().all(|value| value.is_finite()) {
            return Err(Error::InvalidInput(
                "cluster positions must be finite".to_owned(),
            ));
        }
        if species.is_empty() {
            return Err(Error::InvalidInput(
                "a cluster lattice needs at least one species".to_owned(),
            ));
        }
        let mut labels = BTreeSet::new();
        for label in &species {
            if label.trim().is_empty() || !labels.insert(label.clone()) {
                return Err(Error::InvalidInput(
                    "species labels must be non-empty and unique".to_owned(),
                ));
            }
        }
        if species.len() != fractions.len() {
            return Err(Error::InvalidInput(format!(
                "got {} species but {} fractions",
                species.len(),
                fractions.len()
            )));
        }
        if !fractions
            .iter()
            .all(|fraction| fraction.is_finite() && *fraction >= 0.0)
        {
            return Err(Error::InvalidInput(
                "species fractions must be non-negative and finite".to_owned(),
            ));
        }
        let total: f64 = fractions.iter().sum();
        if total <= 0.0 {
            return Err(Error::InvalidInput(
                "species fractions must sum to a positive value".to_owned(),
            ));
        }
        let mut seen = BTreeSet::new();
        for position in &positions {
            let [first, second, third] = wrapped_key(*position);
            if !seen.insert((first, second, third)) {
                return Err(Error::InvalidInput(
                    "cluster positions must be distinct".to_owned(),
                ));
            }
        }
        let normalized = fractions.iter().map(|fraction| fraction / total).collect();
        Ok(Self {
            lattice,
            positions,
            species,
            fractions: normalized,
        })
    }
}

/// Options for [`cluster_correlation_vectors`].
#[derive(Debug, Clone, PartialEq)]
pub struct ClusterExpansionOptions {
    /// Largest pair distance inside a cluster, in Angstrom.
    pub cutoff: f64,
    /// Highest cluster order, from 2 to [`MAX_CLUSTER_ORDER`].
    pub max_order: usize,
    /// Distance tolerance that groups clusters into one shell, in Angstrom.
    pub distance_tolerance: f64,
    /// Largest total number of clusters before the request is refused.
    pub max_clusters: usize,
}

impl Default for ClusterExpansionOptions {
    fn default() -> Self {
        Self {
            cutoff: 6.0,
            max_order: MAX_CLUSTER_ORDER,
            distance_tolerance: 0.1,
            max_clusters: 500_000,
        }
    }
}

/// One exported shell: order, mean distance and cluster count.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct ClusterShell {
    /// Cluster order: 2 pair, 3 triplet, 4 quadruplet.
    pub order: usize,
    /// Mean pair distance of the clusters in the shell, in Angstrom.
    pub distance: f64,
    /// Number of clusters in the shell.
    pub clusters: usize,
}

/// Correlation vector of one decoration.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct ClusterCorrelationVector {
    /// Species assigned to each site, in the lattice site order.
    pub decoration: Vec<String>,
    /// One correlation per pair shell.
    pub pair: Vec<f64>,
    /// One correlation per triplet shell.
    pub triplet: Vec<f64>,
    /// One correlation per quadruplet shell.
    pub quadruplet: Vec<f64>,
    /// Flattened `pair`, `triplet`, `quadruplet` values for external fitting.
    pub vector: Vec<f64>,
}

/// Cluster correlations of a set of decorations.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct ClusterCorrelationSet {
    /// Species labels in index order.
    pub species: Vec<String>,
    /// Normalized species fractions.
    pub fractions: Vec<f64>,
    /// Pair shells in deterministic order.
    pub pair_shells: Vec<ClusterShell>,
    /// Triplet shells in deterministic order.
    pub triplet_shells: Vec<ClusterShell>,
    /// Quadruplet shells in deterministic order.
    pub quadruplet_shells: Vec<ClusterShell>,
    /// One name per entry of [`ClusterCorrelationVector::vector`].
    pub feature_names: Vec<String>,
    /// One vector per input decoration, in input order.
    pub decorations: Vec<ClusterCorrelationVector>,
}

/// Computes cluster correlation vectors for a set of decorations.
///
/// # Errors
///
/// Returns [`Error::InvalidInput`] for bad options, a decoration with the wrong
/// length or an unknown species, or when the cluster count exceeds
/// [`ClusterExpansionOptions::max_clusters`].
pub fn cluster_correlation_vectors(
    lattice: &ClusterLattice,
    decorations: &[Vec<String>],
    options: &ClusterExpansionOptions,
) -> Result<ClusterCorrelationSet, Error> {
    let model = build_cluster_model(lattice, options)?;
    let mut output = Vec::with_capacity(decorations.len());
    for decoration in decorations {
        let assignment = species_indices(&lattice.species, decoration)?;
        if assignment.len() != lattice.positions.len() {
            return Err(Error::InvalidInput(format!(
                "a decoration needs {} species, got {}",
                lattice.positions.len(),
                assignment.len()
            )));
        }
        let correlations = correlations_for_assignment(&model, &assignment)?;
        let mut vector = Vec::with_capacity(
            correlations.pair.len() + correlations.triplet.len() + correlations.quadruplet.len(),
        );
        vector.extend_from_slice(&correlations.pair);
        vector.extend_from_slice(&correlations.triplet);
        vector.extend_from_slice(&correlations.quadruplet);
        output.push(ClusterCorrelationVector {
            decoration: decoration.clone(),
            pair: correlations.pair,
            triplet: correlations.triplet,
            quadruplet: correlations.quadruplet,
            vector,
        });
    }
    Ok(ClusterCorrelationSet {
        species: model.species.clone(),
        fractions: lattice.fractions.clone(),
        pair_shells: model.pair.iter().map(|shell| shell.to_public(2)).collect(),
        triplet_shells: model
            .triplet
            .iter()
            .map(|shell| shell.to_public(3))
            .collect(),
        quadruplet_shells: model
            .quadruplet
            .iter()
            .map(|shell| shell.to_public(4))
            .collect(),
        feature_names: feature_names(&model),
        decorations: output,
    })
}

fn wrapped_key(position: [f64; 3]) -> [i64; 3] {
    position.map(|value| {
        let wrapped = value - value.floor();
        let wrapped = if wrapped >= 1.0 { 0.0 } else { wrapped };
        (wrapped * 1_000_000.0).round() as i64
    })
}
