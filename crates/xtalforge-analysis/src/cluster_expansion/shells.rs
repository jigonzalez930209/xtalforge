//! Shell enumeration for the cluster-expansion model.

use std::collections::BTreeMap;

use xtalforge_core::Error;

use super::{ClusterExpansionOptions, ClusterLattice, ClusterShell, MAX_CLUSTER_ORDER};

/// One shell of the internal model.
#[derive(Debug, Clone, PartialEq)]
pub(crate) struct ClusterShellData {
    /// Mean pair distance.
    pub(crate) distance: f64,
    /// Clusters as site-index lists.
    pub(crate) clusters: Vec<Vec<usize>>,
}

impl ClusterShellData {
    pub(crate) fn to_public(&self, order: usize) -> ClusterShell {
        ClusterShell {
            order,
            distance: self.distance,
            clusters: self.clusters.len(),
        }
    }
}

/// Internal cluster model shared with the multicomponent SQS search.
#[derive(Debug, Clone, PartialEq)]
pub(crate) struct ClusterModel {
    /// Site count.
    pub(crate) n_sites: usize,
    /// Species labels.
    pub(crate) species: Vec<String>,
    /// Normalized fractions.
    pub(crate) fractions: Vec<f64>,
    /// Pair shells.
    pub(crate) pair: Vec<ClusterShellData>,
    /// Triplet shells.
    pub(crate) triplet: Vec<ClusterShellData>,
    /// Quadruplet shells.
    pub(crate) quadruplet: Vec<ClusterShellData>,
}

/// Validates options and builds the internal cluster model.
pub(crate) fn build_cluster_model(
    lattice: &ClusterLattice,
    options: &ClusterExpansionOptions,
) -> Result<ClusterModel, Error> {
    validate_options(options)?;
    let mut total = 0usize;
    let pair = build_pairs(lattice, options, &mut total)?;
    let triplet = if options.max_order >= 3 {
        build_triplets(lattice, options, &mut total)?
    } else {
        Vec::new()
    };
    let quadruplet = if options.max_order >= MAX_CLUSTER_ORDER {
        build_quadruplets(lattice, options, &mut total)?
    } else {
        Vec::new()
    };
    Ok(ClusterModel {
        n_sites: lattice.positions.len(),
        species: lattice.species.clone(),
        fractions: lattice.fractions.clone(),
        pair,
        triplet,
        quadruplet,
    })
}

fn validate_options(options: &ClusterExpansionOptions) -> Result<(), Error> {
    if !options.cutoff.is_finite() || options.cutoff <= 0.0 {
        return Err(Error::InvalidInput(format!(
            "cluster cutoff must be positive and finite, got {}",
            options.cutoff
        )));
    }
    if !(2..=MAX_CLUSTER_ORDER).contains(&options.max_order) {
        return Err(Error::InvalidInput(format!(
            "cluster max_order must be between 2 and {MAX_CLUSTER_ORDER}, got {}",
            options.max_order
        )));
    }
    if !options.distance_tolerance.is_finite() || options.distance_tolerance <= 0.0 {
        return Err(Error::InvalidInput(format!(
            "cluster distance tolerance must be positive and finite, got {}",
            options.distance_tolerance
        )));
    }
    if options.max_clusters == 0 {
        return Err(Error::InvalidInput(
            "cluster max_clusters must be positive".to_owned(),
        ));
    }
    Ok(())
}

fn build_pairs(
    lattice: &ClusterLattice,
    options: &ClusterExpansionOptions,
    total: &mut usize,
) -> Result<Vec<ClusterShellData>, Error> {
    let mut shells: BTreeMap<Vec<i64>, ShellAccumulator> = BTreeMap::new();
    for (first_index, first) in lattice.positions.iter().enumerate() {
        for (second_index, second) in lattice.positions.iter().enumerate().skip(first_index + 1) {
            let distance = lattice.lattice.distance(*first, *second);
            if distance > options.cutoff {
                continue;
            }
            add_cluster(
                &mut shells,
                distance_key(&[distance], options.distance_tolerance),
                distance,
                vec![first_index, second_index],
                total,
                options.max_clusters,
            )?;
        }
    }
    Ok(finish_shells(shells))
}

fn build_triplets(
    lattice: &ClusterLattice,
    options: &ClusterExpansionOptions,
    total: &mut usize,
) -> Result<Vec<ClusterShellData>, Error> {
    let mut shells: BTreeMap<Vec<i64>, ShellAccumulator> = BTreeMap::new();
    for (first_index, first) in lattice.positions.iter().enumerate() {
        for (second_index, second) in lattice.positions.iter().enumerate().skip(first_index + 1) {
            let first_second = lattice.lattice.distance(*first, *second);
            if first_second > options.cutoff {
                continue;
            }
            for (third_index, third) in lattice.positions.iter().enumerate().skip(second_index + 1)
            {
                let first_third = lattice.lattice.distance(*first, *third);
                let second_third = lattice.lattice.distance(*second, *third);
                if first_third > options.cutoff || second_third > options.cutoff {
                    continue;
                }
                let distances = [first_second, first_third, second_third];
                let mean = distances.iter().sum::<f64>() / 3.0;
                add_cluster(
                    &mut shells,
                    distance_key(&distances, options.distance_tolerance),
                    mean,
                    vec![first_index, second_index, third_index],
                    total,
                    options.max_clusters,
                )?;
            }
        }
    }
    Ok(finish_shells(shells))
}

fn build_quadruplets(
    lattice: &ClusterLattice,
    options: &ClusterExpansionOptions,
    total: &mut usize,
) -> Result<Vec<ClusterShellData>, Error> {
    let mut shells: BTreeMap<Vec<i64>, ShellAccumulator> = BTreeMap::new();
    for (first_index, first) in lattice.positions.iter().enumerate() {
        for (second_index, second) in lattice.positions.iter().enumerate().skip(first_index + 1) {
            let first_second = lattice.lattice.distance(*first, *second);
            if first_second > options.cutoff {
                continue;
            }
            for (third_index, third) in lattice.positions.iter().enumerate().skip(second_index + 1)
            {
                let first_third = lattice.lattice.distance(*first, *third);
                let second_third = lattice.lattice.distance(*second, *third);
                if first_third > options.cutoff || second_third > options.cutoff {
                    continue;
                }
                for (fourth_index, fourth) in
                    lattice.positions.iter().enumerate().skip(third_index + 1)
                {
                    let first_fourth = lattice.lattice.distance(*first, *fourth);
                    let second_fourth = lattice.lattice.distance(*second, *fourth);
                    let third_fourth = lattice.lattice.distance(*third, *fourth);
                    if first_fourth > options.cutoff
                        || second_fourth > options.cutoff
                        || third_fourth > options.cutoff
                    {
                        continue;
                    }
                    let distances = [
                        first_second,
                        first_third,
                        first_fourth,
                        second_third,
                        second_fourth,
                        third_fourth,
                    ];
                    let mean = distances.iter().sum::<f64>() / 6.0;
                    add_cluster(
                        &mut shells,
                        distance_key(&distances, options.distance_tolerance),
                        mean,
                        vec![first_index, second_index, third_index, fourth_index],
                        total,
                        options.max_clusters,
                    )?;
                }
            }
        }
    }
    Ok(finish_shells(shells))
}

#[derive(Debug, Clone, Default)]
struct ShellAccumulator {
    distance_sum: f64,
    clusters: Vec<Vec<usize>>,
}

fn add_cluster(
    shells: &mut BTreeMap<Vec<i64>, ShellAccumulator>,
    key: Vec<i64>,
    distance: f64,
    cluster: Vec<usize>,
    total: &mut usize,
    max_clusters: usize,
) -> Result<(), Error> {
    if *total >= max_clusters {
        return Err(Error::InvalidInput(format!(
            "cluster expansion exceeds {max_clusters} clusters; raise max_clusters or lower the cutoff"
        )));
    }
    let entry = shells.entry(key).or_default();
    entry.distance_sum += distance;
    entry.clusters.push(cluster);
    *total += 1;
    Ok(())
}

fn finish_shells(shells: BTreeMap<Vec<i64>, ShellAccumulator>) -> Vec<ClusterShellData> {
    shells
        .into_values()
        .map(|shell| ClusterShellData {
            distance: shell.distance_sum / (shell.clusters.len().max(1) as f64),
            clusters: shell.clusters,
        })
        .collect()
}

fn distance_key(distances: &[f64], tolerance: f64) -> Vec<i64> {
    let mut key: Vec<i64> = distances
        .iter()
        .map(|distance| (distance / tolerance).round() as i64)
        .collect();
    key.sort_unstable();
    key
}
