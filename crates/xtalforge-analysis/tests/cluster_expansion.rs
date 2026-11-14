//! Integration tests for cluster-correlation feature vectors.

#![allow(
    clippy::unwrap_used,
    clippy::expect_used,
    clippy::panic,
    clippy::indexing_slicing
)]

use xtalforge_analysis::{ClusterExpansionOptions, ClusterLattice, cluster_correlation_vectors};
use xtalforge_core::Lattice;

#[test]
fn two_site_pair_correlations_are_exact() {
    let lattice = ClusterLattice::new(
        Lattice::cubic(4.0).unwrap(),
        vec![[0.0, 0.0, 0.0], [0.5, 0.0, 0.0]],
        vec!["A".to_owned(), "B".to_owned()],
        vec![0.5, 0.5],
    )
    .unwrap();
    let options = ClusterExpansionOptions {
        cutoff: 2.5,
        max_order: 4,
        distance_tolerance: 0.1,
        max_clusters: 100,
    };
    let decorations = vec![
        vec!["A".to_owned(), "A".to_owned()],
        vec!["A".to_owned(), "B".to_owned()],
        vec!["B".to_owned(), "A".to_owned()],
        vec!["B".to_owned(), "B".to_owned()],
    ];
    let result = cluster_correlation_vectors(&lattice, &decorations, &options).unwrap();
    assert_eq!(result.pair_shells.len(), 1);
    assert_eq!(result.pair_shells[0].order, 2);
    assert_eq!(result.pair_shells[0].clusters, 1);
    assert!(result.triplet_shells.is_empty());
    assert!(result.quadruplet_shells.is_empty());
    assert_eq!(result.feature_names, vec!["pair[0]".to_owned()]);
    let pair: Vec<f64> = result
        .decorations
        .iter()
        .map(|vector| vector.pair[0])
        .collect();
    assert_eq!(pair, vec![0.5, -0.5, -0.5, 0.5]);
    assert_eq!(result.decorations[0].vector, vec![0.5]);
}

#[test]
fn small_cubic_lattice_exports_pair_triplet_and_quadruplet_shells() {
    let mut positions = Vec::new();
    for first in 0..2 {
        for second in 0..2 {
            for third in 0..2 {
                positions.push([
                    f64::from(first) * 0.5,
                    f64::from(second) * 0.5,
                    f64::from(third) * 0.5,
                ]);
            }
        }
    }
    let lattice = ClusterLattice::new(
        Lattice::cubic(2.0).unwrap(),
        positions,
        vec!["A".to_owned(), "B".to_owned(), "C".to_owned()],
        vec![1.0, 1.0, 1.0],
    )
    .unwrap();
    assert!((lattice.fractions[0] - 1.0 / 3.0).abs() < 1e-12);
    let options = ClusterExpansionOptions {
        cutoff: 1.5,
        max_order: 4,
        distance_tolerance: 0.1,
        max_clusters: 10_000,
    };
    let first = vec![
        "A".to_owned(),
        "B".to_owned(),
        "C".to_owned(),
        "A".to_owned(),
        "B".to_owned(),
        "C".to_owned(),
        "A".to_owned(),
        "B".to_owned(),
    ];
    let second: Vec<String> = first.iter().rev().cloned().collect();
    let result = cluster_correlation_vectors(&lattice, &[first.clone(), second], &options).unwrap();
    assert!(!result.pair_shells.is_empty());
    assert!(!result.triplet_shells.is_empty());
    assert!(!result.quadruplet_shells.is_empty());
    let expected =
        result.pair_shells.len() + result.triplet_shells.len() + result.quadruplet_shells.len();
    assert_eq!(result.feature_names.len(), expected);
    assert_eq!(result.decorations[0].vector.len(), expected);
    assert!(
        result.decorations[0]
            .vector
            .iter()
            .all(|value| value.is_finite())
    );
    let repeated = cluster_correlation_vectors(&lattice, &[first], &options).unwrap();
    assert_eq!(repeated.decorations[0].vector, result.decorations[0].vector);
}

#[test]
fn invalid_cluster_inputs_are_rejected() {
    let lattice = Lattice::cubic(2.0).unwrap();
    assert!(
        ClusterLattice::new(
            lattice,
            vec![[0.0, 0.0, 0.0], [0.0, 0.0, 0.0]],
            vec!["A".to_owned(), "B".to_owned()],
            vec![0.5, 0.5],
        )
        .is_err()
    );
    assert!(
        ClusterLattice::new(
            lattice,
            vec![[0.0, 0.0, 0.0]],
            vec!["A".to_owned(), "B".to_owned()],
            vec![1.0, -1.0],
        )
        .is_err()
    );
    let valid = ClusterLattice::new(
        lattice,
        vec![[0.0, 0.0, 0.0], [0.5, 0.0, 0.0]],
        vec!["A".to_owned(), "B".to_owned()],
        vec![0.5, 0.5],
    )
    .unwrap();
    let options = ClusterExpansionOptions::default();
    assert!(
        cluster_correlation_vectors(&valid, &[vec!["A".to_owned(), "C".to_owned()]], &options,)
            .is_err()
    );
    assert!(cluster_correlation_vectors(&valid, &[vec!["A".to_owned()]], &options).is_err());
}
