//! Integration tests for multicomponent SQS and the binary backward
//! compatibility of the original entry point.

#![allow(
    clippy::unwrap_used,
    clippy::expect_used,
    clippy::panic,
    clippy::indexing_slicing
)]

use xtalforge_analysis::{
    MulticomponentSqsOptions, SqsOptions, multicomponent_sqs, special_quasirandom,
};
use xtalforge_core::{Lattice, Site, Structure, StructureBuilder};

fn binary_alloy() -> Structure {
    StructureBuilder::new("binary", Lattice::cubic(3.0).unwrap())
        .site(
            Site::new("Na1", "Na", [0.0, 0.0, 0.0])
                .unwrap()
                .with_occupancy(0.5)
                .unwrap(),
        )
        .site(
            Site::new("K1", "K", [0.0, 0.0, 0.0])
                .unwrap()
                .with_occupancy(0.5)
                .unwrap(),
        )
        .build()
        .unwrap()
}

fn ternary_alloy() -> Structure {
    StructureBuilder::new("ternary", Lattice::cubic(3.0).unwrap())
        .site(
            Site::new("Na1", "Na", [0.0, 0.0, 0.0])
                .unwrap()
                .with_occupancy(1.0 / 3.0)
                .unwrap(),
        )
        .site(
            Site::new("K1", "K", [0.0, 0.0, 0.0])
                .unwrap()
                .with_occupancy(1.0 / 3.0)
                .unwrap(),
        )
        .site(
            Site::new("Li1", "Li", [0.0, 0.0, 0.0])
                .unwrap()
                .with_occupancy(1.0 / 3.0)
                .unwrap(),
        )
        .build()
        .unwrap()
}

fn ternary_options() -> MulticomponentSqsOptions {
    MulticomponentSqsOptions {
        supercell: [3, 1, 1],
        seed: 42,
        cutoff: 3.1,
        iterations: 500,
        max_order: 4,
        distance_tolerance: 0.1,
        max_clusters: 10_000,
        concentration_tolerance: 1e-6,
        pair_weight: 1.0,
        triplet_weight: 0.25,
        quadruplet_weight: 0.1,
        pair_targets: Vec::new(),
        triplet_targets: Vec::new(),
        quadruplet_targets: Vec::new(),
    }
}

#[test]
fn ternary_sqs_is_reproducible_and_preserves_counts() {
    let first = multicomponent_sqs(&ternary_alloy(), &ternary_options()).unwrap();
    let second = multicomponent_sqs(&ternary_alloy(), &ternary_options()).unwrap();
    assert_eq!(first, second);
    assert_eq!(first.assignment.len(), 3);
    assert_eq!(first.species.len(), 3);
    assert_eq!(first.fractions.len(), 3);
    assert!(!first.pair_shells.is_empty());
    assert!(!first.triplet_shells.is_empty());
    assert!(first.objective.is_finite());
    for label in ["Na", "K", "Li"] {
        assert_eq!(
            first
                .assignment
                .iter()
                .filter(|species| **species == label)
                .count(),
            1
        );
    }
}

#[test]
fn binary_entry_point_is_unchanged() {
    let options = SqsOptions {
        supercell: [2, 2, 2],
        seed: 11,
        cutoff: 5.0,
        iterations: 500,
        ..SqsOptions::default()
    };
    let binary = special_quasirandom(&binary_alloy(), &options).unwrap();
    assert_eq!(binary.assignment.len(), 8);
    assert!(!binary.shells.is_empty());
    assert!(binary.objective.is_finite());
    let multicomponent = multicomponent_sqs(
        &binary_alloy(),
        &MulticomponentSqsOptions {
            supercell: [2, 2, 2],
            seed: 11,
            cutoff: 5.0,
            iterations: 500,
            ..MulticomponentSqsOptions::default()
        },
    )
    .unwrap();
    assert_eq!(multicomponent.assignment.len(), 8);
    assert!(!multicomponent.pair_shells.is_empty());
}

#[test]
fn quadruplet_shells_are_reported() {
    let options = MulticomponentSqsOptions {
        supercell: [2, 2, 2],
        seed: 3,
        cutoff: 4.5,
        iterations: 200,
        max_order: 4,
        max_clusters: 10_000,
        ..MulticomponentSqsOptions::default()
    };
    let result = multicomponent_sqs(&binary_alloy(), &options).unwrap();
    assert!(!result.quadruplet_shells.is_empty());
    assert!(
        result
            .quadruplet_shells
            .iter()
            .all(|shell| shell.correlation.is_finite())
    );
}

#[test]
fn target_lengths_and_option_ranges_are_validated() {
    let mut options = ternary_options();
    options.pair_targets = vec![1.0, 2.0];
    assert!(multicomponent_sqs(&ternary_alloy(), &options).is_err());
    let mut options = ternary_options();
    options.iterations = 0;
    assert!(multicomponent_sqs(&ternary_alloy(), &options).is_err());
    let mut options = ternary_options();
    options.supercell = [2, 1, 1];
    assert!(multicomponent_sqs(&ternary_alloy(), &options).is_err());
}

#[test]
fn an_ordered_structure_is_rejected() {
    let ordered = StructureBuilder::new("ordered", Lattice::cubic(3.0).unwrap())
        .site(Site::new("Na1", "Na", [0.0, 0.0, 0.0]).unwrap())
        .build()
        .unwrap();
    assert!(multicomponent_sqs(&ordered, &ternary_options()).is_err());
}
