//! Integration tests for HNF supercells and derivative-structure enumeration.

#![allow(
    clippy::unwrap_used,
    clippy::expect_used,
    clippy::panic,
    clippy::indexing_slicing
)]

use std::collections::BTreeSet;

use xtalforge_analysis::{DerivativeOptions, derivative_structures, hnf_supercells};
use xtalforge_core::{Lattice, Site, Structure, StructureBuilder};
use xtalforge_symmetry::MatcherOptions;

fn disordered(occupancy_a: f64, occupancy_b: f64) -> Structure {
    StructureBuilder::new("alloy", Lattice::cubic(3.0).unwrap())
        .site(
            Site::new("A1", "Na", [0.0, 0.0, 0.0])
                .unwrap()
                .with_occupancy(occupancy_a)
                .unwrap(),
        )
        .site(
            Site::new("B1", "K", [0.0, 0.0, 0.0])
                .unwrap()
                .with_occupancy(occupancy_b)
                .unwrap(),
        )
        .build()
        .unwrap()
}

fn options() -> DerivativeOptions {
    DerivativeOptions {
        min_index: 8,
        max_index: 8,
        max_structures: 16,
        max_candidates: 10_000,
        concentration_tolerance: 1e-6,
        matrices: Vec::new(),
        matcher: MatcherOptions {
            ltol: 0.05,
            stol: 0.1,
            angle_tol: 1.0,
            allow_supercell: false,
            ..MatcherOptions::default()
        },
    }
}

fn determinant(matrix: [[i32; 3]; 3]) -> i32 {
    let [[a, b, c], [d, e, f], [g, h, i]] = matrix;
    a * (e * i - f * h) - b * (d * i - f * g) + c * (d * h - e * g)
}

fn count(decoration: &[String], species: &str) -> usize {
    decoration.iter().filter(|label| **label == species).count()
}

#[test]
fn hnf_supercells_of_index_eight_are_deterministic_and_unique() {
    let first = hnf_supercells(8).unwrap();
    let second = hnf_supercells(8).unwrap();
    assert_eq!(first, second);
    let unique: BTreeSet<[[i32; 3]; 3]> = first.iter().copied().collect();
    assert_eq!(unique.len(), first.len());
    assert!(first.contains(&[[2, 0, 0], [0, 2, 0], [0, 0, 2]]));
    assert!(first.iter().all(|matrix| determinant(*matrix) == 8));
    assert!(first.iter().all(|matrix| {
        matrix[0][1] < matrix[1][1] && matrix[0][2] < matrix[2][2] && matrix[1][2] < matrix[2][2]
    }));
}

#[test]
fn hnf_index_is_validated() {
    assert!(hnf_supercells(0).is_err());
    assert!(hnf_supercells(129).is_err());
}

#[test]
fn two_by_two_by_two_decorations_keep_the_concentration() {
    let structure = disordered(0.875, 0.125);
    let mut settings = options();
    settings.matrices = vec![[[2, 0, 0], [0, 2, 0], [0, 0, 2]]];
    let first = derivative_structures(&structure, &settings).unwrap();
    let second = derivative_structures(&structure, &settings).unwrap();
    assert_eq!(first, second);
    assert!(!first.is_empty());
    assert!(first.len() <= settings.max_structures);
    for derivative in &first {
        assert_eq!(derivative.matrix, [[2, 0, 0], [0, 2, 0], [0, 0, 2]]);
        assert_eq!(derivative.decoration.len(), 8);
        assert_eq!(derivative.supercell.sites().len(), 8);
        assert_eq!(count(&derivative.decoration, "Na"), 7);
        assert_eq!(count(&derivative.decoration, "K"), 1);
        assert!(derivative.symmetry_count >= 1);
    }
}

#[test]
fn equivalent_supercells_collapse_to_one_derivative() {
    let structure = disordered(0.5, 0.5);
    let mut settings = options();
    settings.matrices = vec![
        [[2, 0, 0], [0, 1, 0], [0, 0, 1]],
        [[1, 0, 0], [0, 2, 0], [0, 0, 1]],
        [[1, 0, 0], [0, 1, 0], [0, 0, 2]],
    ];
    settings.max_structures = 8;
    let derivatives = derivative_structures(&structure, &settings).unwrap();
    assert_eq!(
        derivatives.len(),
        1,
        "{:?}",
        derivatives
            .iter()
            .map(|derivative| derivative.matrix)
            .collect::<Vec<_>>()
    );
    assert_eq!(derivatives[0].decoration.len(), 2);
}

#[test]
fn concentration_limits_skip_a_supercell_that_is_too_small() {
    let structure = disordered(0.8, 0.2);
    let mut settings = options();
    settings.matrices = vec![[[2, 0, 0], [0, 1, 0], [0, 0, 1]]];
    let too_small = derivative_structures(&structure, &settings).unwrap();
    assert!(
        too_small.is_empty(),
        "{:?}",
        too_small
            .iter()
            .map(|derivative| derivative.decoration.clone())
            .collect::<Vec<_>>()
    );
    settings.matrices = vec![[[4, 0, 0], [0, 1, 0], [0, 0, 1]]];
    let derivatives = derivative_structures(&structure, &settings).unwrap();
    assert!(!derivatives.is_empty());
    for derivative in &derivatives {
        assert_eq!(count(&derivative.decoration, "Na"), 3);
        assert_eq!(count(&derivative.decoration, "K"), 1);
    }
}

#[test]
fn a_structure_without_disorder_is_rejected() {
    let ordered = StructureBuilder::new("ordered", Lattice::cubic(3.0).unwrap())
        .site(Site::new("Na1", "Na", [0.0, 0.0, 0.0]).unwrap())
        .build()
        .unwrap();
    assert!(derivative_structures(&ordered, &options()).is_err());
}
