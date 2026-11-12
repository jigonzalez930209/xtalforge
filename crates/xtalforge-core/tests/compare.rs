//! Integration tests for comparison, interpolation and alignment (Task 9.1).

#![allow(
    clippy::unwrap_used,
    clippy::expect_used,
    clippy::panic,
    clippy::indexing_slicing
)]

use xtalforge_core::{
    Lattice, Site, Structure, StructureBuilder, interpolate, kabsch_alignment, periodic_rmsd,
};

fn atom(element: &str, label: &str, frac: [f64; 3]) -> Site {
    Site::new(label, element, frac).unwrap()
}

fn simple() -> Structure {
    StructureBuilder::new("simple", Lattice::cubic(4.0).unwrap())
        .site(atom("C", "C1", [0.0, 0.0, 0.0]))
        .site(atom("C", "C2", [0.5, 0.5, 0.5]))
        .build()
        .unwrap()
}

#[test]
fn kabsch_recovers_a_known_rotation() {
    let points = [
        [1.0, 0.0, 0.0],
        [0.0, 2.0, 0.0],
        [0.0, 0.0, 3.0],
        [1.0, 1.0, 1.0],
    ];
    let rotation = [[0.0, -1.0, 0.0], [1.0, 0.0, 0.0], [0.0, 0.0, 1.0]];
    let translation = [2.0, -1.0, 0.5];
    let target: Vec<[f64; 3]> = points
        .iter()
        .map(|point| {
            [
                rotation[0][0] * point[0]
                    + rotation[0][1] * point[1]
                    + rotation[0][2] * point[2]
                    + translation[0],
                rotation[1][0] * point[0]
                    + rotation[1][1] * point[1]
                    + rotation[1][2] * point[2]
                    + translation[1],
                rotation[2][0] * point[0]
                    + rotation[2][1] * point[1]
                    + rotation[2][2] * point[2]
                    + translation[2],
            ]
        })
        .collect();
    let alignment = kabsch_alignment(&points, &target).unwrap();
    assert!(alignment.rmsd < 1e-9);
    for (actual, expected) in alignment.aligned.iter().zip(&target) {
        for axis in 0..3 {
            assert!((actual[axis] - expected[axis]).abs() < 1e-8);
        }
    }
}

#[test]
fn periodic_rmsd_measures_a_fixed_mapping() {
    let first = simple();
    let shifted = first.origin_shift([0.25, 0.0, 0.0]).unwrap();
    let mapping = [0usize, 1];
    let rms = periodic_rmsd(&first, &shifted, &mapping).unwrap();
    // The shift is 0.25 of a 4 Angstrom cell, so the RMS is 1 Angstrom.
    assert!((rms - 1.0).abs() < 1e-9);
}

#[test]
fn interpolation_returns_the_requested_number_of_frames() {
    let first = simple();
    let second = first.origin_shift([0.1, 0.0, 0.0]).unwrap();
    let frames = interpolate(&first, &second, 4).unwrap();
    assert_eq!(frames.len(), 4);
    let midpoint = &frames[2];
    // The third of four frames is at t = 3/5 of a 0.1 fractional shift.
    assert!((midpoint.sites()[0].frac()[0] - 0.06).abs() < 1e-9);
}

#[test]
fn interpolation_rejects_different_species() {
    let first = simple();
    let second = StructureBuilder::new("other", Lattice::cubic(4.0).unwrap())
        .site(atom("N", "N1", [0.0, 0.0, 0.0]))
        .site(atom("N", "N2", [0.5, 0.5, 0.5]))
        .build()
        .unwrap();
    assert!(interpolate(&first, &second, 2).is_err());
}
