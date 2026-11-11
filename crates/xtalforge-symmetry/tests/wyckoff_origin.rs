//! Wyckoff letters follow the input origin, as in spglib and the CIF
//! (functional audit 2026-11-11).

#![allow(
    clippy::unwrap_used,
    clippy::expect_used,
    clippy::panic,
    clippy::indexing_slicing
)]

use xtalforge_core::{Lattice, Site, StructureBuilder};
use xtalforge_symmetry::{SymmetryOptions, detect_symmetry};

fn silicon(shift: f64) -> xtalforge_core::Structure {
    let points = [
        [0.0, 0.0, 0.0],
        [0.0, 0.5, 0.5],
        [0.5, 0.0, 0.5],
        [0.5, 0.5, 0.0],
        [0.25, 0.25, 0.25],
        [0.25, 0.75, 0.75],
        [0.75, 0.25, 0.75],
        [0.75, 0.75, 0.25],
    ];
    let mut builder = StructureBuilder::new("si", Lattice::cubic(5.431).unwrap());
    for (index, point) in points.iter().enumerate() {
        let moved = point.map(|value| (value + shift).rem_euclid(1.0));
        builder = builder.site(Site::new(format!("Si{index}"), "Si", moved).unwrap());
    }
    builder.build().unwrap()
}

#[test]
fn silicon_at_the_origin_is_8a_and_shifted_silicon_is_8b() {
    let at_origin = detect_symmetry(&silicon(0.0), SymmetryOptions::default()).unwrap();
    assert_eq!(at_origin.number, 227);
    assert!(
        at_origin.wyckoffs.iter().all(|letter| *letter == 'a'),
        "{:?}",
        at_origin.wyckoffs
    );
    let shifted = detect_symmetry(&silicon(0.5), SymmetryOptions::default()).unwrap();
    assert!(
        shifted.wyckoffs.iter().all(|letter| *letter == 'b'),
        "{:?}",
        shifted.wyckoffs
    );
}
