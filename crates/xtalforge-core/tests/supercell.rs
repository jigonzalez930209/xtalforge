//! Integration tests for supercells and origin shifts (Task 1.7).

#![allow(
    clippy::unwrap_used,
    clippy::expect_used,
    clippy::panic,
    clippy::indexing_slicing
)]

use serde_json::json;
use xtalforge_core::{Lattice, Site, StructureBuilder};

fn carbon() -> xtalforge_core::Structure {
    StructureBuilder::new("C", Lattice::cubic(4.0).unwrap())
        .site(Site::new("C1", "C", [0.0, 0.0, 0.0]).unwrap())
        .metadata("source", json!("test"))
        .build()
        .unwrap()
}

#[test]
fn diagonal_supercell_is_p1_and_scales_the_cell() {
    let cell = carbon().supercell(2, 2, 2).unwrap();
    assert_eq!(cell.sites().len(), 8);
    assert!(cell.is_p1());
    assert!(cell.operations().is_empty());
    assert!((cell.lattice().volume() - 512.0).abs() < 1e-9);
    assert_eq!(cell.metadata()["source"], json!("test"));
}

#[test]
fn general_supercell_uses_the_determinant() {
    let cell = carbon()
        .supercell_matrix([[2.0, 1.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]])
        .unwrap();
    assert_eq!(cell.sites().len(), 2);
    assert!((cell.lattice().volume() - 128.0).abs() < 1e-9);
}

#[test]
fn origin_shift_rewraps_the_sites() {
    let lattice = Lattice::cubic(4.0).unwrap();
    let structure = StructureBuilder::new("shift", lattice)
        .site(Site::new("C1", "C", [0.1, 0.2, 0.3]).unwrap())
        .build()
        .unwrap();
    let shifted = structure.origin_shift([0.95, 0.0, 0.0]).unwrap();
    assert!((shifted.sites()[0].frac()[0] - 0.05).abs() < 1e-12);
    assert!((shifted.sites()[0].frac()[1] - 0.2).abs() < 1e-12);
}

#[test]
fn supercell_repeats_must_be_positive() {
    assert!(carbon().supercell(0, 1, 1).is_err());
}

#[test]
fn supercell_matrix_must_be_integer() {
    assert!(
        carbon()
            .supercell_matrix([[1.5, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]])
            .is_err()
    );
}

#[test]
fn hard_limit_refuses_a_huge_supercell() {
    let error = carbon()
        .supercell(1000, 1000, 1000)
        .expect_err("limit must refuse");
    assert_eq!(error.code(), "limit_exceeded");
}
