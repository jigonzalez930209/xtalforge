//! Integration tests for supercells and origin shifts (Task 1.7).

#![allow(
    clippy::unwrap_used,
    clippy::expect_used,
    clippy::panic,
    clippy::indexing_slicing
)]

use serde_json::json;
use xtalforge_core::{Lattice, Site, StructureBuilder, SymOp};

fn carbon() -> xtalforge_core::Structure {
    StructureBuilder::new("C", Lattice::cubic(4.0).unwrap())
        .site(Site::new("C1", "C", [0.0, 0.0, 0.0]).unwrap())
        .metadata("source", json!("test"))
        .build()
        .unwrap()
}

#[test]
fn supercell_expands_the_whole_cell_before_repeating() {
    let inversion = SymOp::new([[-1, 0, 0], [0, -1, 0], [0, 0, -1]], [0.0, 0.0, 0.0]).unwrap();
    let asymmetric = StructureBuilder::new("Si", Lattice::cubic(5.430941).unwrap())
        .operation(SymOp::identity())
        .operation(inversion)
        .site(Site::new("Si1", "Si", [0.1, 0.2, 0.3]).unwrap())
        .build()
        .unwrap();
    assert_eq!(asymmetric.p1().unwrap().n_atoms(), 2);
    let cell = asymmetric.supercell(2, 1, 1).unwrap();
    assert_eq!(cell.sites().len(), 4);
    assert!(cell.is_p1());
    assert!(cell.operations().is_empty());
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
fn large_diagonal_supercell_keeps_every_repeat() {
    let cell = carbon().supercell(8, 8, 8).unwrap();
    assert_eq!(cell.sites().len(), 512);
    assert!((cell.lattice().volume() - 32_768.0).abs() < 1e-6);
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
