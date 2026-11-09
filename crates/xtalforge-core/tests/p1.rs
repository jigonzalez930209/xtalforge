//! Integration tests for symmetry operations and P1 expansion (Task 1.6).

#![allow(
    clippy::unwrap_used,
    clippy::expect_used,
    clippy::panic,
    clippy::indexing_slicing
)]

use xtalforge_core::{Lattice, Site, StructureBuilder, SymOp};

#[test]
fn symop_parses_and_applies_triplets() {
    let operation = SymOp::parse("-x+1/2,y,z").unwrap();
    let point = operation.apply([0.25, 0.1, 0.2]);
    assert!((point[0] - 0.25).abs() < 1e-12);
    assert!((point[1] - 0.1).abs() < 1e-12);
    assert!((point[2] - 0.2).abs() < 1e-12);
    assert_eq!(operation.triplet(), "-x+1/2,y,z");
    assert_eq!(SymOp::parse("x,y,z").unwrap(), SymOp::identity());
}

#[test]
fn symop_composition_and_inverse() {
    let operation = SymOp::parse("-x+1/2,y,z").unwrap();
    let inverse = operation.inverse().unwrap();
    let identity = operation.compose(&inverse).unwrap();
    assert!(identity.is_identity());
    let inversion = SymOp::parse("-x,-y,-z").unwrap();
    assert!(inversion.compose(&inversion).unwrap().is_identity());
}

#[test]
fn p1_expands_with_provenance() {
    let inversion = SymOp::parse("-x,-y,-z").unwrap();
    let lattice = Lattice::cubic(5.0).unwrap();
    let structure = StructureBuilder::new("test", lattice)
        .hall("P -1")
        .operation(SymOp::identity())
        .operation(inversion)
        .site(Site::new("C1", "C", [0.1, 0.2, 0.3]).unwrap())
        .build()
        .unwrap();
    let cell = structure.p1().unwrap();
    assert_eq!(cell.n_atoms(), 2);
    assert_eq!(cell.images_of(0).len(), 2);
    for atom in &cell.atoms {
        let operation = if atom.op_idx == 0 {
            SymOp::identity()
        } else {
            inversion
        };
        let generated = operation.apply(structure.sites()[0].frac());
        let shifted = [
            generated[0] + f64::from(atom.tau[0]),
            generated[1] + f64::from(atom.tau[1]),
            generated[2] + f64::from(atom.tau[2]),
        ];
        for (actual, expected) in shifted.iter().zip(atom.frac.iter()) {
            assert!((actual - expected).abs() < 1e-12);
        }
    }
}

#[test]
fn special_positions_are_deduplicated() {
    let inversion = SymOp::parse("-x,-y,-z").unwrap();
    let lattice = Lattice::cubic(5.0).unwrap();
    let structure = StructureBuilder::new("special", lattice)
        .hall("P -1")
        .operation(SymOp::identity())
        .operation(inversion)
        .site(Site::new("C1", "C", [0.0, 0.0, 0.0]).unwrap())
        .build()
        .unwrap();
    let cell = structure.p1().unwrap();
    assert_eq!(cell.n_atoms(), 1);
}

#[test]
fn p1_without_operations_is_the_identity() {
    let lattice = Lattice::cubic(4.0).unwrap();
    let structure = StructureBuilder::new("p1", lattice)
        .site(Site::new("C1", "C", [0.1, 0.2, 0.3]).unwrap())
        .site(Site::new("O1", "O", [0.4, 0.5, 0.6]).unwrap())
        .build()
        .unwrap();
    let cell = structure.p1().unwrap();
    assert_eq!(cell.n_atoms(), 2);
    assert_eq!(cell.elements(), vec!["C".to_owned(), "O".to_owned()]);
    assert_eq!(cell.occupancy(), vec![1.0, 1.0]);
}

#[test]
fn p1_tolerance_is_validated() {
    let lattice = Lattice::cubic(4.0).unwrap();
    let structure = StructureBuilder::new("p1", lattice).build().unwrap();
    assert!(xtalforge_core::P1Cell::expand(&structure, -1.0).is_err());
    assert!(xtalforge_core::P1Cell::expand(&structure, f64::NAN).is_err());
}

#[test]
fn disorder_sites_at_the_same_position_both_survive() {
    let lattice = Lattice::cubic(10.0).unwrap();
    let structure = StructureBuilder::new("disorder", lattice)
        .site(Site::new("C1", "C", [0.1, 0.2, 0.3]).unwrap())
        .site(Site::new("C2", "C", [0.1, 0.2, 0.3]).unwrap())
        .build()
        .unwrap();
    let p1 = structure.p1().unwrap();
    assert_eq!(p1.n_atoms(), 2);
    assert_eq!(p1.atoms[0].site_idx, 0);
    assert_eq!(p1.atoms[1].site_idx, 1);
}
