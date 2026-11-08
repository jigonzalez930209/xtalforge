//! Integration tests for the later-phase command wrappers (Task 3.10.7).

#![allow(
    clippy::unwrap_used,
    clippy::expect_used,
    clippy::panic,
    clippy::indexing_slicing
)]

use xtalforge_commands::{Command, Interpenetrate, MarkConnectionPoints, SetCharges};
use xtalforge_core::{Lattice, Site, StructureBuilder};

fn nacl() -> xtalforge_core::Structure {
    StructureBuilder::new("NaCl", Lattice::cubic(5.64).unwrap())
        .site(Site::new("Na1", "Na", [0.0, 0.0, 0.0]).unwrap())
        .site(Site::new("Cl1", "Cl", [0.5, 0.5, 0.5]).unwrap())
        .build()
        .unwrap()
}

#[test]
fn charges_are_written_and_validated() {
    let command = SetCharges::new(vec![Some(1.0), Some(-1.0)]);
    let charged = command.apply(&nacl()).unwrap();
    assert_eq!(charged.sites()[0].charge(), Some(1.0));
    assert_eq!(charged.sites()[1].charge(), Some(-1.0));
    assert!(SetCharges::new(vec![Some(1.0)]).apply(&nacl()).is_err());
}

#[test]
fn interpenetration_duplicates_sites() {
    let doubled = Interpenetrate::new([0.0, 0.0, 0.5]).apply(&nacl()).unwrap();
    assert_eq!(doubled.sites().len(), 4);
    assert_eq!(doubled.sites()[2].label(), "Na1'");
    assert!((doubled.sites()[2].frac()[2] - 0.5).abs() < 1e-12);
    assert!(
        Interpenetrate::new([f64::NAN, 0.0, 0.0])
            .apply(&nacl())
            .is_err()
    );
}

#[test]
fn connection_points_add_dummy_atoms() {
    let marked = MarkConnectionPoints::new(vec![[0.25, 0.25, 0.25]])
        .apply(&nacl())
        .unwrap();
    assert_eq!(marked.sites().len(), 3);
    assert_eq!(marked.sites()[2].element(), "X");
    assert!(
        MarkConnectionPoints::new(Vec::new())
            .apply(&nacl())
            .is_err()
    );
}
