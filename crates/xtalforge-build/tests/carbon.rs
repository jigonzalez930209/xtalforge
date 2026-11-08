//! Integration tests for graphene and nanotube builders (Tasks 6.5, 6.6).

#![allow(
    clippy::unwrap_used,
    clippy::expect_used,
    clippy::panic,
    clippy::indexing_slicing
)]

use xtalforge_build::{graphene, nanotube};

#[test]
fn graphene_sheet_has_two_atoms_per_cell() {
    let sheet = graphene(3, 2, 1, 3.35, 15.0).unwrap();
    assert_eq!(sheet.sites().len(), 12);
    assert!(sheet.is_p1());
    let [a, _b, c] = sheet.lattice().lengths();
    assert!((a - 3.0 * 2.4595).abs() < 1e-2);
    assert!((c - 18.35).abs() < 1e-9);
    let graph = sheet.bond_graph().unwrap();
    assert!(!graph.bonds.is_empty());
}

#[test]
fn stacked_graphene_has_layers() {
    let sheet = graphene(1, 1, 3, 3.35, 10.0).unwrap();
    assert_eq!(sheet.sites().len(), 6);
    assert!((sheet.lattice().lengths()[2] - 20.05).abs() < 1e-9);
}

#[test]
fn nanotube_rolls_a_chiral_sheet() {
    let tube = nanotube(5, 5, 2, 10.0).unwrap();
    assert!(tube.sites().len() > 10);
    let [a, _b, c] = tube.lattice().lengths();
    assert!(a > 0.0);
    assert!(c > 10.0);
    assert!(nanotube(0, 0, 1, 10.0).is_err());
    assert!(nanotube(5, 5, 0, 10.0).is_err());
}
