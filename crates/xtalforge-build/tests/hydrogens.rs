//! Integration tests for hydrogen placement (Task 4.6.1).

#![allow(
    clippy::unwrap_used,
    clippy::expect_used,
    clippy::panic,
    clippy::indexing_slicing
)]

use std::collections::BTreeSet;

use xtalforge_build::{HydrogenOptions, add_hydrogens};
use xtalforge_core::{ExplicitBond, Lattice, Site, Structure, StructureBuilder, SymOp};

fn single(element: &str) -> Structure {
    let points = [(element.to_owned(), [0.0, 0.0, 0.0])];
    Structure::molecule("molecule", &points, 6.0).unwrap()
}

fn hydrogens(structure: &Structure) -> Vec<&Site> {
    structure
        .sites()
        .iter()
        .filter(|site| site.element() == "H")
        .collect()
}

fn nearest_heavy(structure: &Structure, frac: [f64; 3]) -> f64 {
    structure
        .sites()
        .iter()
        .filter(|site| site.element() != "H")
        .map(|site| structure.lattice().distance(site.frac(), frac))
        .fold(f64::INFINITY, f64::min)
}

#[test]
fn methane_gets_four_hydrogens_and_four_bonds() {
    let structure = single("C");
    let filled = add_hydrogens(&structure, &HydrogenOptions::default()).unwrap();
    let added = hydrogens(&filled);
    assert_eq!(added.len(), 4);
    assert_eq!(filled.sites().len(), 5);
    assert_eq!(filled.bond_graph().unwrap().bonds.len(), 4);
    for hydrogen in added {
        let distance = nearest_heavy(&filled, hydrogen.frac());
        assert!((distance - 1.09).abs() < 1e-6, "distance {distance}");
    }
}

#[test]
fn water_gets_two_hydrogens() {
    let structure = single("O");
    let filled = add_hydrogens(&structure, &HydrogenOptions::default()).unwrap();
    let added = hydrogens(&filled);
    assert_eq!(added.len(), 2);
    assert_eq!(filled.bond_graph().unwrap().bonds.len(), 2);
    for hydrogen in added {
        let distance = nearest_heavy(&filled, hydrogen.frac());
        assert!((distance - 0.96).abs() < 1e-6, "distance {distance}");
    }
}

#[test]
fn ammonia_gets_three_hydrogens() {
    let structure = single("N");
    let filled = add_hydrogens(&structure, &HydrogenOptions::default()).unwrap();
    let added = hydrogens(&filled);
    assert_eq!(added.len(), 3);
    assert_eq!(filled.bond_graph().unwrap().bonds.len(), 3);
    for hydrogen in added {
        let distance = nearest_heavy(&filled, hydrogen.frac());
        assert!((distance - 1.01).abs() < 1e-6, "distance {distance}");
    }
}

#[test]
fn ethene_fragment_gets_four_hydrogens() {
    let structure = StructureBuilder::new("ethene", Lattice::cubic(20.0).unwrap())
        .site(Site::new("C1", "C", [0.1, 0.1, 0.1]).unwrap())
        .site(Site::new("C2", "C", [0.167, 0.1, 0.1]).unwrap())
        .explicit_bond(ExplicitBond {
            i: 0,
            j: 1,
            image: [0, 0, 0],
            order: 2.0,
            remove: false,
        })
        .build()
        .unwrap();
    let filled = add_hydrogens(&structure, &HydrogenOptions::default()).unwrap();
    assert_eq!(hydrogens(&filled).len(), 4);
    let graph = filled.bond_graph().unwrap();
    assert_eq!(graph.bonds.len(), 5);
    assert_eq!(graph.degree(0), 3);
    assert_eq!(graph.degree(1), 3);
}

#[test]
fn custom_distance_is_used_for_every_hydrogen() {
    let options = HydrogenOptions {
        distance: 1.2,
        ..HydrogenOptions::default()
    };
    let structure = single("C");
    let filled = add_hydrogens(&structure, &options).unwrap();
    for hydrogen in hydrogens(&filled) {
        let distance = nearest_heavy(&filled, hydrogen.frac());
        assert!((distance - 1.2).abs() < 1e-6, "distance {distance}");
    }
}

#[test]
fn existing_hydrogens_are_not_duplicated() {
    let points = [
        ("O".to_owned(), [0.0, 0.0, 0.0]),
        ("H".to_owned(), [0.96, 0.0, 0.0]),
        ("H".to_owned(), [0.0, 0.96, 0.0]),
    ];
    let structure = Structure::molecule("water", &points, 6.0).unwrap();
    let filled = add_hydrogens(&structure, &HydrogenOptions::default()).unwrap();
    assert_eq!(filled.sites().len(), 3);
    assert_eq!(hydrogens(&filled).len(), 2);
}

#[test]
fn new_labels_do_not_collide_with_existing_ones() {
    let points = [
        ("O".to_owned(), [0.0, 0.0, 0.0]),
        ("H".to_owned(), [0.96, 0.0, 0.0]),
    ];
    let structure = Structure::molecule("water", &points, 6.0).unwrap();
    let filled = add_hydrogens(&structure, &HydrogenOptions::default()).unwrap();
    assert_eq!(filled.sites().len(), 3);
    let labels: BTreeSet<_> = filled
        .sites()
        .iter()
        .map(|site| site.label().to_owned())
        .collect();
    assert_eq!(labels.len(), 3);
}

#[test]
fn invalid_options_are_rejected() {
    let structure = single("C");
    let options = HydrogenOptions {
        distance: f64::NAN,
        ..HydrogenOptions::default()
    };
    assert!(add_hydrogens(&structure, &options).is_err());
    let options = HydrogenOptions {
        distance: -1.0,
        ..HydrogenOptions::default()
    };
    assert!(add_hydrogens(&structure, &options).is_err());

    let mut options = HydrogenOptions::default();
    options.angles.insert("C".to_owned(), 200.0);
    assert!(add_hydrogens(&structure, &options).is_err());

    let mut options = HydrogenOptions::default();
    options.angles.insert("Fe".to_owned(), 120.0);
    assert!(add_hydrogens(&structure, &options).is_err());

    let mut options = HydrogenOptions::default();
    options.angles.insert("Xx".to_owned(), 120.0);
    assert!(add_hydrogens(&structure, &options).is_err());
}

#[test]
fn crystal_fields_are_preserved() {
    let lattice = Lattice::cubic(20.0).unwrap();
    let operation = SymOp::identity();
    let structure = StructureBuilder::new("kept", lattice)
        .hall("P 1")
        .operation(operation)
        .metadata("source", serde_json::json!("test"))
        .site(Site::new("C1", "C", [0.5, 0.5, 0.5]).unwrap())
        .build()
        .unwrap();
    let filled = add_hydrogens(&structure, &HydrogenOptions::default()).unwrap();
    assert_eq!(filled.name(), "kept");
    assert_eq!(filled.lattice(), &lattice);
    assert_eq!(filled.hall(), "P 1");
    assert_eq!(filled.operations(), &[operation]);
    assert_eq!(filled.metadata(), structure.metadata());
    assert_eq!(hydrogens(&filled).len(), 4);
}

#[test]
fn directions_follow_the_periodic_minimum_image() {
    let structure = StructureBuilder::new("periodic", Lattice::cubic(4.0).unwrap())
        .site(Site::new("O1", "O", [0.90, 0.5, 0.5]).unwrap())
        .site(Site::new("H2", "H", [0.14, 0.5, 0.5]).unwrap())
        .build()
        .unwrap();
    let filled = add_hydrogens(&structure, &HydrogenOptions::default()).unwrap();
    assert_eq!(hydrogens(&filled).len(), 2);
    let added = filled
        .sites()
        .iter()
        .find(|site| site.label() == "H1")
        .unwrap();
    let distance = filled.lattice().distance([0.90, 0.5, 0.5], added.frac());
    assert!((distance - 0.96).abs() < 1e-6, "distance {distance}");
}

#[test]
fn coincident_bonds_have_no_direction_basis() {
    let structure = StructureBuilder::new("coincident", Lattice::cubic(20.0).unwrap())
        .site(Site::new("C1", "C", [0.1, 0.1, 0.1]).unwrap())
        .site(Site::new("N1", "N", [0.1, 0.1, 0.1]).unwrap())
        .explicit_bond(ExplicitBond {
            i: 0,
            j: 1,
            image: [0, 0, 0],
            order: 1.0,
            remove: false,
        })
        .build()
        .unwrap();
    assert!(add_hydrogens(&structure, &HydrogenOptions::default()).is_err());
}
