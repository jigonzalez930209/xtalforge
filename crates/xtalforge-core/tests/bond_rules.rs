//! Integration tests for bond rules and explicit-bond propagation (Task 14.3.2).

#![allow(
    clippy::unwrap_used,
    clippy::expect_used,
    clippy::panic,
    clippy::indexing_slicing
)]

use serde_json::json;
use xtalforge_core::{
    BondRules, ExplicitBond, Lattice, Site, Structure, StructureBuilder, SymOp,
    explicit_bond_images, perceive_bonds,
};

fn pair(a: &str, b: &str, first: [f64; 3], second: [f64; 3]) -> Structure {
    StructureBuilder::new("pair", Lattice::cubic(10.0).unwrap())
        .site(Site::new("A", a, first).unwrap())
        .site(Site::new("B", b, second).unwrap())
        .build()
        .unwrap()
}

#[test]
fn pair_range_replaces_the_global_criterion() {
    let structure = pair("C", "O", [0.0, 0.0, 0.0], [0.15, 0.0, 0.0]);
    let global = BondRules::new(0.5, 0.4);
    assert!(perceive_bonds(&structure, global).unwrap().bonds.is_empty());

    let rules = BondRules::new(0.5, 0.4).with_pair_range("C", "O", 1.4, 1.6);
    let graph = perceive_bonds(&structure, rules).unwrap();
    assert_eq!(graph.bonds.len(), 1);
    assert!((graph.bonds[0].distance - 1.5).abs() < 1e-12);

    let reversed = BondRules::new(0.5, 0.4).with_pair_range("O", "C", 1.4, 1.6);
    assert_eq!(perceive_bonds(&structure, reversed).unwrap().bonds.len(), 1);

    let narrow = BondRules::new(0.5, 0.4).with_pair_range("C", "O", 1.6, 1.7);
    assert!(perceive_bonds(&structure, narrow).unwrap().bonds.is_empty());
}

#[test]
fn bond_graph_uses_the_rules_stored_on_the_structure() {
    let rules = BondRules::new(0.5, 0.4).with_pair_range("C", "O", 1.8, 2.0);
    let plain = pair("C", "O", [0.0, 0.0, 0.0], [0.19, 0.0, 0.0]);
    assert!(plain.bond_graph().unwrap().bonds.is_empty());
    let stored = StructureBuilder::new("custom", Lattice::cubic(10.0).unwrap())
        .metadata("bond_rules", rules.to_json())
        .site(Site::new("C1", "C", [0.0, 0.0, 0.0]).unwrap())
        .site(Site::new("O1", "O", [0.19, 0.0, 0.0]).unwrap())
        .build()
        .unwrap();
    assert_eq!(stored.bond_graph().unwrap().bonds.len(), 1);
    let broken = StructureBuilder::new("broken", Lattice::cubic(10.0).unwrap())
        .metadata("bond_rules", json!({"scale": "wide"}))
        .site(Site::new("C1", "C", [0.0, 0.0, 0.0]).unwrap())
        .build()
        .unwrap();
    assert!(broken.bond_graph().is_err());
}

#[test]
fn metal_metal_switch_is_exposed() {
    let structure = pair("Cu", "Cu", [0.0, 0.0, 0.0], [0.2, 0.0, 0.0]);
    assert!(structure.bond_graph().unwrap().bonds.is_empty());
    let rules = BondRules::default().with_metal_metal(true);
    let graph = perceive_bonds(&structure, rules).unwrap();
    assert_eq!(graph.bonds.len(), 1);
    assert!((graph.bonds[0].distance - 2.0).abs() < 1e-12);
}

#[test]
fn rules_round_trip_through_json_and_metadata() {
    let rules = BondRules::default()
        .with_pair_range("Fe", "O", 1.6, 2.2)
        .with_metal_metal(true);
    let value = rules.to_json();
    assert_eq!(BondRules::from_json(&value).unwrap(), rules);
    assert_eq!(
        BondRules::from_json(&json!({})).unwrap(),
        BondRules::default()
    );

    let structure = StructureBuilder::new("stored", Lattice::cubic(10.0).unwrap())
        .metadata("bond_rules", value)
        .build()
        .unwrap();
    assert_eq!(BondRules::from_structure(&structure).unwrap(), rules);
    assert_eq!(
        BondRules::from_structure(&pair("C", "C", [0.0; 3], [0.5; 3])).unwrap(),
        BondRules::default()
    );

    assert!(BondRules::from_json(&json!({"scale": "wide"})).is_err());
    assert!(BondRules::from_json(&json!({"allow_metal_metal": "yes"})).is_err());
    assert!(BondRules::from_json(&json!({"pair_ranges": {"CO": [1.0, 2.0]}})).is_err());
    assert!(BondRules::from_json(&json!({"pair_ranges": {"C-O": [2.0, 1.0]}})).is_err());
    assert!(BondRules::new(0.0, 0.4).validate().is_err());
    assert!(BondRules::new(1.2, -1.0).validate().is_err());
}

fn inverted_pair(first: [f64; 3], second: [f64; 3]) -> Structure {
    StructureBuilder::new("inverted-pair", Lattice::cubic(8.0).unwrap())
        .hall("P -1")
        .operation(SymOp::identity())
        .operation(SymOp::parse("-x,-y,-z").unwrap())
        .site(Site::new("C1", "C", first).unwrap())
        .site(Site::new("S1", "S", second).unwrap())
        .explicit_bond(ExplicitBond {
            i: 0,
            j: 2,
            image: [1, 0, 0],
            order: 2.0,
            remove: false,
        })
        .build()
        .unwrap()
}

#[test]
fn explicit_bond_reaches_every_symmetry_image() {
    let structure = inverted_pair([0.1, 0.1, 0.1], [0.6, 0.1, 0.1]);
    let p1 = structure.p1().unwrap();
    assert_eq!(p1.n_atoms(), 4);
    let images = explicit_bond_images(&structure, &p1, &structure.explicit_bonds()[0]);
    assert_eq!(images.len(), 2);
    assert!(images.iter().all(|bond| (bond.order - 2.0).abs() < 1e-12));

    let graph = structure.bond_graph().unwrap();
    assert_eq!(graph.bonds.len(), 2);
    assert!(graph.bonds.iter().all(|bond| bond.explicit));
    assert!(graph.bonds.iter().any(|bond| bond.image == [1, 0, 0]));
    assert!(graph.bonds.iter().any(|bond| bond.image == [-1, 0, 0]));
}

#[test]
fn explicit_removal_reaches_every_symmetry_image() {
    let structure = StructureBuilder::new("removed", Lattice::cubic(8.0).unwrap())
        .hall("P -1")
        .operation(SymOp::identity())
        .operation(SymOp::parse("-x,-y,-z").unwrap())
        .site(Site::new("C1", "C", [0.1, 0.1, 0.1]).unwrap())
        .site(Site::new("C2", "C", [0.16, 0.1, 0.1]).unwrap())
        .explicit_bond(ExplicitBond {
            i: 0,
            j: 2,
            image: [0, 0, 0],
            order: 1.0,
            remove: true,
        })
        .build()
        .unwrap();
    assert_eq!(
        structure.bond_graph().unwrap().bonds.len(),
        0,
        "the suppression must remove both orbit edges"
    );
}

#[test]
fn explicit_images_follow_the_orbit_on_both_ends() {
    let structure = inverted_pair([0.1, 0.1, 0.1], [0.6, 0.1, 0.1]);
    let images = explicit_bond_images(
        &structure,
        &structure.p1().unwrap(),
        &structure.explicit_bonds()[0],
    );
    assert_eq!(images[0].i, 0);
    assert_eq!(images[1].i, 1);
    assert_eq!(images[0].j, 2);
    assert_eq!(images[1].j, 3);
}
