//! Integration tests for bond perception and the connectivity graph.

#![allow(
    clippy::unwrap_used,
    clippy::expect_used,
    clippy::panic,
    clippy::indexing_slicing
)]

use xtalforge_core::{BondRules, Lattice, Site, StructureBuilder, perceive_bonds};

#[test]
fn carbon_carbon_bond_in_diamond() {
    let lattice = Lattice::cubic(3.566_986).unwrap();
    let structure = StructureBuilder::new("diamond", lattice)
        .site(Site::new("C1", "C", [0.0, 0.0, 0.0]).unwrap())
        .site(Site::new("C2", "C", [0.25, 0.25, 0.25]).unwrap())
        .build()
        .unwrap();
    let graph = structure.bond_graph().unwrap();
    assert_eq!(graph.bonds.len(), 1);
    assert!((graph.bonds[0].distance - 1.544_550_245_471_719_7).abs() < 1e-9);
    assert_eq!(graph.degree(0), 1);
    assert_eq!(graph.degree(1), 1);
    assert_eq!(graph.components(), vec![vec![0, 1]]);
}

#[test]
fn methane_has_four_bonds_and_one_component() {
    let lattice = Lattice::cubic(10.0).unwrap();
    let center = [5.0, 5.0, 5.0];
    let mut builder = StructureBuilder::new("methane", lattice)
        .site(Site::new("C1", "C", [0.5, 0.5, 0.5]).unwrap());
    for (index, direction) in [
        [1.0, 0.0, 0.0],
        [-1.0, 0.0, 0.0],
        [0.0, 1.0, 0.0],
        [0.0, -1.0, 0.0],
    ]
    .iter()
    .enumerate()
    {
        let point = [
            center[0] + 1.09 * direction[0],
            center[1] + 1.09 * direction[1],
            center[2] + 1.09 * direction[2],
        ];
        builder = builder.site(
            Site::new(
                format!("H{index}"),
                "H",
                [point[0] / 10.0, point[1] / 10.0, point[2] / 10.0],
            )
            .unwrap(),
        );
    }
    let graph = builder.build().unwrap().bond_graph().unwrap();
    assert_eq!(graph.bonds.len(), 4);
    assert_eq!(graph.degree(0), 4);
    assert_eq!(graph.components(), vec![vec![0, 1, 2, 3, 4]]);
}

#[test]
fn far_atoms_are_not_bonded() {
    let lattice = Lattice::cubic(10.0).unwrap();
    let structure = StructureBuilder::new("far", lattice)
        .site(Site::new("C1", "C", [0.0, 0.0, 0.0]).unwrap())
        .site(Site::new("C2", "C", [0.5, 0.0, 0.0]).unwrap())
        .build()
        .unwrap();
    assert!(structure.bond_graph().unwrap().bonds.is_empty());
}

#[test]
fn dummy_species_never_bonds() {
    let lattice = Lattice::cubic(10.0).unwrap();
    let structure = StructureBuilder::new("dummy", lattice)
        .site(Site::new("C1", "C", [0.5, 0.5, 0.5]).unwrap())
        .site(Site::new("X1", "X", [0.56, 0.5, 0.5]).unwrap())
        .build()
        .unwrap();
    assert!(structure.bond_graph().unwrap().bonds.is_empty());
}

#[test]
fn rules_are_validated() {
    let lattice = Lattice::cubic(4.0).unwrap();
    let structure = StructureBuilder::new("x", lattice).build().unwrap();
    assert!(
        perceive_bonds(
            &structure,
            BondRules {
                scale: 0.0,
                minimum_distance: 0.4,
                ..BondRules::default()
            }
        )
        .is_err()
    );
    assert!(
        perceive_bonds(
            &structure,
            BondRules {
                scale: 1.2,
                minimum_distance: -1.0,
                ..BondRules::default()
            }
        )
        .is_err()
    );
}

#[test]
fn metal_metal_pairs_are_not_bonds() {
    let lattice = Lattice::cubic(10.0).unwrap();
    let structure = StructureBuilder::new("cu2", lattice)
        .site(Site::new("Cu1", "Cu", [0.0, 0.0, 0.0]).unwrap())
        .site(Site::new("Cu2", "Cu", [0.2, 0.0, 0.0]).unwrap())
        .build()
        .unwrap();
    assert!(structure.bond_graph().unwrap().bonds.is_empty());
}

#[test]
fn carbon_hydrogen_uses_the_reference_scale() {
    let lattice = Lattice::cubic(10.0).unwrap();
    let structure = StructureBuilder::new("ch", lattice)
        .site(Site::new("C1", "C", [0.0, 0.0, 0.0]).unwrap())
        .site(Site::new("H1", "H", [0.119, 0.0, 0.0]).unwrap())
        .site(Site::new("H2", "H", [0.121, 0.0, 0.0]).unwrap())
        .build()
        .unwrap();
    let graph = structure.bond_graph().unwrap();
    assert_eq!(graph.bonds.len(), 1);
    assert_eq!((graph.bonds[0].i, graph.bonds[0].j), (0, 1));
}
