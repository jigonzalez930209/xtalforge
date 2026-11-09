//! Integration tests for explicit and topology bonds (Task 1.11.2/1.11.3).

#![allow(
    clippy::unwrap_used,
    clippy::expect_used,
    clippy::panic,
    clippy::indexing_slicing
)]

use xtalforge_core::{ExplicitBond, Lattice, Site, StructureBuilder, TopologyBond};

#[test]
fn explicit_bond_is_added_beyond_perception() {
    let structure = StructureBuilder::new("far", Lattice::cubic(10.0).unwrap())
        .site(Site::new("C1", "C", [0.0, 0.0, 0.0]).unwrap())
        .site(Site::new("C2", "C", [0.5, 0.0, 0.0]).unwrap())
        .explicit_bond(ExplicitBond {
            i: 0,
            j: 1,
            image: [0, 0, 0],
            order: 2.0,
            remove: false,
        })
        .build()
        .unwrap();
    let graph = structure.bond_graph().unwrap();
    assert_eq!(graph.bonds.len(), 1);
    assert!(graph.bonds[0].explicit);
    assert!((graph.bonds[0].order - 2.0).abs() < 1e-12);
    assert!((graph.bonds[0].distance - 5.0).abs() < 1e-12);
}

#[test]
fn explicit_removal_removes_a_perceived_bond() {
    let structure = StructureBuilder::new("diamond", Lattice::cubic(3.566_986).unwrap())
        .site(Site::new("C1", "C", [0.0, 0.0, 0.0]).unwrap())
        .site(Site::new("C2", "C", [0.25, 0.25, 0.25]).unwrap())
        .explicit_bond(ExplicitBond {
            i: 0,
            j: 1,
            image: [0, 0, 0],
            order: 1.0,
            remove: true,
        })
        .build()
        .unwrap();
    assert!(structure.bond_graph().unwrap().bonds.is_empty());
}

#[test]
fn topology_bonds_are_stored_separately() {
    let structure = StructureBuilder::new("net", Lattice::cubic(10.0).unwrap())
        .topology_bond(TopologyBond {
            i: 0,
            j: 1,
            image: [0, 0, 0],
        })
        .build()
        .unwrap();
    assert_eq!(structure.topology_bonds().len(), 1);
    assert!(structure.bond_graph().unwrap().bonds.is_empty());
}

#[test]
fn explicit_bond_uses_its_lattice_image_for_the_length() {
    let structure = StructureBuilder::new("image", Lattice::cubic(10.0).unwrap())
        .site(Site::new("C1", "C", [0.0, 0.0, 0.0]).unwrap())
        .site(Site::new("C2", "C", [0.9, 0.0, 0.0]).unwrap())
        .explicit_bond(ExplicitBond {
            i: 0,
            j: 1,
            image: [-1, 0, 0],
            order: 1.0,
            remove: false,
        })
        .build()
        .unwrap();
    let graph = structure.bond_graph().unwrap();
    let bond = graph
        .bonds
        .iter()
        .find(|bond| bond.explicit)
        .expect("explicit bond");
    assert_eq!((bond.i, bond.j, bond.image), (0, 1, [-1, 0, 0]));
    assert!((bond.distance - 1.0).abs() < 1e-12);
}
