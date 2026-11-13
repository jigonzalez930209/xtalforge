//! Integration tests for structure-to-net contraction (Task 14.4.2).

#![allow(
    clippy::unwrap_used,
    clippy::expect_used,
    clippy::panic,
    clippy::indexing_slicing
)]

use xtalforge_analysis::contraction::{ContractionOptions, NetNodeRule, NodeKind, contract};
use xtalforge_core::{Lattice, Site, Structure, StructureBuilder, SymOp, TopologyBond};

fn carbon_pair() -> Structure {
    StructureBuilder::new("pair", Lattice::cubic(4.0).unwrap())
        .site(Site::new("C1", "C", [0.0, 0.0, 0.0]).unwrap())
        .site(Site::new("C2", "C", [0.35, 0.0, 0.0]).unwrap())
        .build()
        .unwrap()
}

fn carbon_chain() -> Structure {
    StructureBuilder::new("chain", Lattice::cubic(6.0).unwrap())
        .site(Site::new("C1", "C", [0.0, 0.0, 0.0]).unwrap())
        .site(Site::new("C2", "C", [0.25, 0.0, 0.0]).unwrap())
        .site(Site::new("C3", "C", [0.50, 0.0, 0.0]).unwrap())
        .build()
        .unwrap()
}

fn zinc_chain() -> Structure {
    StructureBuilder::new("zn-chain", Lattice::cubic(12.0).unwrap())
        .site(Site::new("Zn1", "Zn", [0.0, 0.0, 0.0]).unwrap())
        .site(Site::new("O1", "O", [0.125, 0.0, 0.0]).unwrap())
        .site(Site::new("C1", "C", [0.208, 0.0, 0.0]).unwrap())
        .site(Site::new("C2", "C", [0.292, 0.0, 0.0]).unwrap())
        .site(Site::new("O2", "O", [0.375, 0.0, 0.0]).unwrap())
        .site(Site::new("Zn2", "Zn", [0.5, 0.0, 0.0]).unwrap())
        .build()
        .unwrap()
}

#[test]
fn drawn_topology_contracts_to_its_bonds() {
    let structure = StructureBuilder::new("drawn", Lattice::cubic(10.0).unwrap())
        .site(Site::new("C1", "C", [0.0, 0.0, 0.0]).unwrap())
        .site(Site::new("C2", "C", [0.5, 0.0, 0.0]).unwrap())
        .topology_bond(TopologyBond {
            i: 0,
            j: 1,
            image: [0, 0, 1],
        })
        .build()
        .unwrap();
    let contracted = contract(
        &structure,
        &ContractionOptions {
            rule: NetNodeRule::Drawn,
            ..ContractionOptions::default()
        },
    )
    .unwrap();
    assert_eq!(contracted.net.nodes.len(), 2);
    assert_eq!(contracted.net.edges.len(), 1);
    assert_eq!(contracted.net.edges[0].shift, [0, 0, 1]);
    assert!(
        contracted
            .nodes
            .iter()
            .all(|node| node.kind == NodeKind::Drawn)
    );
    assert_eq!(contracted.invariants().coordination, vec![1, 1]);
}

#[test]
fn atom_rule_follows_perceived_bonds() {
    let contracted = contract(
        &carbon_pair(),
        &ContractionOptions {
            rule: NetNodeRule::Atoms,
            ..ContractionOptions::default()
        },
    )
    .unwrap();
    assert_eq!(contracted.net.nodes.len(), 2);
    assert_eq!(contracted.net.edges.len(), 1);
    assert!(
        contracted
            .nodes
            .iter()
            .all(|node| node.kind == NodeKind::Atom)
    );
}

#[test]
fn two_connected_atoms_contract_away() {
    let options = ContractionOptions {
        rule: NetNodeRule::Atoms,
        contract_degree_two: true,
        ..ContractionOptions::default()
    };
    let contracted = contract(&carbon_chain(), &options).unwrap();
    assert_eq!(contracted.net.nodes.len(), 2);
    assert_eq!(contracted.net.edges.len(), 1);
    assert_eq!(contracted.stripped.len(), 1);
    assert_eq!(contracted.nodes.len(), contracted.net.nodes.len());
    assert_eq!(contracted.edges.len(), contracted.net.edges.len());
    assert!(
        contracted
            .nodes
            .iter()
            .all(|node| node.kind == NodeKind::Atom)
    );
}

#[test]
fn all_node_rule_contracts_a_metal_linker_chain() {
    let contracted = contract(&zinc_chain(), &ContractionOptions::default()).unwrap();
    assert_eq!(contracted.net.nodes.len(), 2);
    assert_eq!(contracted.net.edges.len(), 1);
    assert!(
        contracted
            .nodes
            .iter()
            .all(|node| node.kind == NodeKind::MetalCluster)
    );
    assert!(!contracted.edges[0].atoms.is_empty());
    assert_eq!(contracted.invariants().coordination, vec![1, 1]);
}

#[test]
fn contraction_helpers_and_errors() {
    let atoms = ContractionOptions {
        rule: NetNodeRule::Atoms,
        ..ContractionOptions::default()
    };
    let contracted = contract(&carbon_pair(), &atoms).unwrap();
    assert!(contracted.vertex_symbols(12).is_ok());
    assert!(contracted.identify().is_ok());

    let drawn = ContractionOptions {
        rule: NetNodeRule::Drawn,
        ..ContractionOptions::default()
    };
    assert!(contract(&carbon_pair(), &drawn).is_err());

    let symmetric = StructureBuilder::new("symmetric", Lattice::cubic(10.0).unwrap())
        .site(Site::new("C1", "C", [0.0, 0.0, 0.0]).unwrap())
        .operation(SymOp::new([[-1, 0, 0], [0, 1, 0], [0, 0, 1]], [0.0, 0.0, 0.0]).unwrap())
        .topology_bond(TopologyBond {
            i: 0,
            j: 0,
            image: [0, 0, 0],
        })
        .build()
        .unwrap();
    assert!(contract(&symmetric, &drawn).is_err());

    let no_metal = ContractionOptions::default();
    assert!(contract(&carbon_pair(), &no_metal).is_err());

    let isolated = StructureBuilder::new("isolated", Lattice::cubic(10.0).unwrap())
        .site(Site::new("C1", "C", [0.0, 0.0, 0.0]).unwrap())
        .build()
        .unwrap();
    let atoms = ContractionOptions {
        rule: NetNodeRule::Atoms,
        ..ContractionOptions::default()
    };
    assert!(contract(&isolated, &atoms).is_err());

    assert!(
        contract(
            &carbon_pair(),
            &ContractionOptions {
                tolerance: 0.0,
                ..ContractionOptions::default()
            }
        )
        .is_err()
    );
}
