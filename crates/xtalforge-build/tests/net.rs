//! Integration tests for the net-as-structure builder (Task 6.10.1).

#![allow(
    clippy::unwrap_used,
    clippy::expect_used,
    clippy::panic,
    clippy::indexing_slicing
)]

use xtalforge_build::build_net;
use xtalforge_io::cgd::{CgdEntry, CgdNode};

#[test]
fn square_net_becomes_sites_and_topology_bonds() {
    let entry = CgdEntry {
        name: "square".to_owned(),
        group: "P1".to_owned(),
        nodes: vec![
            CgdNode {
                label: "A".to_owned(),
                coordination: 4,
                frac: vec![0.0, 0.0, 0.0],
            },
            CgdNode {
                label: "B".to_owned(),
                coordination: 4,
                frac: vec![0.5, 0.0, 0.0],
            },
            CgdNode {
                label: "C".to_owned(),
                coordination: 4,
                frac: vec![0.5, 0.5, 0.0],
            },
            CgdNode {
                label: "D".to_owned(),
                coordination: 4,
                frac: vec![0.0, 0.5, 0.0],
            },
        ],
        edges: vec![
            (vec![0.0, 0.0, 0.0], vec![0.5, 0.0, 0.0]),
            (vec![0.5, 0.5, 0.0], vec![0.0, 0.5, 0.0]),
        ],
        cell: Some(vec![4.0, 4.0, 4.0, 90.0, 90.0, 90.0]),
        line: 1,
    };
    let structure = build_net(&entry, "C").unwrap();
    assert_eq!(structure.sites().len(), 4);
    assert_eq!(structure.topology_bonds().len(), 2);
    assert!((structure.lattice().lengths()[0] - 4.0).abs() < 1e-9);
}
