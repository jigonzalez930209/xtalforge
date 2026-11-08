//! Integration tests for periodic nets (Task 5.6).

#![allow(
    clippy::unwrap_used,
    clippy::expect_used,
    clippy::panic,
    clippy::indexing_slicing
)]

use xtalforge_analysis::quotient_graph;
use xtalforge_io::cgd::{CgdEntry, CgdNode};

fn square_net() -> CgdEntry {
    CgdEntry {
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
            (vec![0.5, 0.0, 0.0], vec![0.5, 0.5, 0.0]),
            (vec![0.5, 0.5, 0.0], vec![0.0, 0.5, 0.0]),
            (vec![0.0, 0.5, 0.0], vec![0.0, 0.0, 0.0]),
            (vec![0.0, 0.0, 0.0], vec![0.0, 0.0, 1.0]),
        ],
        cell: None,
        line: 1,
    }
}

#[test]
fn quotient_graph_matches_nodes_and_shifts() {
    let graph = quotient_graph(&square_net(), 1e-6).unwrap();
    assert_eq!(graph.nodes.len(), 4);
    assert_eq!(graph.edges.len(), 5);
    let invariants = graph.invariants();
    assert_eq!(invariants.cycle_rank, 2);
    assert_eq!(invariants.coordination, vec![4, 2, 2, 2]);
    let shifts = graph.shift_map();
    assert!(shifts.values().any(|values| values.contains(&[0, 0, 1])));
}

#[test]
fn bad_inputs_are_rejected() {
    let mut entry = square_net();
    entry.nodes.clear();
    assert!(quotient_graph(&entry, 1e-6).is_err());
    let mut entry = square_net();
    entry
        .edges
        .push((vec![0.25, 0.25, 0.0], vec![0.0, 0.0, 0.0]));
    assert!(quotient_graph(&entry, 1e-6).is_err());
    assert!(quotient_graph(&square_net(), 0.0).is_err());
}

#[test]
fn net_index_and_components() {
    let file = xtalforge_io::cgd::CgdFile {
        entries: vec![square_net()],
        problems: Vec::new(),
    };
    let index = xtalforge_analysis::NetIndex::from_file(&file, 1e-6);
    assert_eq!(index.entries.len(), 1);
    let graph = quotient_graph(&square_net(), 1e-6).unwrap();
    let matches = index.search(&graph.invariants());
    assert_eq!(matches, vec!["square"]);
    assert_eq!(graph.components(), vec![vec![0, 1, 2, 3]]);
}

#[test]
fn interpenetration_creates_two_components() {
    let graph = quotient_graph(&square_net(), 1e-6).unwrap();
    let double = xtalforge_analysis::interpenetrate(&graph, [0, 0, 1]).unwrap();
    assert_eq!(double.nodes.len(), 8);
    assert_eq!(double.edges.len(), 10);
    assert_eq!(double.components().len(), 2);
    assert!(xtalforge_analysis::interpenetrate(&graph, [0, 0, 1]).is_ok());
}
