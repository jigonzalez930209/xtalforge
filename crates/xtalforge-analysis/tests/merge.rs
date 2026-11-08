//! Integration tests for merged nets (Task 5.7.2).

#![allow(
    clippy::unwrap_used,
    clippy::expect_used,
    clippy::panic,
    clippy::indexing_slicing
)]

use xtalforge_analysis::{NetEdge, QuotientGraph, merge_nets};

fn square() -> QuotientGraph {
    let nodes = vec![
        [0.0, 0.0, 0.0],
        [0.5, 0.0, 0.0],
        [0.5, 0.5, 0.0],
        [0.0, 0.5, 0.0],
    ];
    let labels = vec![
        "A".to_owned(),
        "B".to_owned(),
        "C".to_owned(),
        "D".to_owned(),
    ];
    let edges = vec![
        NetEdge {
            a: 0,
            b: 1,
            shift: [0; 3],
        },
        NetEdge {
            a: 1,
            b: 2,
            shift: [0; 3],
        },
        NetEdge {
            a: 2,
            b: 3,
            shift: [0; 3],
        },
        NetEdge {
            a: 0,
            b: 3,
            shift: [0; 3],
        },
        NetEdge {
            a: 0,
            b: 0,
            shift: [0, 0, 1],
        },
    ];
    QuotientGraph {
        nodes,
        edges,
        labels,
    }
}

#[test]
fn two_squares_keep_their_components_and_double_their_counts() {
    let merged = merge_nets(&square(), &square());
    assert_eq!(merged.nodes.len(), 8);
    assert_eq!(merged.edges.len(), 10);
    let invariants = merged.invariants();
    assert_eq!(invariants.nodes, 8);
    assert_eq!(invariants.edges, 10);
    assert_eq!(merged.components().len(), 2);
    assert_eq!(merged.labels[4], "A'");
    assert_eq!(merged.labels[7], "D'");
    assert_eq!(
        merged
            .edges
            .iter()
            .filter(|edge| edge.a >= 4 && edge.b >= 4)
            .count(),
        5
    );
    assert!(merged.edges.iter().all(|edge| edge.a < 8 && edge.b < 8));
}

#[test]
fn merging_is_deterministic_and_order_sensitive() {
    let first = merge_nets(&square(), &square());
    let second = merge_nets(&square(), &square());
    assert_eq!(first, second);
    let offset = merge_nets(&square(), &square());
    assert!(offset.edges.contains(&NetEdge {
        a: 4,
        b: 5,
        shift: [0; 3],
    }));
}
