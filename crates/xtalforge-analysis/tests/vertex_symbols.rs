//! Integration tests for vertex symbols (Task 14.4.2).

#![allow(
    clippy::unwrap_used,
    clippy::expect_used,
    clippy::panic,
    clippy::indexing_slicing
)]

use xtalforge_analysis::QuotientGraph;
use xtalforge_analysis::vertex_symbols::{vertex_symbol_text, vertex_symbols};

fn key(key: &str) -> QuotientGraph {
    QuotientGraph::from_canonical_key(key).unwrap()
}

#[test]
fn pcu_vertex_symbol_counts_four_and_six_rings() {
    let pcu = key("3:1:0,0,-1,0,0;0,0,0,-1,0;0,0,0,0,-1");
    let symbols = vertex_symbols(&pcu, 12).unwrap();
    assert_eq!(symbols.len(), 1);
    let symbol = &symbols[0];
    assert_eq!(symbol.neighbours.len(), 6);
    assert_eq!(symbol.pairs.len(), 15);
    assert_eq!(symbol.shortest, Some(4));
    assert_eq!(symbol.edge_rings, vec![4; 6]);
    let fours = symbol.pairs.iter().filter(|ring| **ring == Some(4)).count();
    let sixes = symbol.pairs.iter().filter(|ring| **ring == Some(6)).count();
    assert_eq!((fours, sixes), (12, 3));
    assert_eq!(symbol.text().matches('4').count(), 12);
}

#[test]
fn dia_vertex_symbol_is_six_rings() {
    let dia = key("3:2:0,1,0,0,0;0,1,0,0,1;0,1,0,1,0;0,1,1,0,0");
    let symbols = vertex_symbols(&dia, 12).unwrap();
    assert_eq!(symbols.len(), 2);
    for symbol in &symbols {
        assert_eq!(symbol.shortest, Some(6));
        assert!(symbol.pairs.iter().all(|ring| *ring == Some(6)));
    }
    let text = vertex_symbol_text(&dia, 12).unwrap();
    assert_eq!(text.len(), 2);
    assert!(text.iter().all(|line| line.starts_with("6.6.6")));
}

#[test]
fn a_leaf_has_no_rings() {
    let chain = key("3:2:0,1,0,0,0");
    let symbols = vertex_symbols(&chain, 12).unwrap();
    for symbol in &symbols {
        assert!(symbol.pairs.is_empty());
        assert!(symbol.text().is_empty());
        assert_eq!(symbol.shortest, None);
    }
}

#[test]
fn vertex_symbol_errors_are_typed() {
    let empty = QuotientGraph {
        nodes: Vec::new(),
        edges: Vec::new(),
        labels: Vec::new(),
    };
    assert!(vertex_symbols(&empty, 12).is_err());
    let pcu = key("3:1:0,0,-1,0,0;0,0,0,-1,0;0,0,0,0,-1");
    let error = vertex_symbols(&pcu, 2).unwrap_err();
    assert_eq!(error.code(), "invalid_input");
}

#[test]
fn pair_lookup_uses_row_major_order() {
    let pcu = key("3:1:0,0,-1,0,0;0,0,0,-1,0;0,0,0,0,-1");
    let symbol = vertex_symbols(&pcu, 12).unwrap().remove(0);
    assert_eq!(symbol.pair_ring(0, 1), symbol.pairs[0]);
    assert_eq!(symbol.pair_ring(0, 2), Some(4));
    assert_eq!(symbol.pair_ring(1, 0), None);
    assert_eq!(symbol.pair_ring(0, 9), None);
}
