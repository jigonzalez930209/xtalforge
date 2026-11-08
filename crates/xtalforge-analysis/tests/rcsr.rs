//! Integration tests for RCSR identification (Tasks 5.6.1, 5.6.4, 5.6.5).

#![allow(
    clippy::unwrap_used,
    clippy::expect_used,
    clippy::panic,
    clippy::indexing_slicing
)]

use xtalforge_analysis::rcsr::{RCSR_LICENSE_NOTE, RcsrIndex};
use xtalforge_analysis::{NetEdge, QuotientGraph};

fn graph(node_count: usize, edges: &[(usize, usize)]) -> QuotientGraph {
    let nodes = (0..node_count)
        .map(|index| [index as f64 / node_count.max(1) as f64, 0.0, 0.0])
        .collect();
    QuotientGraph {
        nodes,
        edges: edges
            .iter()
            .map(|&(a, b)| NetEdge {
                a,
                b,
                shift: [0; 3],
            })
            .collect(),
        labels: (0..node_count).map(|index| format!("N{index}")).collect(),
    }
}

#[test]
fn common_nets_are_identified() {
    let index = RcsrIndex::from_subset();
    let pcu = graph(1, &[(0, 0), (0, 0), (0, 0)]);
    let dia = graph(2, &[(0, 1), (0, 1), (0, 1), (0, 1)]);
    let sql = graph(1, &[(0, 0), (0, 0)]);
    assert_eq!(index.identify(&pcu), Some("pcu"));
    assert_eq!(index.identify(&dia), Some("dia"));
    assert_eq!(index.identify(&sql), Some("sql"));
}

#[test]
fn an_unknown_graph_returns_none() {
    let index = RcsrIndex::from_subset();
    let unknown = graph(1, &[(0, 0)]);
    assert_eq!(index.identify(&unknown), None);
}

#[test]
fn the_subset_validates_and_is_unambiguous() {
    let index = RcsrIndex::from_subset();
    let count = index.validate().unwrap();
    assert!(count >= 10, "only {count} subset entries");
    assert_eq!(count, index.entries.len());
    assert_eq!(index.invariant_entries().len(), count);
}

#[test]
fn the_licence_note_names_the_external_bundle() {
    assert!(!RCSR_LICENSE_NOTE.trim().is_empty());
    assert!(RCSR_LICENSE_NOTE.contains("external"));
    assert!(RCSR_LICENSE_NOTE.contains("RCSR"));
}
