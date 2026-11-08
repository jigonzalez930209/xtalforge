//! Merged nets (Task 5.7.2).
//!
//! A merged net is two or three edge-transitive nets that share their
//! vertices. [`merge_nets`] keeps the two graphs as disjoint node sets in one
//! [`QuotientGraph`]: the second graph is copied with a node offset and its
//! labels are suffixed, so the connected components of the result are exactly
//! the components of the inputs.

use crate::nets::{NetEdge, QuotientGraph};

/// Concatenates two quotient graphs into one graph with disjoint nodes.
///
/// Nodes of `second` are appended after those of `first`; every edge of
/// `second` is offset by `first.nodes.len()`; labels of `second` get a
/// trailing apostrophe. Edges are sorted, so the output is deterministic.
/// The connected components of the result are the union of the components of
/// the two inputs.
#[must_use]
pub fn merge_nets(first: &QuotientGraph, second: &QuotientGraph) -> QuotientGraph {
    let offset = first.nodes.len();
    let mut nodes = first.nodes.clone();
    nodes.extend_from_slice(&second.nodes);
    let mut labels = first.labels.clone();
    labels.extend(second.labels.iter().map(|label| format!("{label}'")));
    let mut edges = first.edges.clone();
    edges.extend(second.edges.iter().map(|edge| NetEdge {
        a: edge.a + offset,
        b: edge.b + offset,
        shift: edge.shift,
    }));
    edges.sort();
    QuotientGraph {
        nodes,
        edges,
        labels,
    }
}
