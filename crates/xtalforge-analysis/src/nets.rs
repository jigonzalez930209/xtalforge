//! Periodic nets as quotient graphs (Task 5.6).

#![allow(clippy::indexing_slicing)] // fixed 3-coordinate vectors

use std::collections::BTreeMap;

use serde::{Deserialize, Serialize};
use xtalforge_core::Error;
use xtalforge_io::cgd::CgdEntry;

/// One edge of the quotient graph, with the lattice shift of `b`.
#[derive(Debug, Clone, PartialEq, Eq, PartialOrd, Ord, Serialize, Deserialize)]
pub struct NetEdge {
    /// First node.
    pub a: usize,
    /// Second node.
    pub b: usize,
    /// Lattice shift applied to `b`.
    pub shift: [i32; 3],
}

/// A net reduced to its quotient graph.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct QuotientGraph {
    /// Node fractional coordinates.
    pub nodes: Vec<[f64; 3]>,
    /// Edges with shifts.
    pub edges: Vec<NetEdge>,
    /// Declared node labels from the CGD file.
    pub labels: Vec<String>,
}

/// Simple invariants used by net search and identification.
#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
pub struct NetInvariants {
    /// Number of nodes.
    pub nodes: usize,
    /// Number of edges.
    pub edges: usize,
    /// Cycle rank `E - V + 1` for a connected quotient graph.
    pub cycle_rank: usize,
    /// Coordination sequence in node order.
    pub coordination: Vec<usize>,
}

/// Builds the quotient graph of a CGD entry.
///
/// Endpoints are matched to nodes modulo an integer translation within
/// `tolerance`; an unmatched endpoint is an error.
///
/// # Errors
///
/// Returns [`Error::InvalidInput`] for an empty net, a bad tolerance or an
/// endpoint that does not match any node.
pub fn quotient_graph(entry: &CgdEntry, tolerance: f64) -> Result<QuotientGraph, Error> {
    if entry.nodes.is_empty() {
        return Err(Error::InvalidInput(
            "a net needs at least one node".to_owned(),
        ));
    }
    if !tolerance.is_finite() || tolerance <= 0.0 {
        return Err(Error::InvalidInput(
            "net tolerance must be positive".to_owned(),
        ));
    }
    let dimension = entry.nodes[0].frac.len();
    let mut nodes = Vec::with_capacity(entry.nodes.len());
    let mut labels = Vec::with_capacity(entry.nodes.len());
    for node in &entry.nodes {
        if node.frac.len() != dimension {
            return Err(Error::InvalidInput(
                "net nodes disagree on dimension".to_owned(),
            ));
        }
        nodes.push(pad(&node.frac));
        labels.push(node.label.clone());
    }
    let mut edges = Vec::with_capacity(entry.edges.len());
    for (first, second) in &entry.edges {
        let (a, _) = match_node(&nodes, first, tolerance)?;
        let (b, shift) = match_node(&nodes, second, tolerance)?;
        edges.push(NetEdge { a, b, shift });
    }
    edges.sort();
    Ok(QuotientGraph {
        nodes,
        edges,
        labels,
    })
}

impl QuotientGraph {
    /// Coordination numbers and counts.
    #[must_use]
    pub fn invariants(&self) -> NetInvariants {
        let mut coordination = vec![0usize; self.nodes.len()];
        for edge in &self.edges {
            if let Some(value) = coordination.get_mut(edge.a) {
                *value += 1;
            }
            if let Some(value) = coordination.get_mut(edge.b) {
                *value += 1;
            }
        }
        let cycle_rank = self
            .edges
            .len()
            .saturating_sub(self.nodes.len())
            .saturating_add(1);
        NetInvariants {
            nodes: self.nodes.len(),
            edges: self.edges.len(),
            cycle_rank,
            coordination,
        }
    }

    /// Lattice shifts per edge, keyed by node pair.
    #[must_use]
    pub fn shift_map(&self) -> BTreeMap<(usize, usize), Vec<[i32; 3]>> {
        let mut map: BTreeMap<(usize, usize), Vec<[i32; 3]>> = BTreeMap::new();
        for edge in &self.edges {
            map.entry((edge.a, edge.b)).or_default().push(edge.shift);
        }
        map
    }
}

fn match_node(
    nodes: &[[f64; 3]],
    point: &[f64],
    tolerance: f64,
) -> Result<(usize, [i32; 3]), Error> {
    let point = pad(point);
    for (index, node) in nodes.iter().enumerate() {
        let mut shift = [0i32; 3];
        let mut matched = true;
        for axis in 0..3 {
            let delta = point[axis] - node[axis];
            let rounded = delta.round();
            if (delta - rounded).abs() > tolerance {
                matched = false;
                break;
            }
            shift[axis] = rounded as i32;
        }
        if matched {
            return Ok((index, shift));
        }
    }
    Err(Error::InvalidInput(
        "a net endpoint does not match any node".to_owned(),
    ))
}

fn pad(point: &[f64]) -> [f64; 3] {
    [
        point.first().copied().unwrap_or(0.0),
        point.get(1).copied().unwrap_or(0.0),
        point.get(2).copied().unwrap_or(0.0),
    ]
}
