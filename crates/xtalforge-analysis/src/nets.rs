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

    /// The subgraph induced by `nodes`, renumbered in the given order.
    ///
    /// Edges with an endpoint outside `nodes` are dropped. Repeated nodes
    /// make the result empty rather than silently duplicating a vertex.
    ///
    /// # Errors
    ///
    /// Returns [`Error::InvalidInput`] for a node index outside the graph or
    /// a repeated node.
    pub fn subgraph(&self, nodes: &[usize]) -> Result<Self, Error> {
        let mut renumber = BTreeMap::new();
        for (position, node) in nodes.iter().enumerate() {
            if *node >= self.nodes.len() || renumber.insert(*node, position).is_some() {
                return Err(Error::InvalidInput(
                    "a subgraph node is outside the graph or repeated".to_owned(),
                ));
            }
        }
        let mut sub_nodes = Vec::with_capacity(nodes.len());
        let mut labels = Vec::with_capacity(nodes.len());
        for node in nodes {
            let frac = self.nodes.get(*node).ok_or_else(|| {
                Error::InvalidInput("a subgraph node is outside the graph".to_owned())
            })?;
            sub_nodes.push(*frac);
            labels.push(
                self.labels
                    .get(*node)
                    .cloned()
                    .unwrap_or_else(|| format!("N{node}")),
            );
        }
        let mut edges = Vec::new();
        for edge in &self.edges {
            let (Some(&a), Some(&b)) = (renumber.get(&edge.a), renumber.get(&edge.b)) else {
                continue;
            };
            edges.push(NetEdge {
                a,
                b,
                shift: edge.shift,
            });
        }
        edges.sort();
        Ok(Self {
            nodes: sub_nodes,
            edges,
            labels,
        })
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

/// Builds a searchable index of quotient graphs from a CGD file.
#[derive(Debug, Clone, Default, PartialEq, Serialize, Deserialize)]
pub struct NetIndex {
    /// One invariant set per successfully read net.
    pub entries: Vec<(String, NetInvariants)>,
}

impl NetIndex {
    /// Indexes every entry that has a valid quotient graph.
    #[must_use]
    pub fn from_file(file: &xtalforge_io::cgd::CgdFile, tolerance: f64) -> Self {
        let mut entries = Vec::new();
        for entry in &file.entries {
            if let Ok(graph) = quotient_graph(entry, tolerance) {
                entries.push((entry.name.clone(), graph.invariants()));
            }
        }
        Self { entries }
    }

    /// Nets with the same node/edge counts, cycle rank and coordination.
    #[must_use]
    pub fn search(&self, invariants: &NetInvariants) -> Vec<&str> {
        self.entries
            .iter()
            .filter(|(_, candidate)| {
                candidate.nodes == invariants.nodes
                    && candidate.edges == invariants.edges
                    && candidate.cycle_rank == invariants.cycle_rank
                    && candidate.coordination == invariants.coordination
            })
            .map(|(name, _)| name.as_str())
            .collect()
    }
}

impl QuotientGraph {
    /// Connected components of the quotient graph, ignoring lattice shifts.
    #[must_use]
    pub fn components(&self) -> Vec<Vec<usize>> {
        let mut adjacency = vec![Vec::new(); self.nodes.len()];
        for edge in &self.edges {
            if let Some(list) = adjacency.get_mut(edge.a) {
                list.push(edge.b);
            }
            if let Some(list) = adjacency.get_mut(edge.b) {
                list.push(edge.a);
            }
        }
        let mut seen = vec![false; self.nodes.len()];
        let mut components = Vec::new();
        for start in 0..self.nodes.len() {
            if seen.get(start).copied().unwrap_or(true) {
                continue;
            }
            let mut stack = vec![start];
            let mut component = Vec::new();
            while let Some(node) = stack.pop() {
                if seen.get(node).copied().unwrap_or(true) {
                    continue;
                }
                if let Some(slot) = seen.get_mut(node) {
                    *slot = true;
                }
                component.push(node);
                if let Some(neighbours) = adjacency.get(node) {
                    stack.extend(neighbours.iter().copied());
                }
            }
            component.sort_unstable();
            components.push(component);
        }
        components
    }
}

/// Duplicates a quotient graph at a lattice shift, producing two interpenetrating
/// components. The result keeps the two copies as disconnected components.
///
/// # Errors
///
/// Returns [`Error::InvalidInput`] for a non-finite shift.
pub fn interpenetrate(graph: &QuotientGraph, shift: [i32; 3]) -> Result<QuotientGraph, Error> {
    let count = graph.nodes.len();
    let mut nodes = graph.nodes.clone();
    let mut labels = graph.labels.clone();
    for (index, node) in graph.nodes.iter().enumerate() {
        nodes.push([
            (node[0] + f64::from(shift[0])).rem_euclid(1.0),
            (node[1] + f64::from(shift[1])).rem_euclid(1.0),
            (node[2] + f64::from(shift[2])).rem_euclid(1.0),
        ]);
        labels.push(format!(
            "{}{}",
            graph.labels.get(index).cloned().unwrap_or_default(),
            "'"
        ));
    }
    let mut edges = graph.edges.clone();
    for edge in &graph.edges {
        edges.push(NetEdge {
            a: edge.a + count,
            b: edge.b + count,
            shift: edge.shift,
        });
    }
    edges.sort();
    Ok(QuotientGraph {
        nodes,
        edges,
        labels,
    })
}
