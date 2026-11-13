//! Vertex symbols: the shortest rings through every pair of edges at a node
//! (Task 14.4.2).
//!
//! A point symbol ([`crate::topology::fingerprint`]) aggregates the ring sizes
//! around a vertex into one Schläfli string. A vertex symbol keeps the
//! per-pair answer, which is what distinguishes nets that walk and ring
//! differently but aggregate to the same counts.
//!
//! The ring search is the same one the point symbol uses: a periodic
//! breadth-first walk through `(node, lattice shift)` pairs, so a ring that
//! closes across a cell boundary counts once. The ring search is the one the
//! point symbol uses, shared with it inside the crate.
//!
//! # Deliberate subset
//!
//! The conventional (O'Keeffe) vertex symbol lists the pairs of *adjacent*
//! edges in the cyclic order of the coordination figure around the vertex.
//! A quotient graph carries no angular order, so this module lists every
//! unordered pair of incident edges in the order `(i, j)` with `i < j`.
//! The full-ring property is preserved; only the printed order and the
//! omission of the diagonals differ from the conventional symbol when the
//! geometry is known. Callers that need the conventional cyclic form can
//! reorder [`VertexSymbol::pairs`] with their own geometry.

#![allow(clippy::indexing_slicing)] // bounded local index arrays

use serde::{Deserialize, Serialize};
use xtalforge_core::Error;

use crate::nets::QuotientGraph;
use crate::topology::{adjacency, smallest_ring};

/// The shortest rings around one node of a net.
#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
pub struct VertexSymbol {
    /// Node index in the quotient graph.
    pub node: usize,
    /// Neighbour node for each incident edge, in adjacency order.
    pub neighbours: Vec<usize>,
    /// Lattice shift of each incident edge, in adjacency order.
    pub shifts: Vec<[i32; 3]>,
    /// Shortest ring through each incident edge; `0` means none within
    /// `max_ring`.
    pub edge_rings: Vec<usize>,
    /// For each unordered pair of incident edges, the shortest ring through
    /// both (row-major over `i < j`); `None` means none within `max_ring`.
    pub pairs: Vec<Option<usize>>,
    /// The smallest ring at this node, `None` when there is none.
    pub shortest: Option<usize>,
}

impl VertexSymbol {
    /// Ring size of the pair `(first, second)`, if the pair exists.
    #[must_use]
    pub fn pair_ring(&self, first: usize, second: usize) -> Option<usize> {
        if first >= second || second >= self.neighbours.len() {
            return None;
        }
        let offset = first * self.neighbours.len() - first * (first + 1) / 2 + (second - first - 1);
        self.pairs.get(offset).copied().flatten()
    }

    /// The per-pair symbol: ring sizes in pair order, `*` for no ring, joined
    /// by `.` (for example `4.4.4.4.4.4` for **pcu**).
    #[must_use]
    pub fn text(&self) -> String {
        if self.pairs.is_empty() {
            return String::new();
        }
        self.pairs
            .iter()
            .map(|ring| ring.map_or_else(|| "*".to_owned(), |size| size.to_string()))
            .collect::<Vec<_>>()
            .join(".")
    }
}

impl std::fmt::Display for VertexSymbol {
    fn fmt(&self, formatter: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        formatter.write_str(&self.text())
    }
}

/// Computes the vertex symbol of every node.
///
/// `max_ring` bounds the ring search; the point-symbol convention uses
/// [`crate::topology::MAX_RING`]. The graph must have at least one node.
///
/// # Errors
///
/// Returns [`Error::InvalidInput`] for an empty graph or a `max_ring` below
/// three.
pub fn vertex_symbols(graph: &QuotientGraph, max_ring: usize) -> Result<Vec<VertexSymbol>, Error> {
    if graph.nodes.is_empty() {
        return Err(Error::InvalidInput(
            "a vertex symbol needs at least one node".to_owned(),
        ));
    }
    if max_ring < 3 {
        return Err(Error::InvalidInput(
            "a ring needs at least three edges".to_owned(),
        ));
    }
    let adjacency = adjacency(graph);
    let mut symbols = Vec::with_capacity(graph.nodes.len());
    for node in 0..graph.nodes.len() {
        symbols.push(symbol_of(&adjacency, node, max_ring));
    }
    Ok(symbols)
}

/// The vertex symbols as text, one per node in node order.
///
/// # Errors
///
/// Returns the errors of [`vertex_symbols`].
pub fn vertex_symbol_text(graph: &QuotientGraph, max_ring: usize) -> Result<Vec<String>, Error> {
    Ok(vertex_symbols(graph, max_ring)?
        .iter()
        .map(VertexSymbol::text)
        .collect())
}

fn symbol_of(adjacency: &[Vec<(usize, [i32; 3])>], node: usize, max_ring: usize) -> VertexSymbol {
    let partners = adjacency.get(node).cloned().unwrap_or_default();
    let origin = (node, [0, 0, 0]);
    let count = partners.len();
    let mut pairs = Vec::with_capacity(count.saturating_mul(count.saturating_sub(1)) / 2);
    let mut edge_rings = vec![0usize; count];
    let mut shortest: Option<usize> = None;
    for first in 0..count {
        for second in (first + 1)..count {
            let ring = smallest_ring(
                adjacency,
                origin,
                partners[first],
                partners[second],
                max_ring,
            );
            if ring != 0 {
                for index in [first, second] {
                    if let Some(slot) = edge_rings.get_mut(index)
                        && (*slot == 0 || ring < *slot)
                    {
                        *slot = ring;
                    }
                }
                shortest = Some(shortest.map_or(ring, |value| value.min(ring)));
                pairs.push(Some(ring));
            } else {
                pairs.push(None);
            }
        }
    }
    VertexSymbol {
        node,
        neighbours: partners.iter().map(|(other, _)| *other).collect(),
        shifts: partners.iter().map(|(_, shift)| *shift).collect(),
        edge_rings,
        pairs,
        shortest,
    }
}
