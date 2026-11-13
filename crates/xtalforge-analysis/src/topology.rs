//! Net fingerprints for RCSR identification (Task 5.6.3 and 5.6.4).
//!
//! A fingerprint is the periodicity and the sets of coordination sequences
//! and Schläfli point symbols of the net's distinct vertices. Both invariants
//! are computed on the infinite net, so a supercell or a different cell
//! setting cannot change them.
//!
//! The functions here port Crystal-Builder `xtal/core/bonding.py` (MIT):
//! the periodic breadth-first walk counts `(vertex, lattice shift)` pairs, not
//! atoms, and the point symbol searches the shortest cycle through every pair
//! of edges at a vertex.

#![allow(clippy::indexing_slicing)] // bounded local edge and index arrays

use std::collections::{BTreeMap, BTreeSet, HashSet};

use serde::{Deserialize, Serialize};

use crate::{Error, NetEdge, QuotientGraph};

/// Coordination-sequence depth used by the RCSR index.
pub const DEPTH: usize = 10;
/// Largest ring the point symbol searches for.
pub const MAX_RING: usize = 12;

/// The invariants that name a net.
#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
pub struct NetFingerprint {
    /// Periodicity dimension (0, 1, 2 or 3).
    pub periodicity: usize,
    /// Distinct coordination sequences, sorted.
    pub sequences: Vec<Vec<usize>>,
    /// Distinct point symbols, sorted.
    pub symbols: Vec<String>,
}

impl NetFingerprint {
    /// The catalogue token (`periodicity|walks|symbols`).
    #[must_use]
    pub fn token(&self) -> String {
        format!(
            "{}|{}|{}",
            self.periodicity,
            self.sequences
                .iter()
                .map(|sequence| {
                    sequence
                        .iter()
                        .map(usize::to_string)
                        .collect::<Vec<_>>()
                        .join(",")
                })
                .collect::<Vec<_>>()
                .join(";"),
            self.symbols.join(";")
        )
    }

    /// The token without the point symbols, used for nearest candidates.
    #[must_use]
    pub fn sequence_token(&self) -> String {
        format!(
            "{}|{}",
            self.periodicity,
            self.sequences
                .iter()
                .map(|sequence| {
                    sequence
                        .iter()
                        .map(usize::to_string)
                        .collect::<Vec<_>>()
                        .join(",")
                })
                .collect::<Vec<_>>()
                .join(";")
        )
    }

    /// FNV-1a 64-bit hash of [`NetFingerprint::token`].
    #[must_use]
    pub fn hash(&self) -> u64 {
        fnv1a(self.token().as_bytes())
    }

    /// Salted 32-bit FNV-1a hash, the embedded index's lookup key.
    ///
    /// The generator picks the salt so that two different tokens cannot
    /// share a hash; the index header carries it.
    #[must_use]
    pub fn hash_with_salt(&self, salt: u32) -> u32 {
        let mut hash = 0xcbf2_9ce4_8422_2325_u64;
        for byte in salt.to_le_bytes().iter().chain(self.token().as_bytes()) {
            hash ^= u64::from(*byte);
            hash = hash.wrapping_mul(0x0000_0100_0000_01b3);
        }
        (hash ^ (hash >> 32)) as u32
    }

    /// FNV-1a 64-bit hash of [`NetFingerprint::sequence_token`].
    #[must_use]
    pub fn sequence_hash(&self) -> u64 {
        fnv1a(self.sequence_token().as_bytes())
    }

    /// Distinct coordination numbers, sorted.
    #[must_use]
    pub fn coordination(&self) -> Vec<usize> {
        let numbers: BTreeSet<usize> = self
            .sequences
            .iter()
            .filter_map(|sequence| sequence.first().copied())
            .collect();
        numbers.into_iter().collect()
    }
}

/// FNV-1a, the stable hash written into the embedded index.
#[must_use]
pub fn fnv1a(bytes: &[u8]) -> u64 {
    let mut hash = 0xcbf2_9ce4_8422_2325_u64;
    for byte in bytes {
        hash ^= u64::from(*byte);
        hash = hash.wrapping_mul(0x0000_0100_0000_01b3);
    }
    hash
}

/// Adjacency with the lattice shift of each neighbour.
#[must_use]
pub fn adjacency(graph: &QuotientGraph) -> Vec<Vec<(usize, [i32; 3])>> {
    let mut out = vec![Vec::new(); graph.nodes.len()];
    for edge in &graph.edges {
        if let Some(list) = out.get_mut(edge.a) {
            list.push((edge.b, edge.shift));
        }
        if let Some(list) = out.get_mut(edge.b) {
            list.push((edge.a, edge.shift.map(|value| -value)));
        }
    }
    for list in &mut out {
        list.sort_unstable();
    }
    out
}

/// How many vertices lie exactly `k` edges from `vertex`.
#[must_use]
pub fn coordination_sequence(
    adjacency: &[Vec<(usize, [i32; 3])>],
    vertex: usize,
    depth: usize,
) -> Vec<usize> {
    let origin = (vertex, [0, 0, 0]);
    let mut seen: HashSet<(usize, [i32; 3])> = HashSet::from([origin]);
    let mut frontier: HashSet<(usize, [i32; 3])> = HashSet::from([origin]);
    let mut out = Vec::with_capacity(depth);
    for _ in 0..depth {
        let mut next = HashSet::new();
        for (node, offset) in &frontier {
            let Some(neighbours) = adjacency.get(*node) else {
                continue;
            };
            for (other, shift) in neighbours {
                let step = (
                    *other,
                    [
                        offset[0] + shift[0],
                        offset[1] + shift[1],
                        offset[2] + shift[2],
                    ],
                );
                if seen.insert(step) {
                    next.insert(step);
                }
            }
        }
        out.push(next.len());
        if next.is_empty() {
            break;
        }
        frontier = next;
    }
    out
}

/// The Schläfli point symbol at one vertex.
#[must_use]
pub fn point_symbol(
    adjacency: &[Vec<(usize, [i32; 3])>],
    vertex: usize,
    max_ring: usize,
) -> String {
    let Some(partners) = adjacency.get(vertex) else {
        return String::new();
    };
    let origin = (vertex, [0, 0, 0]);
    let mut counts: BTreeMap<(u8, usize), usize> = BTreeMap::new();
    for first in 0..partners.len() {
        for second in (first + 1)..partners.len() {
            let size = smallest_ring(
                adjacency,
                origin,
                partners[first],
                partners[second],
                max_ring,
            );
            let key = if size == 0 { (1, 0) } else { (0, size) };
            *counts.entry(key).or_default() += 1;
        }
    }
    counts
        .into_iter()
        .map(|((star, size), count)| {
            let name = if star == 0 {
                size.to_string()
            } else {
                "*".to_owned()
            };
            if count > 1 {
                format!("{name}^{count}")
            } else {
                name
            }
        })
        .collect::<Vec<_>>()
        .join(".")
}

pub(crate) fn smallest_ring(
    adjacency: &[Vec<(usize, [i32; 3])>],
    origin: (usize, [i32; 3]),
    first: (usize, [i32; 3]),
    second: (usize, [i32; 3]),
    max_ring: usize,
) -> usize {
    if first == second {
        return 0;
    }
    let mut seen: HashSet<(usize, [i32; 3])> = HashSet::from([origin, first]);
    let mut frontier = vec![first];
    for steps in 1..max_ring.saturating_sub(1) {
        let mut next = Vec::new();
        for (node, offset) in &frontier {
            let Some(neighbours) = adjacency.get(*node) else {
                continue;
            };
            for (other, shift) in neighbours {
                let step = (
                    *other,
                    [
                        offset[0] + shift[0],
                        offset[1] + shift[1],
                        offset[2] + shift[2],
                    ],
                );
                if step == second {
                    return steps + 2;
                }
                if step == origin || !seen.insert(step) {
                    continue;
                }
                next.push(step);
            }
        }
        if next.is_empty() {
            break;
        }
        frontier = next;
    }
    0
}

/// Rank of the lattice the net's own cycles generate.
#[must_use]
pub fn periodicity(graph: &QuotientGraph) -> usize {
    let adjacency = adjacency(graph);
    let mut offset: Vec<Option<[i64; 3]>> = vec![None; graph.nodes.len()];
    let mut voltages: Vec<[i64; 3]> = Vec::new();
    for start in 0..graph.nodes.len() {
        if offset.get(start).is_some_and(Option::is_some) {
            continue;
        }
        if let Some(slot) = offset.get_mut(start) {
            *slot = Some([0, 0, 0]);
        }
        let mut stack = vec![start];
        while let Some(node) = stack.pop() {
            let base = offset.get(node).and_then(|value| *value).unwrap_or([0; 3]);
            let Some(neighbours) = adjacency.get(node) else {
                continue;
            };
            for (other, shift) in neighbours {
                let image = [
                    base[0] + i64::from(shift[0]),
                    base[1] + i64::from(shift[1]),
                    base[2] + i64::from(shift[2]),
                ];
                if offset.get(*other).is_some_and(Option::is_some) {
                    let existing = offset
                        .get(*other)
                        .and_then(|value| *value)
                        .unwrap_or([0; 3]);
                    if existing != image {
                        voltages.push([
                            image[0] - existing[0],
                            image[1] - existing[1],
                            image[2] - existing[2],
                        ]);
                    }
                } else {
                    if let Some(slot) = offset.get_mut(*other) {
                        *slot = Some(image);
                    }
                    stack.push(*other);
                }
            }
        }
    }
    integer_rank(&voltages)
}

fn integer_rank(vectors: &[[i64; 3]]) -> usize {
    let mut rows: Vec<[i128; 3]> = vectors
        .iter()
        .map(|vector| vector.map(i128::from))
        .collect();
    let mut rank = 0usize;
    for column in 0..3 {
        let pivot = (rank..rows.len()).find(|row| rows[*row][column] != 0);
        let Some(pivot) = pivot else {
            continue;
        };
        rows.swap(rank, pivot);
        for row in (rank + 1)..rows.len() {
            if rows[row][column] != 0 {
                let pivot_row = rows[rank];
                let factor = rows[row][column];
                let lead = pivot_row[column];
                rows[row] = [
                    lead * rows[row][0] - factor * pivot_row[0],
                    lead * rows[row][1] - factor * pivot_row[1],
                    lead * rows[row][2] - factor * pivot_row[2],
                ];
            }
        }
        rank += 1;
    }
    rank
}

/// Computes the fingerprint of a quotient graph.
///
/// # Errors
///
/// Returns [`Error::InvalidInput`] for a graph with no nodes.
pub fn fingerprint(graph: &QuotientGraph) -> Result<NetFingerprint, Error> {
    if graph.nodes.is_empty() {
        return Err(Error::InvalidInput(
            "a net fingerprint needs at least one node".to_owned(),
        ));
    }
    let adjacency = adjacency(graph);
    let mut sequences = BTreeSet::new();
    let mut symbols = BTreeSet::new();
    for vertex in 0..graph.nodes.len() {
        sequences.insert(coordination_sequence(&adjacency, vertex, DEPTH));
        symbols.insert(point_symbol(&adjacency, vertex, MAX_RING));
    }
    Ok(NetFingerprint {
        periodicity: periodicity(graph),
        sequences: sequences.into_iter().collect(),
        symbols: symbols.into_iter().collect(),
    })
}

impl QuotientGraph {
    /// Parses one canonical key (`d:n:a,b,sx,sy,sz;...`).
    ///
    /// The key stores the minimal cell the RCSR index was built from. Two
    /// dimensional keys get a zero `z` shift.
    ///
    /// # Errors
    ///
    /// Returns [`Error::InvalidInput`] for a malformed key or an endpoint
    /// outside the node count.
    pub fn from_canonical_key(key: &str) -> Result<Self, Error> {
        let mut parts = key.split(':');
        let dimension = parts
            .next()
            .and_then(|value| value.parse::<usize>().ok())
            .filter(|value| *value == 2 || *value == 3)
            .ok_or_else(|| Error::InvalidInput(format!("invalid net key {key:?}")))?;
        let nodes = parts
            .next()
            .and_then(|value| value.parse::<usize>().ok())
            .ok_or_else(|| Error::InvalidInput(format!("invalid net key {key:?}")))?;
        let body = parts.next().unwrap_or("");
        let mut edges = Vec::new();
        for chunk in body.split(';') {
            if chunk.is_empty() {
                continue;
            }
            let values: Vec<i64> = chunk
                .split(',')
                .filter_map(|value| value.parse().ok())
                .collect();
            let needed = dimension + 2;
            if values.len() != needed {
                return Err(Error::InvalidInput(format!(
                    "net key edge {chunk:?} does not have {needed} numbers"
                )));
            }
            let a = usize::try_from(values[0])
                .map_err(|_| Error::InvalidInput("negative net key node".to_owned()))?;
            let b = usize::try_from(values[1])
                .map_err(|_| Error::InvalidInput("negative net key node".to_owned()))?;
            if a >= nodes || b >= nodes {
                return Err(Error::InvalidInput(format!(
                    "net key edge {chunk:?} leaves the node range"
                )));
            }
            let shift = if dimension == 3 {
                [
                    i32::try_from(values[2]).unwrap_or(0),
                    i32::try_from(values[3]).unwrap_or(0),
                    i32::try_from(values[4]).unwrap_or(0),
                ]
            } else {
                [
                    i32::try_from(values[2]).unwrap_or(0),
                    i32::try_from(values[3]).unwrap_or(0),
                    0,
                ]
            };
            edges.push(NetEdge { a, b, shift });
        }
        Ok(Self {
            nodes: vec![[0.0; 3]; nodes],
            edges,
            labels: (0..nodes).map(|index| format!("N{index}")).collect(),
        })
    }
}
