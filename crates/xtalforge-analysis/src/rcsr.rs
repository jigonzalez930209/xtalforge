//! RCSR identification against a small curated net subset (Tasks 5.6.1,
//! 5.6.4 and 5.6.5).
//!
//! The full RCSR database (about 1.9 MB, `RCSRnets-2019-06-01.cgd.gz`) is
//! **not** bundled with this crate. It is an external data bundle, as recorded
//! in [`RCSR_LICENSE_NOTE`]. What ships here is `data/rcsr_subset.json`: 17
//! common nets whose primitive quotient graphs were derived from the RCSR
//! 2019-06-01 release. Identification uses only the node count, edge count,
//! cycle rank and the multiset of coordination numbers. Nets that share those
//! four invariants (for example the 2D `kgm` and the 3D `nbo` quotient graphs)
//! cannot be separated by counts alone, so the subset keeps one representative
//! per signature and the full catalogue resolves the others by canonical key.

use std::collections::BTreeSet;

use serde::{Deserialize, Serialize};
use xtalforge_core::Error;

use crate::nets::{NetEdge, NetInvariants, QuotientGraph};

/// Licence and provenance note for RCSR data.
///
/// The RCSR (Reticular Chemistry Structure Resource, <https://rcsr.ou.edu>)
/// is a curated database of periodic nets. Only the small derived subset in
/// `data/rcsr_subset.json` travels with this crate; the complete 2019-06-01
/// snapshot is an external bundle that users obtain under the RCSR terms of
/// use. Redistributing the full snapshot requires the permission of the RCSR
/// maintainers and must be recorded in `THIRD_PARTY.md`.
pub const RCSR_LICENSE_NOTE: &str = "RCSR data (Reticular Chemistry Structure \
Resource, rcsr.ou.edu) are used as a small derived subset of net names and \
primitive quotient graphs. The full 2019-06-01 database is not bundled with \
xtalforge; it is an external data bundle distributed under the RCSR terms of \
use. Verify those terms before redistributing the full snapshot.";

/// The embedded curated subset, exactly as shipped.
const RCSR_SUBSET: &str = include_str!("../data/rcsr_subset.json");

/// One named net in the subset index.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct RcsrEntry {
    /// RCSR name (`pcu`, `dia`, `sql`, ...).
    pub name: String,
    /// Primitive quotient graph as stored in the subset.
    pub graph: QuotientGraph,
    /// Provenance of the entry.
    pub source: String,
}

/// A small searchable RCSR subset.
#[derive(Debug, Clone, Default, PartialEq, Serialize, Deserialize)]
pub struct RcsrIndex {
    /// Entries in file order (alphabetical by name).
    pub entries: Vec<RcsrEntry>,
}

/// One row of `data/rcsr_subset.json`.
#[derive(Debug, Deserialize)]
struct SubsetNet {
    name: String,
    nodes: usize,
    edges: Vec<[usize; 2]>,
    source: String,
}

impl RcsrIndex {
    /// Loads the curated subset embedded in the crate.
    ///
    /// The data file is compiled into the binary. If it is ever edited into
    /// invalid JSON this returns an empty index; call [`RcsrIndex::validate`]
    /// in tests to catch that before shipping.
    #[must_use]
    pub fn from_subset() -> Self {
        let nets: Vec<SubsetNet> = serde_json::from_str(RCSR_SUBSET).unwrap_or_default();
        let entries = nets
            .into_iter()
            .map(|net| {
                let graph = subset_graph(&net);
                RcsrEntry {
                    name: net.name,
                    graph,
                    source: net.source,
                }
            })
            .collect();
        Self { entries }
    }

    /// Names the subset entry whose invariants match `graph`.
    ///
    /// Matching is deliberately simple: node count, edge count, cycle rank and
    /// the sorted coordination numbers. It identifies every entry of the
    /// curated subset, whose signatures are unique (checked by
    /// [`RcsrIndex::validate`]); the full RCSR catalogue resolves the remaining
    /// nets by canonical key.
    #[must_use]
    pub fn identify(&self, graph: &QuotientGraph) -> Option<&str> {
        let invariants = graph.invariants();
        let mut coordination = invariants.coordination.clone();
        coordination.sort_unstable();
        self.entries
            .iter()
            .find(|entry| {
                let candidate = entry.graph.invariants();
                let mut candidate_coordination = candidate.coordination;
                candidate_coordination.sort_unstable();
                candidate.nodes == invariants.nodes
                    && candidate.edges == invariants.edges
                    && candidate.cycle_rank == invariants.cycle_rank
                    && candidate_coordination == coordination
            })
            .map(|entry| entry.name.as_str())
    }

    /// Checks every entry and returns how many are valid.
    ///
    /// A valid entry has a non-empty name, at least one node and one edge,
    /// finite node coordinates, one label per node, endpoints inside the node
    /// range, and a connected quotient graph. Names and invariant signatures
    /// must also be unique so that [`RcsrIndex::identify`] is unambiguous.
    ///
    /// # Errors
    ///
    /// Returns [`Error::InvalidInput`] naming the first entry that fails.
    pub fn validate(&self) -> Result<usize, Error> {
        if self.entries.is_empty() {
            return Err(Error::InvalidInput("the RCSR subset is empty".to_owned()));
        }
        let mut names = BTreeSet::new();
        let mut signatures = BTreeSet::new();
        for entry in &self.entries {
            let graph = &entry.graph;
            if entry.name.trim().is_empty() {
                return Err(Error::InvalidInput(
                    "an RCSR subset entry has an empty name".to_owned(),
                ));
            }
            if graph.nodes.is_empty() {
                return Err(invalid(&entry.name, "it has no nodes"));
            }
            if graph.labels.len() != graph.nodes.len() {
                return Err(invalid(
                    &entry.name,
                    "its label count differs from its nodes",
                ));
            }
            if graph.edges.is_empty() {
                return Err(invalid(&entry.name, "it has no edges"));
            }
            if !graph.nodes.iter().flatten().all(|value| value.is_finite()) {
                return Err(invalid(&entry.name, "it has a non-finite node coordinate"));
            }
            for edge in &graph.edges {
                if edge.a >= graph.nodes.len() || edge.b >= graph.nodes.len() {
                    return Err(invalid(
                        &entry.name,
                        "an edge endpoint is outside the graph",
                    ));
                }
            }
            if graph.components().len() != 1 {
                return Err(invalid(&entry.name, "it is not connected"));
            }
            if !names.insert(entry.name.as_str()) {
                return Err(invalid(&entry.name, "the name is repeated"));
            }
            let invariants = graph.invariants();
            let mut coordination = invariants.coordination;
            coordination.sort_unstable();
            let signature = (
                invariants.nodes,
                invariants.edges,
                invariants.cycle_rank,
                coordination,
            );
            if !signatures.insert(signature) {
                return Err(invalid(
                    &entry.name,
                    "another entry has the same invariant signature",
                ));
            }
        }
        Ok(self.entries.len())
    }

    /// The subset's invariants, in file order.
    #[must_use]
    pub fn invariant_entries(&self) -> Vec<(&str, NetInvariants)> {
        self.entries
            .iter()
            .map(|entry| (entry.name.as_str(), entry.graph.invariants()))
            .collect()
    }
}

fn subset_graph(net: &SubsetNet) -> QuotientGraph {
    let count = net.nodes.max(1);
    let nodes = (0..net.nodes)
        .map(|index| [index as f64 / count as f64, 0.0, 0.0])
        .collect();
    let mut edges: Vec<NetEdge> = net
        .edges
        .iter()
        .map(|pair| {
            let [a, b] = *pair;
            NetEdge {
                a,
                b,
                shift: [0; 3],
            }
        })
        .collect();
    edges.sort();
    let labels = (0..net.nodes)
        .map(|index| format!("{}{}", net.name, index + 1))
        .collect();
    QuotientGraph {
        nodes,
        edges,
        labels,
    }
}

fn invalid(name: &str, reason: &str) -> Error {
    Error::InvalidInput(format!("RCSR subset entry {name} is invalid: {reason}"))
}
