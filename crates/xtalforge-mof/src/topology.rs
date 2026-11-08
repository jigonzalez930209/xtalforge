//! The curated topology table (Tasks 6.11 and 6.12).
//!
//! Each topology is a finite periodic quotient: nodes are fractional
//! coordinates in the cell, and edges are undirected index pairs. The cell is
//! `[a, b, c, alpha, beta, gamma]` in Angstrom and degrees; when it is absent
//! the builder uses a 20 Angstrom cubic box.

use std::collections::BTreeSet;

use serde::{Deserialize, Serialize};
use xtalforge_core::Lattice;

use crate::Error;

/// One net in the curated topology table.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct Topology {
    /// Net name, for example `pcu`.
    pub name: String,
    /// Fractional node coordinates in the unit cell.
    pub nodes: Vec<[f64; 3]>,
    /// Undirected edges as node index pairs.
    pub edges: Vec<[usize; 2]>,
    /// Cell parameters `[a, b, c, alpha, beta, gamma]`; `None` means a 20
    /// Angstrom cubic default.
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub cell: Option<[f64; 6]>,
    /// Provenance of the curated entry.
    pub source: String,
}

impl Topology {
    /// Number of edges incident to one vertex.
    #[must_use]
    pub fn degree(&self, vertex: usize) -> usize {
        self.edges
            .iter()
            .filter(|edge| {
                let [first, second] = **edge;
                first == vertex || second == vertex
            })
            .count()
    }

    /// The cell of this topology, or a 20 Angstrom cubic default.
    ///
    /// # Errors
    ///
    /// Returns [`Error::Core`] when the stored cell parameters are invalid.
    pub fn lattice(&self) -> Result<Lattice, Error> {
        match self.cell {
            Some([a, b, c, alpha, beta, gamma]) => {
                Ok(Lattice::from_parameters(a, b, c, alpha, beta, gamma)?)
            }
            None => Ok(Lattice::cubic(20.0)?),
        }
    }

    /// Checks names, coordinates, indices, duplicate edges and the cell.
    ///
    /// # Errors
    ///
    /// Returns [`Error::InvalidInput`] for an empty name or node list,
    /// non-finite coordinates, an out-of-range or repeated index, a duplicate
    /// edge, or a cell that [`Lattice`] rejects.
    pub fn validate(&self) -> Result<(), Error> {
        if self.name.trim().is_empty() {
            return Err(Error::InvalidInput(
                "topology name must not be empty".to_owned(),
            ));
        }
        if self.nodes.is_empty() {
            return Err(Error::InvalidInput(format!(
                "topology `{}` has no nodes",
                self.name
            )));
        }
        if self
            .nodes
            .iter()
            .flatten()
            .any(|coordinate| !coordinate.is_finite())
        {
            return Err(Error::InvalidInput(format!(
                "topology `{}` has a non-finite node coordinate",
                self.name
            )));
        }
        let mut seen = BTreeSet::new();
        for edge in &self.edges {
            let [first, second] = *edge;
            if first == second {
                return Err(Error::InvalidInput(format!(
                    "topology `{}` has a self-loop at node {first}",
                    self.name
                )));
            }
            let limit = self.nodes.len();
            if first >= limit || second >= limit {
                return Err(Error::InvalidInput(format!(
                    "topology `{}` edge {first}-{second} is outside 0..{limit}",
                    self.name
                )));
            }
            let key = if first < second {
                (first, second)
            } else {
                (second, first)
            };
            if !seen.insert(key) {
                return Err(Error::InvalidInput(format!(
                    "topology `{}` repeats edge {first}-{second}",
                    self.name
                )));
            }
        }
        self.lattice()?;
        Ok(())
    }
}

/// All embedded topologies, validated.
///
/// # Errors
///
/// Returns [`Error::Data`] when the embedded JSON is invalid, has a duplicate
/// name, or fails [`Topology::validate`].
pub fn topologies() -> Result<Vec<Topology>, Error> {
    let entries: Vec<Topology> = serde_json::from_str(include_str!("../data/topologies.json"))
        .map_err(|error| Error::Data(format!("topology data: {error}")))?;
    let mut names = BTreeSet::new();
    for topology in &entries {
        topology
            .validate()
            .map_err(|error| Error::Data(format!("topology data: {error}")))?;
        if !names.insert(topology.name.to_ascii_lowercase()) {
            return Err(Error::Data(format!(
                "topology data repeats name `{}`",
                topology.name
            )));
        }
    }
    Ok(entries)
}

/// Finds a topology by case-insensitive name.
///
/// # Errors
///
/// Returns the errors of [`topologies`].
pub fn find_topology(name: &str) -> Result<Option<Topology>, Error> {
    Ok(topologies()?
        .into_iter()
        .find(|topology| topology.name.eq_ignore_ascii_case(name)))
}
