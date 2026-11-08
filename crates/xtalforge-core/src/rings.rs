//! Ring perception on the periodic bond graph (Task 1.12).
//!
//! The search returns chordless simple cycles up to a size limit: a cycle with
//! a chord is a sum of smaller rings, not a primitive ring. This is the
//! Franzblau primitive-ring criterion for the molecular graph; the aromaticity
//! flag and the full framework statistics follow.

use crate::{Error, Structure};

/// One ring as a sorted list of P1 atom indices.
#[derive(Debug, Clone, PartialEq, Eq, PartialOrd, Ord)]
pub struct Ring {
    /// Atom indices in the ring.
    pub atoms: Vec<usize>,
}

impl Ring {
    /// Ring size.
    #[must_use]
    pub fn size(&self) -> usize {
        self.atoms.len()
    }
}

/// Finds chordless rings up to `max_size` atoms.
///
/// # Errors
///
/// Returns [`Error::InvalidInput`] for `max_size < 3` or a failed expansion.
pub fn rings(structure: &Structure, max_size: usize) -> Result<Vec<Ring>, Error> {
    if max_size < 3 {
        return Err(Error::InvalidInput(
            "ring size limit must be at least 3".to_owned(),
        ));
    }
    let graph = structure.bond_graph()?;
    let mut found: Vec<Ring> = Vec::new();
    for start in 0..graph.n_atoms() {
        let mut path = vec![start];
        let mut visited = vec![false; graph.n_atoms()];
        if let Some(slot) = visited.get_mut(start) {
            *slot = true;
        }
        search(
            &graph,
            start,
            start,
            &mut path,
            &mut visited,
            max_size,
            &mut found,
        );
    }
    found.sort();
    found.dedup();
    Ok(found)
}

/// Counts rings by size.
///
/// # Errors
///
/// Returns the errors of [`rings`].
pub fn ring_counts(structure: &Structure, max_size: usize) -> Result<Vec<(usize, usize)>, Error> {
    let mut counts: Vec<(usize, usize)> = Vec::new();
    for ring in rings(structure, max_size)? {
        if let Some(entry) = counts.iter_mut().find(|(size, _)| *size == ring.size()) {
            entry.1 += 1;
        } else {
            counts.push((ring.size(), 1));
        }
    }
    Ok(counts)
}

fn search(
    graph: &crate::BondGraph,
    start: usize,
    current: usize,
    path: &mut Vec<usize>,
    visited: &mut [bool],
    max_size: usize,
    found: &mut Vec<Ring>,
) {
    if path.len() > max_size {
        return;
    }
    for bond_index in graph.bonds_of(current) {
        let Some(bond) = graph.bonds.get(*bond_index) else {
            continue;
        };
        let neighbour = if bond.i == current { bond.j } else { bond.i };
        if neighbour == start {
            if path.len() >= 3 {
                let mut atoms = path.clone();
                atoms.sort_unstable();
                if is_chordless(graph, &atoms) {
                    found.push(Ring { atoms });
                }
            }
            continue;
        }
        if neighbour <= start || visited.get(neighbour).copied().unwrap_or(true) {
            continue;
        }
        if let Some(slot) = visited.get_mut(neighbour) {
            *slot = true;
        }
        path.push(neighbour);
        search(graph, start, neighbour, path, visited, max_size, found);
        path.pop();
        if let Some(slot) = visited.get_mut(neighbour) {
            *slot = false;
        }
    }
}

fn is_chordless(graph: &crate::BondGraph, atoms: &[usize]) -> bool {
    for (position, atom) in atoms.iter().enumerate() {
        for (offset, other) in atoms.iter().enumerate() {
            if offset <= position + 1 || (position == 0 && offset + 1 == atoms.len()) {
                continue;
            }
            if adjacent(graph, *atom, *other) {
                return false;
            }
        }
    }
    true
}

fn adjacent(graph: &crate::BondGraph, first: usize, second: usize) -> bool {
    graph.bonds_of(first).iter().any(|bond_index| {
        graph.bonds.get(*bond_index).is_some_and(|bond| {
            (bond.i == first && bond.j == second) || (bond.i == second && bond.j == first)
        })
    })
}
