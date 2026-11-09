//! Ring perception on the periodic bond graph (Task 1.12).
//!
//! The search returns primitive (shortest-path) rings of Franzblau: a cycle is
//! kept only if no two of its atoms are closer through the graph than they are
//! round the ring. Rings carry the lattice image of every atom, so a ring that
//! closes through a cell face is still one ring and a small cell can show the
//! same atom more than once.

mod search;

use std::collections::BTreeSet;

use crate::{BondGraph, Error, Structure};

/// One atom of a ring path, with its lattice image relative to the first atom.
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub struct RingStep {
    /// Atom index in the P1 cell.
    pub atom: usize,
    /// Lattice image applied to this atom in the path.
    pub image: [i32; 3],
}

/// One primitive ring.
#[derive(Debug, Clone, PartialEq, Eq, PartialOrd, Ord)]
pub struct Ring {
    /// Unique atom indices in the ring, sorted.
    pub atoms: Vec<usize>,
    /// Atoms in ring order, each with its lattice image.
    pub path: Vec<RingStep>,
}

impl Ring {
    /// Number of atoms in the ring path (a small cell can repeat an atom).
    #[must_use]
    pub fn size(&self) -> usize {
        self.path.len()
    }
}

/// Finds primitive rings up to `max_size` atoms.
///
/// # Errors
///
/// Returns [`Error::InvalidInput`] for `max_size < 3`, a failed P1 expansion,
/// or a bond graph too dense to search within the Franzblau budget.
pub fn rings(structure: &Structure, max_size: usize) -> Result<Vec<Ring>, Error> {
    if max_size < 3 {
        return Err(Error::InvalidInput(
            "ring size limit must be at least 3".to_owned(),
        ));
    }
    search::find(structure, max_size)
}

/// Finds primitive rings using only the stored explicit bonds.
///
/// Builders state their bonds; this is the census those bonds describe,
/// without distance perception adding a close contact across a pore.
///
/// # Errors
///
/// Returns [`Error::InvalidInput`] for `max_size < 3`, a non-P1 structure
/// or a graph too dense to search.
pub fn rings_explicit(structure: &Structure, max_size: usize) -> Result<Vec<Ring>, Error> {
    if max_size < 3 {
        return Err(Error::InvalidInput(
            "ring size limit must be at least 3".to_owned(),
        ));
    }
    let graph = BondGraph::from_explicit(structure)?;
    search::find_graph(&graph, max_size)
}

/// Counts explicit-bond rings by size.
///
/// # Errors
///
/// Returns the errors of [`rings_explicit`].
pub fn ring_counts_explicit(
    structure: &Structure,
    max_size: usize,
) -> Result<Vec<(usize, usize)>, Error> {
    let mut counts: Vec<(usize, usize)> = Vec::new();
    for ring in rings_explicit(structure, max_size)? {
        if let Some(entry) = counts.iter_mut().find(|(size, _)| *size == ring.size()) {
            entry.1 += 1;
        } else {
            counts.push((ring.size(), 1));
        }
    }
    Ok(counts)
}

/// Ring paths with the lattice image of every atom, in ring order.
///
/// This is the periodic form of [`rings`]; each path can be laid out in
/// Cartesian space by adding `image` to the atom's fractional position.
///
/// # Errors
///
/// Returns the errors of [`rings`].
pub fn ring_paths(structure: &Structure, max_size: usize) -> Result<Vec<Vec<RingStep>>, Error> {
    Ok(rings(structure, max_size)?
        .into_iter()
        .map(|ring| ring.path)
        .collect())
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

/// Number of rings each P1 atom belongs to (`ring_membership`).
///
/// An atom repeated in a periodic ring is counted once for that ring.
///
/// # Errors
///
/// Returns the errors of [`rings`].
pub fn ring_membership(structure: &Structure, max_size: usize) -> Result<Vec<usize>, Error> {
    let p1 = structure.p1()?;
    let mut membership = vec![0usize; p1.atoms.len()];
    for ring in rings(structure, max_size)? {
        for atom in ring.atoms {
            if let Some(slot) = membership.get_mut(atom) {
                *slot += 1;
            }
        }
    }
    Ok(membership)
}

/// Heuristic aromaticity for a ring: five- or six-membered and all C/N/O.
///
/// A real aromaticity check needs bond orders and planarity; this flag is a
/// first approximation for the molecule builder and is labelled as such.
///
/// # Errors
///
/// Returns the errors of the P1 expansion.
pub fn is_aromatic(structure: &Structure, ring: &Ring) -> Result<bool, Error> {
    if !matches!(ring.size(), 5 | 6) {
        return Ok(false);
    }
    let p1 = structure.p1()?;
    Ok(ring.atoms.iter().all(|index| {
        p1.atoms
            .get(*index)
            .is_some_and(|atom| matches!(atom.element.as_str(), "C" | "N" | "O"))
    }))
}

/// Canonical ring key: translate the path so its least step is in the home
/// cell, then forget the direction. Used by the search to list each ring once.
pub(crate) fn translated_home(path: &[RingStep]) -> BTreeSet<RingStep> {
    let Some(minimum) = path.iter().min() else {
        return BTreeSet::new();
    };
    path.iter()
        .map(|step| RingStep {
            atom: step.atom,
            image: [
                step.image[0] - minimum.image[0],
                step.image[1] - minimum.image[1],
                step.image[2] - minimum.image[2],
            ],
        })
        .collect()
}
