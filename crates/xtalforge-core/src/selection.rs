//! Selection model and atom tracking (Task 3.11).

use std::collections::{BTreeMap, BTreeSet};

use serde::{Deserialize, Serialize};

/// A set of atom indices and bond index pairs.
#[derive(Debug, Clone, Default, PartialEq, Eq, Serialize, Deserialize)]
pub struct Selection {
    /// Selected atoms in the current cell.
    pub atoms: BTreeSet<usize>,
    /// Selected bonds as ordered index pairs.
    pub bonds: BTreeSet<(usize, usize)>,
}

impl Selection {
    /// Empty selection.
    #[must_use]
    pub fn new() -> Self {
        Self::default()
    }

    /// Selects one atom.
    pub fn add_atom(&mut self, index: usize) {
        self.atoms.insert(index);
    }

    /// Deselects one atom.
    pub fn remove_atom(&mut self, index: usize) {
        self.atoms.remove(&index);
    }

    /// Toggles one atom.
    pub fn toggle_atom(&mut self, index: usize) {
        if !self.atoms.insert(index) {
            self.atoms.remove(&index);
        }
    }

    /// Selects a sorted block of atoms.
    pub fn add_atoms(&mut self, indices: impl IntoIterator<Item = usize>) {
        self.atoms.extend(indices);
    }

    /// Selects one bond; the pair is normalised to `(min, max)`.
    pub fn add_bond(&mut self, first: usize, second: usize) {
        self.bonds.insert((first.min(second), first.max(second)));
    }

    /// Whether an atom is selected.
    #[must_use]
    pub fn contains(&self, index: usize) -> bool {
        self.atoms.contains(&index)
    }

    /// Number of selected atoms.
    #[must_use]
    pub fn len(&self) -> usize {
        self.atoms.len()
    }

    /// Whether nothing is selected.
    #[must_use]
    pub fn is_empty(&self) -> bool {
        self.atoms.is_empty() && self.bonds.is_empty()
    }

    /// Clears atoms and bonds.
    pub fn clear(&mut self) {
        self.atoms.clear();
        self.bonds.clear();
    }

    /// Union with another selection.
    pub fn union(&mut self, other: &Self) {
        self.atoms.extend(other.atoms.iter().copied());
        self.bonds.extend(other.bonds.iter().copied());
    }

    /// Intersection with another selection.
    pub fn intersect(&mut self, other: &Self) {
        self.atoms.retain(|index| other.atoms.contains(index));
        self.bonds.retain(|bond| other.bonds.contains(bond));
    }
}

/// A stable atom identifier that survives command edits.
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash, Serialize, Deserialize)]
pub struct AtomId(pub u64);

/// Maps hero atom ids to current site indices.
#[derive(Debug, Clone, Default, PartialEq, Eq, Serialize, Deserialize)]
pub struct TrackingTable {
    /// Current site index by id.
    pub by_id: BTreeMap<AtomId, usize>,
    next: u64,
}

impl TrackingTable {
    /// Empty table.
    #[must_use]
    pub fn new() -> Self {
        Self::default()
    }

    /// Assigns a new id to a site and returns it.
    pub fn assign(&mut self, site: usize) -> AtomId {
        let id = AtomId(self.next);
        self.next += 1;
        self.by_id.insert(id, site);
        id
    }

    /// Current site of an id.
    #[must_use]
    pub fn site_of(&self, id: AtomId) -> Option<usize> {
        self.by_id.get(&id).copied()
    }

    /// Id assigned to a site, if any.
    #[must_use]
    pub fn id_of(&self, site: usize) -> Option<AtomId> {
        self.by_id
            .iter()
            .find(|(_, index)| **index == site)
            .map(|(id, _)| *id)
    }

    /// Removes an id.
    pub fn remove(&mut self, id: AtomId) {
        self.by_id.remove(&id);
    }

    /// Rewrites every site index through a mapping; missing entries are
    /// removed.
    pub fn remap(&mut self, mapping: &BTreeMap<usize, usize>) {
        self.by_id.retain(|_, site| {
            if let Some(new) = mapping.get(site) {
                *site = *new;
                true
            } else {
                false
            }
        });
    }
}
