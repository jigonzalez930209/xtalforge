//! The P1 view the force field evaluates: real atoms only, and their periodic
//! bond geometry.
//!
//! Dummy `X` atoms are filtered here once; every later stage indexes the
//! compact list and [`AtomView::to_full`] maps forces back. A bond graph over
//! the full P1 cell supplies the periodic images, so a bond crossing a cell
//! boundary is stored as one bond with its translation, exactly as in
//! Crystal-Builder `xtal/ff/uff/calculator.py`.

use std::collections::BTreeMap;

use serde_json::Value;
use xtalforge_core::{BondGraph, P1Cell, Structure};

use crate::Error;
use crate::uff::math::{add, norm, shift_vector, sub};

/// One perceived bond mapped into the compact atom view.
#[derive(Debug, Clone, Copy, PartialEq)]
pub(crate) struct LocalBond {
    pub i: usize,
    pub j: usize,
    pub image: [i32; 3],
    pub order: f64,
}

/// A compact, dummy-free view of a P1 cell.
#[derive(Debug, Clone)]
pub(crate) struct AtomView {
    pub elements: Vec<String>,
    pub positions: Vec<[f64; 3]>,
    pub to_full: Vec<usize>,
    pub overrides: Vec<Option<String>>,
    pub matrix: [[f64; 3]; 3],
    pub n_full: usize,
    full_to_local: Vec<Option<usize>>,
}

/// The chemical element a P1 atom should be typed as.
///
/// Deuterium shares hydrogen's UFF surface; the reference does the same in
/// `xtal.ff.isotopes.as_elements`.
pub(crate) fn chemical_element(symbol: &str) -> String {
    match symbol {
        "D" => "H".to_owned(),
        other => other.to_owned(),
    }
}

fn site_override(structure: &Structure, site_index: usize) -> Option<String> {
    structure
        .sites()
        .get(site_index)
        .and_then(|site| site.props().get("uff_type"))
        .and_then(Value::as_str)
        .map(str::to_owned)
}

impl AtomView {
    pub(crate) fn from_p1(structure: &Structure, p1: &P1Cell) -> Result<Self, Error> {
        let mut elements = Vec::new();
        let mut positions = Vec::new();
        let mut to_full = Vec::new();
        let mut overrides = Vec::new();
        let mut full_to_local = vec![None; p1.atoms.len()];
        for (full, atom) in p1.atoms.iter().enumerate() {
            if atom.element == "X" {
                continue;
            }
            if let Some(slot) = full_to_local.get_mut(full) {
                *slot = Some(elements.len());
            }
            elements.push(chemical_element(&atom.element));
            positions.push(p1.lattice.to_cart(atom.frac));
            to_full.push(full);
            overrides.push(site_override(structure, atom.site_idx));
        }
        if elements.is_empty() {
            return Err(Error::InvalidInput(
                "there are no atoms to compute an energy for".to_owned(),
            ));
        }
        Ok(Self {
            elements,
            positions,
            to_full,
            overrides,
            matrix: p1.lattice.matrix(),
            n_full: p1.atoms.len(),
            full_to_local,
        })
    }

    pub(crate) fn n_atoms(&self) -> usize {
        self.elements.len()
    }

    pub(crate) fn local(&self, full: usize) -> Option<usize> {
        self.full_to_local.get(full).copied().flatten()
    }

    /// Cartesian vector from local atom `i` to the image of `j`.
    pub(crate) fn vector(&self, i: usize, j: usize, image: [i32; 3]) -> Result<[f64; 3], Error> {
        let start = self
            .positions
            .get(i)
            .ok_or_else(|| Error::Topology(format!("atom index {i} is out of range")))?;
        let end = self
            .positions
            .get(j)
            .ok_or_else(|| Error::Topology(format!("atom index {j} is out of range")))?;
        Ok(sub(add(*end, shift_vector(image, self.matrix)), *start))
    }

    /// Local forces scattered back onto the full P1 atom list, zero for `X`.
    pub(crate) fn full_forces(&self, local: &[[f64; 3]]) -> Vec<[f64; 3]> {
        let mut full = vec![[0.0; 3]; self.n_full];
        for (index, force) in local.iter().enumerate() {
            if let Some(target) = self
                .to_full
                .get(index)
                .copied()
                .and_then(|full_index| full.get_mut(full_index))
            {
                *target = *force;
            }
        }
        full
    }
}

/// A neighbour reached through one bond and its lattice image.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub(crate) struct Neighbor {
    pub index: usize,
    pub shift: [i32; 3],
    pub bond: usize,
}

/// Bond adjacency in the compact view.
#[derive(Debug, Clone, Default)]
pub(crate) struct Geometry {
    pub neighbors: Vec<Vec<Neighbor>>,
}

impl Geometry {
    pub(crate) fn new(n_atoms: usize, bonds: &[LocalBond]) -> Self {
        let mut neighbors = vec![Vec::new(); n_atoms];
        for (index, bond) in bonds.iter().enumerate() {
            if let Some(list) = neighbors.get_mut(bond.i) {
                list.push(Neighbor {
                    index: bond.j,
                    shift: bond.image,
                    bond: index,
                });
            }
            if let Some(list) = neighbors.get_mut(bond.j) {
                list.push(Neighbor {
                    index: bond.i,
                    shift: bond.image.map(|value| -value),
                    bond: index,
                });
            }
        }
        Self { neighbors }
    }

    pub(crate) fn partners(&self, index: usize) -> &[Neighbor] {
        self.neighbors.get(index).map_or(&[], Vec::as_slice)
    }

    pub(crate) fn coordination(&self, index: usize) -> usize {
        self.partners(index).len()
    }

    pub(crate) fn distances(&self, view: &AtomView, index: usize) -> Vec<f64> {
        self.partners(index)
            .iter()
            .filter_map(|neighbor| {
                view.vector(index, neighbor.index, neighbor.shift)
                    .ok()
                    .map(norm)
            })
            .collect()
    }

    /// Every neighbour vector of one atom, relative to that atom.
    pub(crate) fn vectors(&self, view: &AtomView, index: usize) -> Vec<[f64; 3]> {
        self.partners(index)
            .iter()
            .filter_map(|neighbor| view.vector(index, neighbor.index, neighbor.shift).ok())
            .collect()
    }

    /// Largest angle between two neighbours, in degrees; zero when fewer
    /// than two neighbours exist.
    pub(crate) fn max_angle(&self, view: &AtomView, index: usize) -> f64 {
        let vectors = self.vectors(view, index);
        let mut best: f64 = 0.0;
        for (first, vector_a) in vectors.iter().enumerate() {
            for vector_b in vectors.iter().skip(first + 1) {
                best = best.max(crate::uff::math::angle_between(*vector_a, *vector_b));
            }
        }
        best.to_degrees()
    }

    /// Sum of all neighbour-neighbour angles at one atom, in degrees.
    pub(crate) fn angle_sum(&self, view: &AtomView, index: usize) -> f64 {
        let vectors = self.vectors(view, index);
        let mut sum = 0.0;
        for (first, vector_a) in vectors.iter().enumerate() {
            for vector_b in vectors.iter().skip(first + 1) {
                sum += crate::uff::math::angle_between(*vector_a, *vector_b).to_degrees();
            }
        }
        sum
    }
}

/// Remaps a bond graph over the full P1 cell onto the compact view.
pub(crate) fn local_bonds(graph: &BondGraph, view: &AtomView) -> Vec<LocalBond> {
    let mut bonds = Vec::new();
    for bond in &graph.bonds {
        let (Some(i), Some(j)) = (view.local(bond.i), view.local(bond.j)) else {
            continue;
        };
        bonds.push(LocalBond {
            i,
            j,
            image: bond.image,
            order: bond.order,
        });
    }
    bonds
}

/// Rings mapped from full P1 indices, dropping any that touch a dummy atom.
///
/// Ring perception is the core's default-rules search (Task 1.12); a custom
/// [`xtalforge_core::BondRules`] in [`crate::CalculatorOptions`] changes the
/// bond graph and the terms but not which rings are called aromatic. That is
/// the documented subset of the reference, where the core graph and the ring
/// search share one rule set.
pub(crate) fn local_rings(
    structure: &Structure,
    view: &AtomView,
) -> Result<Vec<Vec<usize>>, Error> {
    let mut rings = Vec::new();
    for ring in xtalforge_core::rings(structure, 6)? {
        let mapped: Option<Vec<usize>> =
            ring.atoms.iter().map(|index| view.local(*index)).collect();
        if let Some(mapped) = mapped {
            rings.push(mapped);
        }
    }
    Ok(rings)
}

/// Per-bond translations keyed the same way as the neighbour search.
///
/// Lower atom index first; for a self-pair the lexicographically positive
/// translation is kept. This mirrors `_pair_key` in the reference calculator.
pub(crate) type PairKey = (usize, usize, [i32; 3]);

pub(crate) fn pair_key(i: usize, j: usize, shift: [i32; 3]) -> PairKey {
    if i < j {
        return (i, j, shift);
    }
    if i > j {
        return (j, i, shift.map(|value| -value));
    }
    let flipped = shift.map(|value| -value);
    (i, j, if shift >= flipped { shift } else { flipped })
}

/// Exclusions for the van der Waals pair list: every bonded 1-2 pair and
/// every 1-3 pair implied by an angle.
pub(crate) fn exclusions(bonds: &[LocalBond], geometry: &Geometry) -> BTreeMap<PairKey, ()> {
    let mut excluded = BTreeMap::new();
    for bond in bonds {
        excluded.insert(pair_key(bond.i, bond.j, bond.image), ());
    }
    for partners in &geometry.neighbors {
        for (first, neighbor_a) in partners.iter().enumerate() {
            for neighbor_b in partners.iter().skip(first + 1) {
                let shift = [
                    neighbor_b.shift[0] - neighbor_a.shift[0],
                    neighbor_b.shift[1] - neighbor_a.shift[1],
                    neighbor_b.shift[2] - neighbor_a.shift[2],
                ];
                excluded.insert(pair_key(neighbor_a.index, neighbor_b.index, shift), ());
            }
        }
    }
    excluded
}
