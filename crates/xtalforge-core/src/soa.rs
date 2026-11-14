//! Structure-of-arrays views for P1 cells and neighbor lists (Task 15.2.1).
//!
//! The existing [`crate::P1Cell`] and [`crate::Neighbor`] types stay the
//! readable, allocation-per-atom API. This module adds additive views that
//! store the numeric fields in separate contiguous columns, which is what a
//! million-atom pipeline wants: the coordinates of one axis are a single
//! slice, and a neighbor list groups every neighbor of one atom in one
//! contiguous row range.
//!
//! The views are exact: [`P1Arrays::to_p1`] rebuilds the same [`P1Cell`], and
//! [`NeighborArrays::iter`] yields the same [`Neighbor`] values in the same
//! order as the source list.

use std::ops::Range;

use crate::neighbors::{Neighbor, search_columns};
use crate::p1::{P1Atom, P1Cell};
use crate::{Error, Lattice, Structure};

/// Structure-of-arrays view of a P1 cell.
///
/// Every accessor returns the column for one field. The columns always have
/// the same length; a view built from a [`P1Cell`] or a [`Structure`] is
/// self-consistent by construction.
#[derive(Debug, Clone, PartialEq)]
pub struct P1Arrays {
    lattice: Lattice,
    elements: Vec<String>,
    labels: Vec<String>,
    x: Vec<f64>,
    y: Vec<f64>,
    z: Vec<f64>,
    occupancy: Vec<f64>,
    site: Vec<usize>,
    operation: Vec<usize>,
    tau: Vec<[i32; 3]>,
    disorder: Vec<Option<u32>>,
}

impl P1Arrays {
    /// Copies a [`P1Cell`] into the column layout.
    #[must_use]
    pub fn from_p1(cell: &P1Cell) -> Self {
        let mut view = Self {
            lattice: cell.lattice,
            elements: Vec::with_capacity(cell.atoms.len()),
            labels: Vec::with_capacity(cell.atoms.len()),
            x: Vec::with_capacity(cell.atoms.len()),
            y: Vec::with_capacity(cell.atoms.len()),
            z: Vec::with_capacity(cell.atoms.len()),
            occupancy: Vec::with_capacity(cell.atoms.len()),
            site: Vec::with_capacity(cell.atoms.len()),
            operation: Vec::with_capacity(cell.atoms.len()),
            tau: Vec::with_capacity(cell.atoms.len()),
            disorder: Vec::with_capacity(cell.atoms.len()),
        };
        for atom in &cell.atoms {
            view.push(atom);
        }
        view
    }

    /// Expands a structure and copies the result into the column layout.
    ///
    /// # Errors
    ///
    /// Returns the errors of [`P1Cell::expand`].
    pub fn from_structure(structure: &Structure, tolerance: f64) -> Result<Self, Error> {
        Ok(Self::from_p1(&P1Cell::expand(structure, tolerance)?))
    }

    /// Expands a structure with the default P1 tolerance.
    ///
    /// # Errors
    ///
    /// Returns the errors of [`P1Cell::expand`].
    pub fn from_structure_default(structure: &Structure) -> Result<Self, Error> {
        Ok(Self::from_p1(&structure.p1()?))
    }

    /// Rebuilds the equivalent [`P1Cell`].
    #[must_use]
    pub fn to_p1(&self) -> P1Cell {
        let mut atoms = Vec::with_capacity(self.len());
        for index in 0..self.len() {
            let (Some(element), Some(label), Some(frac)) = (
                self.elements.get(index),
                self.labels.get(index),
                self.frac(index),
            ) else {
                continue;
            };
            atoms.push(P1Atom {
                element: element.clone(),
                label: label.clone(),
                frac,
                occupancy: self.occupancy.get(index).copied().unwrap_or(1.0),
                site_idx: self.site.get(index).copied().unwrap_or(0),
                op_idx: self.operation.get(index).copied().unwrap_or(0),
                tau: self.tau.get(index).copied().unwrap_or([0, 0, 0]),
                disorder_group: self.disorder.get(index).copied().flatten(),
            });
        }
        P1Cell {
            lattice: self.lattice,
            atoms,
        }
    }

    /// Cell lattice.
    #[must_use]
    pub fn lattice(&self) -> &Lattice {
        &self.lattice
    }

    /// Number of atoms.
    #[must_use]
    pub fn len(&self) -> usize {
        self.x.len()
    }

    /// Whether the cell has no atoms.
    #[must_use]
    pub fn is_empty(&self) -> bool {
        self.x.is_empty()
    }

    /// Element symbols in atom order.
    #[must_use]
    pub fn elements(&self) -> &[String] {
        &self.elements
    }

    /// Site labels in atom order.
    #[must_use]
    pub fn labels(&self) -> &[String] {
        &self.labels
    }

    /// Fractional `x` column.
    #[must_use]
    pub fn x(&self) -> &[f64] {
        &self.x
    }

    /// Fractional `y` column.
    #[must_use]
    pub fn y(&self) -> &[f64] {
        &self.y
    }

    /// Fractional `z` column.
    #[must_use]
    pub fn z(&self) -> &[f64] {
        &self.z
    }

    /// All three fractional columns.
    #[must_use]
    pub fn columns(&self) -> (&[f64], &[f64], &[f64]) {
        (&self.x, &self.y, &self.z)
    }

    /// Fractional coordinates of one atom.
    #[must_use]
    pub fn frac(&self, index: usize) -> Option<[f64; 3]> {
        Some([
            *self.x.get(index)?,
            *self.y.get(index)?,
            *self.z.get(index)?,
        ])
    }

    /// Occupancies in atom order.
    #[must_use]
    pub fn occupancy(&self) -> &[f64] {
        &self.occupancy
    }

    /// Parent asymmetric-site indices in atom order.
    #[must_use]
    pub fn site_indices(&self) -> &[usize] {
        &self.site
    }

    /// Operation indices in atom order.
    #[must_use]
    pub fn operation_indices(&self) -> &[usize] {
        &self.operation
    }

    /// Integer lattice shifts in atom order.
    #[must_use]
    pub fn taus(&self) -> &[[i32; 3]] {
        &self.tau
    }

    /// Disorder groups in atom order.
    #[must_use]
    pub fn disorder_groups(&self) -> &[Option<u32>] {
        &self.disorder
    }

    /// Builds the same neighbor list [`crate::neighbors()`] returns, without
    /// expanding the cell again.
    ///
    /// # Errors
    ///
    /// Returns [`Error::InvalidInput`] for a non-positive or non-finite
    /// cutoff.
    pub fn neighbors(&self, cutoff: f64) -> Result<Vec<Neighbor>, Error> {
        search_columns(&self.lattice, &self.x, &self.y, &self.z, cutoff)
    }

    /// Neighbor list in the column layout.
    ///
    /// # Errors
    ///
    /// Returns the errors of [`P1Arrays::neighbors`] or
    /// [`NeighborArrays::new`].
    pub fn neighbor_arrays(&self, cutoff: f64) -> Result<NeighborArrays, Error> {
        NeighborArrays::new(&self.neighbors(cutoff)?, self.len())
    }

    fn push(&mut self, atom: &P1Atom) {
        self.elements.push(atom.element.clone());
        self.labels.push(atom.label.clone());
        self.x.push(atom.frac[0]);
        self.y.push(atom.frac[1]);
        self.z.push(atom.frac[2]);
        self.occupancy.push(atom.occupancy);
        self.site.push(atom.site_idx);
        self.operation.push(atom.op_idx);
        self.tau.push(atom.tau);
        self.disorder.push(atom.disorder_group);
    }
}

/// Structure-of-arrays view of a neighbor list.
///
/// Rows are grouped by the first atom index `i`: every neighbor of atom `k`
/// occupies `row_range(k)`, a contiguous `Range<usize>` into the columns.
/// Values inside a row keep the order of the source list.
#[derive(Debug, Clone, PartialEq)]
pub struct NeighborArrays {
    n_atoms: usize,
    i: Vec<usize>,
    j: Vec<usize>,
    image_x: Vec<i32>,
    image_y: Vec<i32>,
    image_z: Vec<i32>,
    distance: Vec<f64>,
    starts: Vec<usize>,
}

impl NeighborArrays {
    /// Copies a neighbor list into the column layout and groups it by `i`.
    ///
    /// The input order inside one row is preserved, so a list returned by
    /// [`crate::neighbors()`] keeps its deterministic `(i, j, image)` order.
    ///
    /// # Errors
    ///
    /// Returns [`Error::InvalidInput`] when an atom index is outside
    /// `0..n_atoms`.
    pub fn new(neighbors: &[Neighbor], n_atoms: usize) -> Result<Self, Error> {
        let mut counts = vec![0_usize; n_atoms];
        for neighbor in neighbors {
            if neighbor.i >= n_atoms || neighbor.j >= n_atoms {
                return Err(Error::InvalidInput(format!(
                    "neighbor ({}, {}) is outside a cell with {n_atoms} atoms",
                    neighbor.i, neighbor.j
                )));
            }
            if let Some(count) = counts.get_mut(neighbor.i) {
                *count += 1;
            }
        }
        let mut starts = Vec::with_capacity(n_atoms.saturating_add(1));
        starts.push(0);
        for count in &counts {
            let next = starts
                .last()
                .copied()
                .unwrap_or(0_usize)
                .saturating_add(*count);
            starts.push(next);
        }
        let mut cursor = starts.clone();
        let mut view = Self {
            n_atoms,
            i: vec![0; neighbors.len()],
            j: vec![0; neighbors.len()],
            image_x: vec![0; neighbors.len()],
            image_y: vec![0; neighbors.len()],
            image_z: vec![0; neighbors.len()],
            distance: vec![0.0; neighbors.len()],
            starts,
        };
        for neighbor in neighbors {
            let Some(position) = cursor.get_mut(neighbor.i) else {
                continue;
            };
            let slot = *position;
            *position += 1;
            if let (Some(i), Some(j), Some(x), Some(y), Some(z), Some(distance)) = (
                view.i.get_mut(slot),
                view.j.get_mut(slot),
                view.image_x.get_mut(slot),
                view.image_y.get_mut(slot),
                view.image_z.get_mut(slot),
                view.distance.get_mut(slot),
            ) {
                *i = neighbor.i;
                *j = neighbor.j;
                *x = neighbor.image[0];
                *y = neighbor.image[1];
                *z = neighbor.image[2];
                *distance = neighbor.distance;
            }
        }
        Ok(view)
    }

    /// Number of atoms the rows are indexed by.
    #[must_use]
    pub fn n_atoms(&self) -> usize {
        self.n_atoms
    }

    /// Number of neighbor entries.
    #[must_use]
    pub fn len(&self) -> usize {
        self.i.len()
    }

    /// Whether the list has no entries.
    #[must_use]
    pub fn is_empty(&self) -> bool {
        self.i.is_empty()
    }

    /// First-atom column.
    #[must_use]
    pub fn i(&self) -> &[usize] {
        &self.i
    }

    /// Second-atom column.
    #[must_use]
    pub fn j(&self) -> &[usize] {
        &self.j
    }

    /// Lattice-image columns of `j`.
    #[must_use]
    pub fn images(&self) -> (&[i32], &[i32], &[i32]) {
        (&self.image_x, &self.image_y, &self.image_z)
    }

    /// Distances in Angstrom.
    #[must_use]
    pub fn distances(&self) -> &[f64] {
        &self.distance
    }

    /// One neighbor entry.
    #[must_use]
    pub fn get(&self, index: usize) -> Option<Neighbor> {
        Some(Neighbor {
            i: *self.i.get(index)?,
            j: *self.j.get(index)?,
            image: [
                *self.image_x.get(index)?,
                *self.image_y.get(index)?,
                *self.image_z.get(index)?,
            ],
            distance: *self.distance.get(index)?,
        })
    }

    /// Contiguous column range of every neighbor of `atom`.
    #[must_use]
    pub fn row_range(&self, atom: usize) -> Option<Range<usize>> {
        let start = *self.starts.get(atom)?;
        let end = *self.starts.get(atom.saturating_add(1))?;
        Some(start..end)
    }

    /// Entries in file order, which is the source order.
    pub fn iter(&self) -> impl Iterator<Item = Neighbor> + '_ {
        (0..self.len()).filter_map(|index| self.get(index))
    }
}

/// Neighbor list in the column layout.
///
/// # Errors
///
/// Returns the errors of [`crate::neighbors()`] or [`NeighborArrays::new`].
pub fn neighbors_arrays(structure: &Structure, cutoff: f64) -> Result<NeighborArrays, Error> {
    P1Arrays::from_structure_default(structure)?.neighbor_arrays(cutoff)
}
