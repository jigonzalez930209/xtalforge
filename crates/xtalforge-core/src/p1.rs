//! P1 expansion with provenance back to the asymmetric site and operation.

use std::collections::HashMap;

use serde::{Deserialize, Serialize};

use crate::Error;
use crate::Lattice;
use crate::Structure;
use crate::symop::SymOp;

/// One atom of an expanded P1 cell.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct P1Atom {
    /// Element symbol.
    pub element: String,
    /// Site label inherited from the asymmetric unit.
    pub label: String,
    /// Fractional coordinates in the cell.
    pub frac: [f64; 3],
    /// Site occupancy.
    pub occupancy: f64,
    /// Index of the parent site.
    pub site_idx: usize,
    /// Index of the operation that generated this atom.
    pub op_idx: usize,
    /// Integer lattice shift applied after the operation.
    pub tau: [i32; 3],
}

/// The expanded P1 cell.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct P1Cell {
    /// Cell lattice.
    pub lattice: Lattice,
    /// Expanded atoms.
    pub atoms: Vec<P1Atom>,
}

impl P1Cell {
    /// Expands a structure with its explicit operations, or the identity when
    /// it has none.
    ///
    /// # Errors
    ///
    /// Returns [`Error::InvalidInput`] for a negative or non-finite tolerance.
    pub fn expand(structure: &Structure, tolerance: f64) -> Result<Self, Error> {
        if !tolerance.is_finite() || tolerance < 0.0 {
            return Err(Error::InvalidInput(
                "P1 expansion tolerance must be non-negative and finite".to_owned(),
            ));
        }
        let identity = [SymOp::identity()];
        let operations = if structure.operations().is_empty() {
            identity.as_slice()
        } else {
            structure.operations()
        };
        let reciprocal = structure.lattice().reciprocal(false)?;
        let widths = reciprocal
            .matrix()
            .map(|row| 1.0 / row.iter().map(|value| value * value).sum::<f64>().sqrt());
        let counts = if tolerance > 0.0 {
            widths.map(|width| {
                let bins = (width / tolerance).floor();
                if bins >= 1.0 { bins as usize } else { 1 }
            })
        } else {
            [1, 1, 1]
        };
        let [count_x, count_y, count_z] = counts;
        let offsets_x = offsets(count_x);
        let offsets_y = offsets(count_y);
        let offsets_z = offsets(count_z);
        let mut grid: HashMap<(&str, [usize; 3]), Vec<usize>> = HashMap::new();
        let mut atoms: Vec<P1Atom> = Vec::new();
        for (site_idx, site) in structure.sites().iter().enumerate() {
            for (op_idx, operation) in operations.iter().enumerate() {
                let base = operation.apply(site.frac());
                let [bx, by, bz] = base;
                let floor_x = bx.floor();
                let floor_y = by.floor();
                let floor_z = bz.floor();
                let frac = [bx - floor_x, by - floor_y, bz - floor_z];
                let tau = [-floor_x as i32, -floor_y as i32, -floor_z as i32];
                let [frac_x, frac_y, frac_z] = frac;
                let bin = [
                    bin_index(frac_x, count_x),
                    bin_index(frac_y, count_y),
                    bin_index(frac_z, count_z),
                ];
                let [bin_x, bin_y, bin_z] = bin;
                let mut duplicate = false;
                'candidates: for &offset_x in offsets_x {
                    let target_x = wrap(bin_x, offset_x, count_x);
                    for &offset_y in offsets_y {
                        let target_y = wrap(bin_y, offset_y, count_y);
                        for &offset_z in offsets_z {
                            let target_z = wrap(bin_z, offset_z, count_z);
                            let Some(bucket) =
                                grid.get(&(site.element(), [target_x, target_y, target_z]))
                            else {
                                continue;
                            };
                            for &existing_index in bucket {
                                let Some(existing) = atoms.get(existing_index) else {
                                    continue;
                                };
                                if structure.lattice().distance(existing.frac, frac) <= tolerance {
                                    duplicate = true;
                                    break 'candidates;
                                }
                            }
                        }
                    }
                }
                if duplicate {
                    continue;
                }
                let atom_index = atoms.len();
                atoms.push(P1Atom {
                    element: site.element().to_owned(),
                    label: site.label().to_owned(),
                    frac,
                    occupancy: site.occupancy(),
                    site_idx,
                    op_idx,
                    tau,
                });
                grid.entry((site.element(), bin))
                    .or_default()
                    .push(atom_index);
            }
        }
        Ok(Self {
            lattice: *structure.lattice(),
            atoms,
        })
    }

    /// Number of atoms in the P1 cell.
    #[must_use]
    pub fn n_atoms(&self) -> usize {
        self.atoms.len()
    }

    /// Element symbols in atom order.
    #[must_use]
    pub fn elements(&self) -> Vec<String> {
        self.atoms.iter().map(|atom| atom.element.clone()).collect()
    }

    /// Fractional coordinates in atom order.
    #[must_use]
    pub fn frac(&self) -> Vec<[f64; 3]> {
        self.atoms.iter().map(|atom| atom.frac).collect()
    }

    /// Occupancies in atom order.
    #[must_use]
    pub fn occupancy(&self) -> Vec<f64> {
        self.atoms.iter().map(|atom| atom.occupancy).collect()
    }

    /// Atoms generated by one asymmetric site.
    #[must_use]
    pub fn images_of(&self, site_idx: usize) -> Vec<&P1Atom> {
        self.atoms
            .iter()
            .filter(|atom| atom.site_idx == site_idx)
            .collect()
    }
}

impl Structure {
    /// Expands this structure to P1 with the default tolerance.
    ///
    /// # Errors
    ///
    /// Returns the errors of [`P1Cell::expand`].
    pub fn p1(&self) -> Result<P1Cell, Error> {
        P1Cell::expand(self, 0.05)
    }
}

fn bin_index(value: f64, count: usize) -> usize {
    let scaled = value.rem_euclid(1.0) * count as f64;
    let index = scaled.floor() as usize;
    index.min(count.saturating_sub(1))
}

fn offsets(count: usize) -> &'static [isize] {
    match count {
        0 | 1 => &[0],
        2 => &[-1, 0],
        _ => &[-1, 0, 1],
    }
}

fn wrap(index: usize, offset: isize, count: usize) -> usize {
    let count = count as isize;
    (index as isize + offset).rem_euclid(count) as usize
}
