//! Carrying parent bonds onto the slab copies (Task 14.3.1).
//!
//! A bond joins atom `i` to the point `frac[j] + image`. Both ends have a copy
//! in every layer, and the bond survives only when the partner's copy is in the
//! slab too. A bond whose partner is outside is severed by the cut and counted;
//! the atoms are left as they are because capping them is chemistry.

use crate::bonding::{BondGraph, ExplicitBond, explicit_bond_images};
use crate::{P1Cell, Structure};

/// Everything one slab cut needs to map parent bonds onto its copies.
pub(super) struct SlabCarry<'a> {
    /// Parent structure, for its stored suppressions.
    pub(super) structure: &'a Structure,
    /// P1 cell the graph was perceived on.
    pub(super) p1: &'a P1Cell,
    /// Perceived bonds of the parent cell.
    pub(super) graph: &'a BondGraph,
    /// Number of layers.
    pub(super) layers: u32,
    /// Whole surface cells each atom sits away from its first-layer copy.
    pub(super) whole: &'a [[i32; 3]],
    /// In-plane whole cell of every `(layer, atom)` copy.
    pub(super) cells: &'a [[i32; 2]],
    /// Inverse of the surface basis.
    pub(super) inverse: [[f64; 3]; 3],
    /// Atoms per layer.
    pub(super) n_atoms: usize,
}

impl SlabCarry<'_> {
    /// Carried bonds and the number of chemical bonds the cut severed.
    pub(super) fn bonds(&self) -> (Vec<ExplicitBond>, usize) {
        let dummy: Vec<bool> = self
            .p1
            .atoms
            .iter()
            .map(|atom| atom.element == "X")
            .collect();
        let carrying = Carrying {
            whole: self.whole,
            cells: self.cells,
            inverse: self.inverse,
            n_atoms: self.n_atoms,
            layers: self.layers,
        };
        let mut carried: Vec<ExplicitBond> = Vec::new();
        let mut cut = 0;
        for bond in &self.graph.bonds {
            let chemical = !dummy.get(bond.i).copied().unwrap_or(false)
                && !dummy.get(bond.j).copied().unwrap_or(false);
            let reverse = [-bond.image[0], -bond.image[1], -bond.image[2]];
            for layer in 0..self.layers {
                let forward = carrying.carried(bond.i, bond.j, bond.image, layer);
                if let Some((i, j, image)) = forward {
                    carried.push(ExplicitBond {
                        i,
                        j,
                        image,
                        order: bond.order,
                        remove: false,
                    });
                } else if chemical {
                    cut += 1;
                }
                if chemical && carrying.carried(bond.j, bond.i, reverse, layer).is_none() {
                    cut += 1;
                }
            }
        }
        for stored in self.structure.explicit_bonds() {
            if !stored.remove {
                continue;
            }
            for removal in explicit_bond_images(self.structure, self.p1, stored) {
                for layer in 0..self.layers {
                    if let Some((i, j, image)) =
                        carrying.carried(removal.i, removal.j, removal.image, layer)
                    {
                        carried.push(ExplicitBond {
                            i,
                            j,
                            image,
                            order: removal.order,
                            remove: true,
                        });
                    }
                }
            }
        }
        (carried, cut)
    }
}

/// Maps a bond of the parent cell onto the slab copies.
struct Carrying<'a> {
    whole: &'a [[i32; 3]],
    cells: &'a [[i32; 2]],
    inverse: [[f64; 3]; 3],
    n_atoms: usize,
    layers: u32,
}

impl Carrying<'_> {
    /// Maps the bond `i -> j + image` onto the copy whose first atom is in
    /// `layer`.
    ///
    /// Returns `None` when the partner lies outside the slab: that bond is one
    /// the cut severs. The third component of the returned image is always zero
    /// because the layer is named by the atom index, not by a translation.
    fn carried(
        &self,
        i: usize,
        j: usize,
        image: [i32; 3],
        layer: u32,
    ) -> Option<(usize, usize, [i32; 3])> {
        let start = self.whole.get(i)?;
        let end = self.whole.get(j)?;
        let shifted = super::multiply_row(image.map(f64::from), self.inverse);
        let step = [
            (f64::from(end[0] - start[0]) + shifted[0]).round() as i32,
            (f64::from(end[1] - start[1]) + shifted[1]).round() as i32,
            (f64::from(end[2] - start[2]) + shifted[2]).round() as i32,
        ];
        let target = i64::from(layer) + i64::from(step[2]);
        if target < 0 || target >= i64::from(self.layers) {
            return None;
        }
        let target = target as usize;
        let from = layer as usize * self.n_atoms + i;
        let to = target * self.n_atoms + j;
        let cell_from = self.cells.get(from)?;
        let cell_to = self.cells.get(to)?;
        Some((
            from,
            to,
            [
                step[0] + cell_to[0] - cell_from[0],
                step[1] + cell_to[1] - cell_from[1],
                0,
            ],
        ))
    }
}
