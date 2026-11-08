//! Distance-based bond perception and the connectivity graph (Task 1.11).

use std::collections::VecDeque;

use crate::elements::element_by_symbol;
use crate::{Error, Neighbor, Structure, neighbors};

/// One perceived or explicit bond.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Bond {
    /// First atom.
    pub i: usize,
    /// Second atom.
    pub j: usize,
    /// Lattice image of `j`.
    pub image: [i32; 3],
    /// Distance in Angstrom.
    pub distance: f64,
    /// Bond order.
    pub order: f64,
    /// Whether the bond was stored explicitly rather than perceived.
    pub explicit: bool,
}

/// Perception rules.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct BondRules {
    /// Multiplier on the sum of covalent radii.
    pub scale: f64,
    /// Bonds shorter than this are ignored.
    pub minimum_distance: f64,
}

impl Default for BondRules {
    fn default() -> Self {
        Self {
            scale: 1.2,
            minimum_distance: 0.4,
        }
    }
}

/// Connectivity graph over the P1 cell.
#[derive(Debug, Clone, PartialEq)]
pub struct BondGraph {
    /// Perceived bonds.
    pub bonds: Vec<Bond>,
    /// Adjacency; entry `i` lists bond indices incident to atom `i`.
    adjacency: Vec<Vec<usize>>,
}

impl BondGraph {
    /// Number of atoms.
    #[must_use]
    pub fn n_atoms(&self) -> usize {
        self.adjacency.len()
    }

    /// Coordination number of one atom.
    #[must_use]
    pub fn degree(&self, index: usize) -> usize {
        self.adjacency.get(index).map_or(0, Vec::len)
    }

    /// Bond indices incident to one atom.
    #[must_use]
    pub fn bonds_of(&self, index: usize) -> &[usize] {
        self.adjacency.get(index).map_or(&[], Vec::as_slice)
    }

    /// Connected components as sorted atom-index lists.
    #[must_use]
    pub fn components(&self) -> Vec<Vec<usize>> {
        let mut seen = vec![false; self.n_atoms()];
        let mut components = Vec::new();
        for start in 0..self.n_atoms() {
            if seen.get(start).copied().unwrap_or(true) {
                continue;
            }
            let mut component = Vec::new();
            let mut queue = VecDeque::from([start]);
            while let Some(index) = queue.pop_front() {
                if seen.get(index).copied().unwrap_or(true) {
                    continue;
                }
                if let Some(slot) = seen.get_mut(index) {
                    *slot = true;
                }
                component.push(index);
                for bond_index in self.bonds_of(index) {
                    if let Some(bond) = self.bonds.get(*bond_index) {
                        let other = if bond.i == index { bond.j } else { bond.i };
                        if !seen.get(other).copied().unwrap_or(true) {
                            queue.push_back(other);
                        }
                    }
                }
            }
            component.sort_unstable();
            components.push(component);
        }
        components
    }
}

/// Perceives distance-based bonds on the P1 cell.
///
/// # Errors
///
/// Returns [`Error::InvalidInput`] for invalid rules or a failed P1 expansion.
pub fn perceive_bonds(structure: &Structure, rules: BondRules) -> Result<BondGraph, Error> {
    if !rules.scale.is_finite() || rules.scale <= 0.0 {
        return Err(Error::InvalidInput(
            "bond scale must be positive and finite".to_owned(),
        ));
    }
    if !rules.minimum_distance.is_finite() || rules.minimum_distance < 0.0 {
        return Err(Error::InvalidInput(
            "bond minimum distance must be non-negative and finite".to_owned(),
        ));
    }
    let p1 = structure.p1()?;
    let mut radii = Vec::with_capacity(p1.atoms.len());
    let mut dummy = Vec::with_capacity(p1.atoms.len());
    for atom in &p1.atoms {
        radii.push(element_by_symbol(&atom.element)?.covalent_radius);
        dummy.push(atom.element == "X");
    }
    let maximum = radii.iter().copied().fold(0.0, f64::max);
    let cutoff = 2.0 * maximum * rules.scale;
    if cutoff <= 0.0 {
        return Ok(BondGraph {
            bonds: Vec::new(),
            adjacency: vec![Vec::new(); p1.atoms.len()],
        });
    }
    let found = neighbors(&p1_to_structure(&p1)?, cutoff)?;
    let mut bonds = Vec::new();
    for Neighbor {
        i,
        j,
        image,
        distance,
    } in found
    {
        if i == j
            || distance < rules.minimum_distance
            || dummy.get(i).copied().unwrap_or(true)
            || dummy.get(j).copied().unwrap_or(true)
        {
            continue;
        }
        let Some(radius_i) = radii.get(i) else {
            continue;
        };
        let Some(radius_j) = radii.get(j) else {
            continue;
        };
        if distance <= (radius_i + radius_j) * rules.scale {
            bonds.push(Bond {
                i,
                j,
                image,
                distance,
                order: 1.0,
                explicit: false,
            });
        }
    }
    let mut adjacency = vec![Vec::new(); p1.atoms.len()];
    for (index, bond) in bonds.iter().enumerate() {
        if let Some(list) = adjacency.get_mut(bond.i) {
            list.push(index);
        }
        if let Some(list) = adjacency.get_mut(bond.j) {
            list.push(index);
        }
    }
    Ok(BondGraph { bonds, adjacency })
}

impl Structure {
    /// Perceives bonds with the default rules.
    ///
    /// # Errors
    ///
    /// Returns the errors of [`perceive_bonds`].
    pub fn bond_graph(&self) -> Result<BondGraph, Error> {
        perceive_bonds(self, BondRules::default())
    }
}

fn p1_to_structure(cell: &crate::P1Cell) -> Result<Structure, Error> {
    let mut builder = crate::StructureBuilder::new("p1", cell.lattice);
    for (index, atom) in cell.atoms.iter().enumerate() {
        builder = builder.site(
            crate::Site::new(
                format!("{}{}", atom.element, index + 1),
                &atom.element,
                atom.frac,
            )?
            .with_occupancy(atom.occupancy)?,
        );
    }
    builder.build()
}
