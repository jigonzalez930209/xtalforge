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

/// An explicit bond stored by the user.
#[derive(Debug, Clone, Copy, PartialEq, serde::Serialize, serde::Deserialize)]
pub struct ExplicitBond {
    /// First atom index in the asymmetric unit or P1 cell.
    pub i: usize,
    /// Second atom index.
    pub j: usize,
    /// Lattice image of `j`.
    pub image: [i32; 3],
    /// Bond order.
    pub order: f64,
    /// When true, remove a perceived bond instead of adding one.
    pub remove: bool,
}

/// A topology edge between node centres; never a chemical bond.
#[derive(Debug, Clone, Copy, PartialEq, serde::Serialize, serde::Deserialize)]
pub struct TopologyBond {
    /// First node.
    pub i: usize,
    /// Second node.
    pub j: usize,
    /// Lattice image of `j`.
    pub image: [i32; 3],
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
            scale: 1.15,
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
    let mut metal = Vec::with_capacity(p1.atoms.len());
    for atom in &p1.atoms {
        let element = element_by_symbol(&atom.element)?;
        radii.push(element.covalent_radius);
        dummy.push(atom.element == "X");
        metal.push(element.is_metal);
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
            || (metal.get(i).copied().unwrap_or(false) && metal.get(j).copied().unwrap_or(false))
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
    for explicit in structure.explicit_bonds() {
        let distance = explicit_distance(&p1, explicit);
        if explicit.remove {
            bonds.retain(|bond| {
                !same_bond(bond, explicit.i, explicit.j, explicit.image)
                    && !same_contact(bond, explicit, distance)
            });
            continue;
        }
        if let Some(existing) = bonds
            .iter_mut()
            .find(|bond| same_bond(bond, explicit.i, explicit.j, explicit.image))
        {
            existing.order = explicit.order;
            existing.explicit = true;
            continue;
        }
        if let Some(existing) = bonds
            .iter_mut()
            .find(|bond| same_contact(bond, explicit, distance))
        {
            existing.order = explicit.order;
            existing.explicit = true;
            continue;
        }
        if !structure.is_p1() {
            // A non-P1 record names a site and an operation; the P1 mapping
            // may land on a different image representative. Only an existing
            // perceived contact is safe to annotate; an unmatched record is
            // kept in memory but does not invent a bond.
            continue;
        }
        if let (Some(distance), Some(_)) = (distance, p1.atoms.get(explicit.j)) {
            let (i, j, image) = canonical(explicit.i, explicit.j, explicit.image);
            bonds.push(Bond {
                i,
                j,
                image,
                distance,
                order: explicit.order,
                explicit: true,
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

fn same_bond(bond: &Bond, i: usize, j: usize, image: [i32; 3]) -> bool {
    (bond.i == i && bond.j == j && bond.image == image)
        || (bond.i == j && bond.j == i && bond.image == [-image[0], -image[1], -image[2]])
}

fn canonical(i: usize, j: usize, image: [i32; 3]) -> (usize, usize, [i32; 3]) {
    if i > j
        || (i == j
            && (image[0] < 0
                || (image[0] == 0 && image[1] < 0)
                || (image[0] == 0 && image[1] == 0 && image[2] < 0)))
    {
        (j, i, [-image[0], -image[1], -image[2]])
    } else {
        (i, j, image)
    }
}

fn explicit_distance(cell: &crate::P1Cell, explicit: &ExplicitBond) -> Option<f64> {
    let first = cell.atoms.get(explicit.i)?;
    let second = cell.atoms.get(explicit.j)?;
    let delta = [
        second.frac[0] + f64::from(explicit.image[0]) - first.frac[0],
        second.frac[1] + f64::from(explicit.image[1]) - first.frac[1],
        second.frac[2] + f64::from(explicit.image[2]) - first.frac[2],
    ];
    let matrix = cell.lattice.matrix();
    let cart = [
        delta[0] * matrix[0][0] + delta[1] * matrix[1][0] + delta[2] * matrix[2][0],
        delta[0] * matrix[0][1] + delta[1] * matrix[1][1] + delta[2] * matrix[2][1],
        delta[0] * matrix[0][2] + delta[1] * matrix[1][2] + delta[2] * matrix[2][2],
    ];
    Some(cart.iter().map(|value| value * value).sum::<f64>().sqrt())
}

/// The same atom pair with a length that agrees to a milliangstrom.
///
/// A non-P1 CIF stores the operation and image of the second site, and the
/// P1 expansion may choose the neighbouring cell's representative for one
/// end. The pair and length identify the contact; the perceived image is kept
/// because it is the one the graph is drawn with.
fn same_contact(bond: &Bond, explicit: &ExplicitBond, distance: Option<f64>) -> bool {
    let same_pair = (bond.i == explicit.i && bond.j == explicit.j)
        || (bond.i == explicit.j && bond.j == explicit.i);
    same_pair && distance.is_some_and(|value| (bond.distance - value).abs() <= 1e-3)
}
