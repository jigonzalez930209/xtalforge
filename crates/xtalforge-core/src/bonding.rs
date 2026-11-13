//! Distance-based bond perception and the connectivity graph (Task 1.11).

mod explicit;
mod rules;

use std::collections::VecDeque;

pub use explicit::explicit_bond_images;
pub use rules::BondRules;

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
///
/// The indices name atoms of the P1 cell; perception expands the bond onto every
/// symmetry image of both ends, so drawing one Fe-O bond draws the whole orbit.
#[derive(Debug, Clone, Copy, PartialEq, serde::Serialize, serde::Deserialize)]
pub struct ExplicitBond {
    /// First atom index in the P1 cell.
    pub i: usize,
    /// Second atom index in the P1 cell.
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

/// Connectivity graph over the P1 cell.
#[derive(Debug, Clone, PartialEq)]
pub struct BondGraph {
    /// Perceived bonds.
    pub bonds: Vec<Bond>,
    /// Adjacency; entry `i` lists bond indices incident to atom `i`.
    adjacency: Vec<Vec<usize>>,
}

impl BondGraph {
    /// Connectivity from the explicit bonds only.
    ///
    /// This is the graph a builder states: distances are not perceived, so
    /// a close contact across a pore is not silently made a bond. The
    /// structure must be P1 (the carbon builders are), because an explicit
    /// bond names an atom index in the P1 cell.
    ///
    /// # Errors
    ///
    /// Returns [`Error::InvalidInput`] for a non-P1 structure or a failed
    /// P1 expansion.
    pub fn from_explicit(structure: &Structure) -> Result<Self, Error> {
        if !structure.is_p1() {
            return Err(Error::InvalidInput(
                "an explicit-only bond graph needs a P1 structure".to_owned(),
            ));
        }
        let p1 = structure.p1()?;
        let mut bonds = Vec::new();
        for explicit in structure.explicit_bonds() {
            if explicit.remove {
                continue;
            }
            let distance = explicit_distance(&p1, explicit).unwrap_or(0.0);
            bonds.push(Bond {
                i: explicit.i,
                j: explicit.j,
                image: explicit.image,
                distance,
                order: explicit.order,
                explicit: true,
            });
        }
        let mut adjacency = vec![Vec::new(); p1.atoms.len()];
        for (index, bond) in bonds.iter().enumerate() {
            if let Some(list) = adjacency.get_mut(bond.i) {
                list.push(index);
            }
            if bond.j != bond.i
                && let Some(list) = adjacency.get_mut(bond.j)
            {
                list.push(index);
            }
        }
        Ok(Self { bonds, adjacency })
    }

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
/// Distance perception uses [`BondRules::range`] for each element pair, so a
/// pair-specific range and the metal-metal switch are honoured. Every stored
/// explicit bond is then expanded over the symmetry group with
/// [`explicit_bond_images`], which is how a bond drawn on one image of an atom
/// reaches its whole orbit; a suppression removes every image it draws.
///
/// # Errors
///
/// Returns [`Error::InvalidInput`] for invalid rules or a failed P1 expansion.
pub fn perceive_bonds(structure: &Structure, rules: BondRules) -> Result<BondGraph, Error> {
    rules.validate()?;
    let p1 = structure.p1()?;
    let symbols: Vec<String> = p1.atoms.iter().map(|atom| atom.element.clone()).collect();
    let elements: Vec<&'static crate::elements::Element> = symbols
        .iter()
        .map(|symbol| element_by_symbol(symbol))
        .collect::<Result<_, _>>()?;
    let dummy: Vec<bool> = elements.iter().map(|element| element.is_dummy()).collect();
    let cutoff = rules.max_cutoff(&symbols)?;
    let mut bonds = Vec::new();
    if cutoff > 0.0 && p1.atoms.len() > 1 {
        let found = neighbors(&p1_to_structure(&p1)?, cutoff)?;
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
            let (Some(first), Some(second)) = (elements.get(i), elements.get(j)) else {
                continue;
            };
            if !rules.allows_elements(first, second) {
                continue;
            }
            let (minimum, maximum) = rules.range_of(first, second);
            if minimum <= distance && distance <= maximum {
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
    }
    for stored in structure.explicit_bonds() {
        for explicit in explicit_bond_images(structure, &p1, stored) {
            apply_explicit(&mut bonds, &p1, &explicit);
        }
    }
    let mut adjacency = vec![Vec::new(); p1.atoms.len()];
    for (index, bond) in bonds.iter().enumerate() {
        if let Some(list) = adjacency.get_mut(bond.i) {
            list.push(index);
        }
        if bond.j != bond.i
            && let Some(list) = adjacency.get_mut(bond.j)
        {
            list.push(index);
        }
    }
    Ok(BondGraph { bonds, adjacency })
}

impl Structure {
    /// Perceives bonds with the rules stored on the structure.
    ///
    /// The rules come from the `bond_rules` metadata [`BondRules::from_structure`]
    /// reads, so a pair range or the metal-metal switch set by the editor applies
    /// here too; a structure with no stored rules uses the defaults.
    ///
    /// # Errors
    ///
    /// Returns the errors of [`perceive_bonds`].
    pub fn bond_graph(&self) -> Result<BondGraph, Error> {
        perceive_bonds(self, BondRules::from_structure(self)?)
    }
}

fn apply_explicit(bonds: &mut Vec<Bond>, p1: &crate::P1Cell, explicit: &ExplicitBond) {
    let distance = explicit_distance(p1, explicit);
    if explicit.remove {
        let exact = bonds
            .iter()
            .any(|bond| same_bond(bond, explicit.i, explicit.j, explicit.image));
        let target = contact_vector(p1, explicit.i, explicit.j, explicit.image);
        bonds.retain(|bond| {
            if exact {
                !same_bond(bond, explicit.i, explicit.j, explicit.image)
            } else {
                !same_contact(bond, p1, target)
            }
        });
        return;
    }
    if let Some(existing) = bonds
        .iter_mut()
        .find(|bond| same_bond(bond, explicit.i, explicit.j, explicit.image))
    {
        existing.order = explicit.order;
        existing.explicit = true;
        return;
    }
    let target = contact_vector(p1, explicit.i, explicit.j, explicit.image);
    if let Some(existing) = bonds.iter_mut().find(|bond| same_contact(bond, p1, target)) {
        existing.order = explicit.order;
        existing.explicit = true;
        return;
    }
    let Some(distance) = distance else {
        return;
    };
    let (i, j, image) = canonical(explicit.i, explicit.j, explicit.image);
    if bonds.iter().any(|bond| same_bond(bond, i, j, image)) {
        return;
    }
    bonds.push(Bond {
        i,
        j,
        image,
        distance,
        order: explicit.order,
        explicit: true,
    });
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

pub(crate) fn canonical(i: usize, j: usize, image: [i32; 3]) -> (usize, usize, [i32; 3]) {
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
fn same_contact(bond: &Bond, p1: &crate::P1Cell, target: Option<[f64; 3]>) -> bool {
    let Some(target) = target else {
        return false;
    };
    let Some(vector) = contact_vector(p1, bond.i, bond.j, bond.image) else {
        return false;
    };
    vector
        .iter()
        .zip(target)
        .all(|(a, b)| (a - b).abs() <= 1e-4)
}

/// Canonical fractional vector of a periodic contact.
///
/// Two lattice images of the same physical contact produce the same vector
/// after wrapping each component to `[-0.5, 0.5)` and orienting the pair.
fn contact_vector(p1: &crate::P1Cell, i: usize, j: usize, image: [i32; 3]) -> Option<[f64; 3]> {
    let first = p1.atoms.get(i)?;
    let second = p1.atoms.get(j)?;
    let mut delta = [
        second.frac[0] + f64::from(image[0]) - first.frac[0],
        second.frac[1] + f64::from(image[1]) - first.frac[1],
        second.frac[2] + f64::from(image[2]) - first.frac[2],
    ];
    for value in &mut delta {
        *value -= value.round();
    }
    if i > j
        || (i == j
            && (delta[0] < 0.0
                || (delta[0] == 0.0 && delta[1] < 0.0)
                || (delta[0] == 0.0 && delta[1] == 0.0 && delta[2] < 0.0)))
    {
        delta = [-delta[0], -delta[1], -delta[2]];
    }
    Some(delta)
}
