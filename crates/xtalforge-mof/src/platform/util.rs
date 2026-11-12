//! Shared helpers for the MOF platform: P1 structures, bond graphs and
//! connected components.

use std::collections::BTreeMap;

use serde::{Deserialize, Serialize};
use xtalforge_core::{BondGraph, Error, P1Cell, Site, Structure, StructureBuilder, element_info};

/// Converts a P1 cell back into a P1 structure.
///
/// # Errors
///
/// Returns the errors of [`Site::new`].
pub(crate) fn p1_structure(cell: &P1Cell) -> Result<Structure, Error> {
    let mut builder = StructureBuilder::new("p1", cell.lattice);
    for atom in &cell.atoms {
        builder = builder.site(
            Site::new(atom.label.clone(), &atom.element, atom.frac)?
                .with_occupancy(atom.occupancy)?,
        );
    }
    builder.build()
}

/// The P1 expansion of a structure and its perceived bond graph.
///
/// # Errors
///
/// Returns the errors of the P1 expansion or of the bond perception.
pub(crate) fn p1_and_bonds(structure: &Structure) -> Result<(Structure, BondGraph), Error> {
    let p1 = p1_structure(&structure.p1()?)?;
    let graph = xtalforge_core::perceive_bonds(&p1, xtalforge_core::BondRules::default())?;
    Ok((p1, graph))
}

/// Element counts of a set of atom indices.
///
/// # Errors
///
/// Returns the errors of the element lookup.
pub(crate) fn element_counts(
    structure: &Structure,
    atoms: &[usize],
) -> Result<BTreeMap<String, usize>, Error> {
    let mut counts = BTreeMap::new();
    for &index in atoms {
        if let Some(site) = structure.sites().get(index) {
            let symbol = element_info(site.element())?.symbol;
            *counts.entry(symbol).or_insert(0) += 1;
        }
    }
    Ok(counts)
}

/// Formula string from element counts, Hill order.
#[must_use]
pub(crate) fn formula(counts: &BTreeMap<String, usize>) -> String {
    let mut order: Vec<&String> = counts.keys().collect();
    order.sort_by(|first, second| {
        let carbon_first = *first == "C";
        let carbon_second = *second == "C";
        let hydrogen_first = *first == "H";
        let hydrogen_second = *second == "H";
        if carbon_first != carbon_second {
            return if carbon_first {
                std::cmp::Ordering::Less
            } else {
                std::cmp::Ordering::Greater
            };
        }
        if hydrogen_first != hydrogen_second {
            return if hydrogen_first {
                std::cmp::Ordering::Less
            } else {
                std::cmp::Ordering::Greater
            };
        }
        first.cmp(second)
    });
    let mut text = String::new();
    for symbol in order {
        let count = counts.get(symbol).copied().unwrap_or(0);
        text.push_str(symbol);
        if count > 1 {
            text.push_str(&count.to_string());
        }
    }
    text
}

/// Connected components of an undirected adjacency list.
#[must_use]
pub(crate) fn components(adjacency: &[Vec<usize>]) -> Vec<Vec<usize>> {
    let mut seen = vec![false; adjacency.len()];
    let mut output = Vec::new();
    for start in 0..adjacency.len() {
        if seen.get(start).copied().unwrap_or(true) {
            continue;
        }
        let mut component = Vec::new();
        let mut stack = vec![start];
        if let Some(slot) = seen.get_mut(start) {
            *slot = true;
        }
        while let Some(node) = stack.pop() {
            component.push(node);
            for &neighbour in adjacency.get(node).map(Vec::as_slice).unwrap_or(&[]) {
                if seen.get(neighbour).copied().unwrap_or(true) {
                    continue;
                }
                if let Some(slot) = seen.get_mut(neighbour) {
                    *slot = true;
                }
                stack.push(neighbour);
            }
        }
        component.sort_unstable();
        output.push(component);
    }
    output
}

/// Adjacency list from a bond graph, adding both directions.
#[must_use]
pub(crate) fn bond_adjacency(graph: &BondGraph, atom_count: usize) -> Vec<Vec<usize>> {
    let mut adjacency = vec![Vec::new(); atom_count];
    for bond in &graph.bonds {
        if bond.i == bond.j {
            continue;
        }
        if let Some(list) = adjacency.get_mut(bond.i)
            && !list.contains(&bond.j)
        {
            list.push(bond.j);
        }
        if let Some(list) = adjacency.get_mut(bond.j)
            && !list.contains(&bond.i)
        {
            list.push(bond.i);
        }
    }
    for list in &mut adjacency {
        list.sort_unstable();
    }
    adjacency
}

/// One atom of a deconstructed unit.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct UnitAtom {
    /// P1 atom index in the parent structure.
    pub index: usize,
    /// Element symbol.
    pub element: String,
    /// Fractional coordinates.
    pub frac: [f64; 3],
}

/// One connection from a unit atom to another unit.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct UnitConnection {
    /// Atom index inside the unit.
    pub atom: usize,
    /// Parent P1 index of the atom on the other side.
    pub other: usize,
    /// Element of the atom on the other side.
    pub other_element: String,
    /// Lattice image of the other atom.
    pub image: [i32; 3],
}
