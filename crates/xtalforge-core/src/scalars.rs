//! Named per-atom and per-bond scalars (Task 1.14.2).

use crate::elements::element_by_symbol;
use crate::{Error, Structure};

/// Bond lengths in the perceived graph.
///
/// # Errors
///
/// Returns the errors of the P1 expansion or bond perception.
pub fn bond_lengths(structure: &Structure) -> Result<Vec<f64>, Error> {
    Ok(structure
        .bond_graph()?
        .bonds
        .iter()
        .map(|bond| bond.distance)
        .collect())
}

/// Coordination number of every P1 atom.
///
/// # Errors
///
/// Returns the errors of the P1 expansion or bond perception.
pub fn coordination_numbers(structure: &Structure) -> Result<Vec<usize>, Error> {
    let graph = structure.bond_graph()?;
    Ok((0..graph.n_atoms())
        .map(|index| graph.degree(index))
        .collect())
}

/// Formal charge of every asymmetric site, when present.
#[must_use]
pub fn atom_charges(structure: &Structure) -> Vec<Option<f64>> {
    structure.sites().iter().map(|site| site.charge()).collect()
}

/// Molar mass of every asymmetric site.
///
/// # Errors
///
/// Returns an error if an element is missing from the table.
pub fn atom_masses(structure: &Structure) -> Result<Vec<f64>, Error> {
    structure
        .sites()
        .iter()
        .map(|site| Ok(element_by_symbol(site.element())?.mass))
        .collect()
}
