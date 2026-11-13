//! OPTIMADE `structures` entry writer (Task 13.3.3).

use serde_json::{Value, json};

use xtalforge_core::{Error, Structure};

/// Serializes a structure as one OPTIMADE `structures` entry.
///
/// # Errors
///
/// Returns [`Error::InvalidInput`] when the P1 expansion fails or the JSON
/// encoder rejects a value.
pub fn write_optimade(structure: &Structure) -> Result<String, Error> {
    let p1 = structure.p1()?;
    let lattice = p1.lattice.matrix();
    let mut positions = Vec::with_capacity(p1.atoms.len());
    let mut species = Vec::with_capacity(p1.atoms.len());
    let mut elements = Vec::<String>::new();
    for atom in &p1.atoms {
        let cart = p1.lattice.to_cart(atom.frac);
        positions.push([cart[0], cart[1], cart[2]]);
        species.push(atom.element.clone());
        if !elements.contains(&atom.element) {
            elements.push(atom.element.clone());
        }
    }
    let document: Value = json!({
        "id": format!("xtalforge-{}", structure.name()),
        "type": "structures",
        "attributes": {
            "chemical_formula_reduced": structure.formula(),
            "lattice_vectors": lattice,
            "cartesian_site_positions": positions,
            "species_at_sites": species,
            "nsites": p1.atoms.len(),
            "elements": elements,
        },
    });
    serde_json::to_string(&document)
        .map_err(|error| Error::InvalidInput(format!("could not encode OPTIMADE: {error}")))
}
