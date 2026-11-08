//! pymatgen JSON structures without importing pymatgen (Task 3.3.1).

#![allow(clippy::indexing_slicing)] // fixed-format parser, lengths checked above
use serde_json::{Value, json};
use xtalforge_core::{Error, Lattice, Site, Structure, StructureBuilder};

/// Reads a pymatgen `Structure` JSON document.
///
/// # Errors
///
/// Returns [`Error::InvalidInput`] for missing lattice/sites or bad values.
pub fn read_pmg_json(text: &str) -> Result<Structure, Error> {
    let value: Value = serde_json::from_str(text)
        .map_err(|error| Error::InvalidInput(format!("invalid pymatgen JSON: {error}")))?;
    let matrix = value["lattice"]["matrix"]
        .as_array()
        .ok_or_else(|| Error::InvalidInput("pymatgen JSON has no lattice matrix".to_owned()))?;
    let mut rows = [[0.0; 3]; 3];
    for (slot, row) in rows.iter_mut().zip(matrix.iter()) {
        let values = row
            .as_array()
            .ok_or_else(|| Error::InvalidInput("lattice row is not an array".to_owned()))?;
        for (axis, value) in values.iter().take(3).enumerate() {
            *slot.get_mut(axis).ok_or_else(|| {
                Error::InvalidInput("lattice row has more than three entries".to_owned())
            })? = value
                .as_f64()
                .ok_or_else(|| Error::InvalidInput("lattice entry is not a number".to_owned()))?;
        }
    }
    let lattice = Lattice::new(rows)?;
    let sites = value["sites"]
        .as_array()
        .ok_or_else(|| Error::InvalidInput("pymatgen JSON has no sites".to_owned()))?;
    let mut builder = StructureBuilder::new("pmg", lattice);
    for (index, site) in sites.iter().enumerate() {
        let abc = site["abc"]
            .as_array()
            .ok_or_else(|| Error::InvalidInput("site has no abc coordinates".to_owned()))?;
        let mut frac = [0.0; 3];
        for (axis, value) in abc.iter().take(3).enumerate() {
            *frac.get_mut(axis).ok_or_else(|| {
                Error::InvalidInput("site abc has more than three entries".to_owned())
            })? = value
                .as_f64()
                .ok_or_else(|| Error::InvalidInput("site coordinate is not a number".to_owned()))?;
        }
        let (element, occupancy) = species(site)?;
        let label = format!("{element}{}", index + 1);
        builder = builder.site(Site::new(label, &element, frac)?.with_occupancy(occupancy)?);
    }
    builder.build()
}

fn species(site: &Value) -> Result<(String, f64), Error> {
    let entries = site["species"]
        .as_array()
        .ok_or_else(|| Error::InvalidInput("site has no species".to_owned()))?;
    let first = entries
        .first()
        .ok_or_else(|| Error::InvalidInput("site species is empty".to_owned()))?;
    let element = first["element"]
        .as_str()
        .or_else(|| first["name"].as_str())
        .ok_or_else(|| Error::InvalidInput("species has no element".to_owned()))?;
    let occupancy = first["occu"].as_f64().unwrap_or(1.0).clamp(0.0, 1.0);
    Ok((element.to_owned(), occupancy))
}

/// Writes a P1 pymatgen `Structure` JSON document.
///
/// # Errors
///
/// Returns the errors of the P1 expansion.
pub fn write_pmg_json(structure: &Structure) -> Result<String, Error> {
    let p1 = structure.p1()?;
    let matrix: Vec<Vec<f64>> = p1.lattice.matrix().iter().map(|row| row.to_vec()).collect();
    let sites: Vec<Value> = p1
        .atoms
        .iter()
        .map(|atom| {
            json!({
                "species": [{"element": atom.element, "occu": atom.occupancy}],
                "abc": atom.frac,
                "xyz": p1.lattice.to_cart(atom.frac),
                "label": atom.label,
            })
        })
        .collect();
    let document = json!({
        "@module": "pymatgen.core.structure",
        "@class": "Structure",
        "lattice": {"matrix": matrix, "pbc": [true, true, true]},
        "sites": sites,
    });
    serde_json::to_string_pretty(&document)
        .map_err(|error| Error::InvalidInput(format!("could not serialize pymatgen JSON: {error}")))
}
