//! DFTB+ `.gen` reader and writer (Task 3.3.3).

#![allow(clippy::indexing_slicing)] // fixed-format parser, lengths checked above
use std::collections::BTreeMap;

use xtalforge_core::{Error, Lattice, Site, Structure, StructureBuilder};

/// Reads a DFTB+ `.gen` geometry (`C`, `S` or `F`).
///
/// # Errors
///
/// Returns [`Error::InvalidInput`] for a bad header, species index or line
/// count.
pub fn read_gen(text: &str) -> Result<Structure, Error> {
    let rows: Vec<&str> = text
        .lines()
        .map(|line| line.split('#').next().unwrap_or("").trim())
        .filter(|line| !line.is_empty())
        .collect();
    if rows.len() < 3 {
        return Err(Error::InvalidInput(
            "a .gen file needs a count line, a species line and atoms".to_owned(),
        ));
    }
    let head: Vec<&str> = rows[0].split_whitespace().collect();
    let count: usize = head
        .first()
        .and_then(|value| value.parse().ok())
        .ok_or_else(|| Error::InvalidInput("bad .gen atom count".to_owned()))?;
    let kind = head.get(1).copied().unwrap_or("C").to_ascii_uppercase();
    if !matches!(kind.as_str(), "C" | "S" | "F") {
        return Err(Error::InvalidInput(format!(
            "{kind:?} is not a .gen geometry kind; use C, S or F"
        )));
    }
    let species: Vec<String> = rows[1].split_whitespace().map(str::to_owned).collect();
    if species.is_empty() {
        return Err(Error::InvalidInput(
            "the .gen species line is empty".to_owned(),
        ));
    }
    let mut points = Vec::with_capacity(count);
    for index in 0..count {
        let line = rows
            .get(2 + index)
            .ok_or_else(|| Error::InvalidInput(".gen ended before all atoms".to_owned()))?;
        let parts: Vec<&str> = line.split_whitespace().collect();
        if parts.len() < 5 {
            return Err(Error::InvalidInput(format!(
                "atom {} needs an index, species number and three coordinates",
                index + 1
            )));
        }
        let species_index: usize = parts[1]
            .parse::<usize>()
            .map_err(|_| Error::InvalidInput("bad .gen species number".to_owned()))?
            .checked_sub(1)
            .ok_or_else(|| Error::InvalidInput(".gen species numbers are one-based".to_owned()))?;
        let element = species
            .get(species_index)
            .ok_or_else(|| Error::InvalidInput("bad .gen species index".to_owned()))?;
        let coordinates = [
            parts[2]
                .parse::<f64>()
                .map_err(|_| Error::InvalidInput("bad .gen x".to_owned()))?,
            parts[3]
                .parse::<f64>()
                .map_err(|_| Error::InvalidInput("bad .gen y".to_owned()))?,
            parts[4]
                .parse::<f64>()
                .map_err(|_| Error::InvalidInput("bad .gen z".to_owned()))?,
        ];
        points.push((element.clone(), coordinates));
    }
    if kind == "C" {
        return Structure::molecule("gen", &points, 5.0);
    }
    let origin = rows
        .get(2 + count)
        .ok_or_else(|| Error::InvalidInput(".gen needs an origin and a lattice".to_owned()))?;
    let mut matrix = [[0.0; 3]; 3];
    for (row_index, row) in matrix.iter_mut().enumerate() {
        let Some(line) = rows.get(3 + count + row_index) else {
            return Err(Error::InvalidInput(".gen lattice is incomplete".to_owned()));
        };
        let values: Vec<f64> = line
            .split_whitespace()
            .filter_map(|v| v.parse().ok())
            .collect();
        if values.len() < 3 {
            return Err(Error::InvalidInput(
                "a .gen lattice row needs three values".to_owned(),
            ));
        }
        *row = [values[0], values[1], values[2]];
    }
    let _ = origin;
    let lattice = Lattice::new(matrix)?;
    let mut builder = StructureBuilder::new("gen", lattice);
    for (index, (element, coordinates)) in points.iter().enumerate() {
        let frac = if kind == "F" {
            *coordinates
        } else {
            lattice.to_frac(*coordinates)
        };
        builder = builder.site(Site::new(format!("{element}{}", index + 1), element, frac)?);
    }
    builder.build()
}

/// Writes a P1 `.gen` with fractional coordinates.
///
/// # Errors
///
/// Returns the errors of the P1 expansion or an empty structure.
pub fn write_gen(structure: &Structure) -> Result<String, Error> {
    let p1 = structure.p1()?;
    if p1.atoms.is_empty() {
        return Err(Error::InvalidInput(
            "a .gen file cannot hold a structure with no atoms".to_owned(),
        ));
    }
    let mut index: BTreeMap<&str, usize> = BTreeMap::new();
    for atom in &p1.atoms {
        let next = index.len() + 1;
        index.entry(atom.element.as_str()).or_insert(next);
    }
    let mut out = format!("{} F\n", p1.atoms.len());
    let mut order: Vec<(&str, usize)> = index.iter().map(|(name, index)| (*name, *index)).collect();
    order.sort_by_key(|(_, index)| *index);
    out.push_str(
        &order
            .iter()
            .map(|(name, _)| *name)
            .collect::<Vec<_>>()
            .join(" "),
    );
    out.push('\n');
    for (atom_index, atom) in p1.atoms.iter().enumerate() {
        let kind = index.get(atom.element.as_str()).copied().unwrap_or(1);
        let [x, y, z] = atom.frac;
        out.push_str(&format!(
            "{:6} {:3} {x:20.12} {y:20.12} {z:20.12}\n",
            atom_index + 1,
            kind
        ));
    }
    out.push_str(&format!("{:20.12} {:20.12} {:20.12}\n", 0.0, 0.0, 0.0));
    for row in structure.lattice().matrix() {
        out.push_str(&format!(
            "{:20.12} {:20.12} {:20.12}\n",
            row[0], row[1], row[2]
        ));
    }
    Ok(out)
}
