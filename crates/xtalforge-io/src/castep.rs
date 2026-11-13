//! CASTEP `.cell` reader and writer (Task 13.2.2).
//!
//! Reads `%BLOCK LATTICE_CART`, `%BLOCK LATTICE_ABC`, `%BLOCK POSITIONS_FRAC`
//! and `%BLOCK POSITIONS_CART` (singular `POSITION_*` is accepted too).
//! Occupancies and labels after the coordinates are kept. The writer emits an
//! angstrom `LATTICE_CART` block, `POSITIONS_FRAC` and a placeholder
//! `SPECIES_POT` block.

#![allow(clippy::indexing_slicing)] // fixed-format parser, lengths checked

use serde_json::Value;
use xtalforge_core::elements::parse_symbol;
use xtalforge_core::{Error, Lattice, Site, Structure, StructureBuilder};

use crate::molecule::parse_number;

/// Reads a CASTEP `.cell` file.
///
/// # Errors
///
/// Returns [`Error::InvalidInput`] for a missing lattice or position block, a
/// non-orthogonal Cartesian cell or a malformed coordinate.
pub fn read_castep(text: &str) -> Result<Structure, Error> {
    let lines: Vec<&str> = text.lines().collect();
    let cart = find_block(&lines, "LATTICE_CART");
    let abc = find_block(&lines, "LATTICE_ABC");
    let lattice = if let Some((start, end)) = cart {
        let mut matrix = [[0.0; 3]; 3];
        let rows = block_rows(&lines, start, end);
        for (slot, line) in matrix.iter_mut().zip(rows.iter()) {
            let values: Vec<f64> = line
                .split_whitespace()
                .filter_map(|value| value.parse().ok())
                .collect();
            if values.len() < 3 {
                return Err(Error::InvalidInput(
                    "a CASTEP LATTICE_CART row needs three values".to_owned(),
                ));
            }
            *slot = [values[0], values[1], values[2]];
        }
        Lattice::new(matrix)?
    } else if let Some((start, end)) = abc {
        let rows = block_rows(&lines, start, end);
        let line = rows.first().copied().unwrap_or("");
        let values: Vec<f64> = line
            .split_whitespace()
            .filter_map(|value| value.parse().ok())
            .collect();
        if values.len() < 6 {
            return Err(Error::InvalidInput(
                "CASTEP LATTICE_ABC needs a b c alpha beta gamma".to_owned(),
            ));
        }
        Lattice::from_parameters(
            values[0], values[1], values[2], values[3], values[4], values[5],
        )?
    } else {
        return Err(Error::InvalidInput(
            "CASTEP .cell has no LATTICE_CART or LATTICE_ABC block".to_owned(),
        ));
    };
    let (positions, is_cartesian) = if let Some((start, end)) =
        find_block(&lines, "POSITIONS_FRAC").or_else(|| find_block(&lines, "POSITION_FRAC"))
    {
        (block_rows(&lines, start, end), false)
    } else if let Some((start, end)) =
        find_block(&lines, "POSITIONS_CART").or_else(|| find_block(&lines, "POSITION_CART"))
    {
        (block_rows(&lines, start, end), true)
    } else {
        return Err(Error::InvalidInput(
            "CASTEP .cell has no POSITIONS_FRAC or POSITIONS_CART block".to_owned(),
        ));
    };
    let mut builder = StructureBuilder::new("castep", lattice);
    for (index, line) in positions.iter().enumerate() {
        let fields: Vec<&str> = line.split_whitespace().collect();
        if fields.len() < 4 {
            return Err(Error::InvalidInput(
                "a CASTEP position line needs species and three coordinates".to_owned(),
            ));
        }
        let element = parse_symbol(fields[0])?;
        let cart = [
            parse_number(fields[1], "CASTEP x")?,
            parse_number(fields[2], "CASTEP y")?,
            parse_number(fields[3], "CASTEP z")?,
        ];
        let frac = if is_cartesian {
            lattice.to_frac(cart)
        } else {
            cart
        };
        let mut site = Site::new(format!("{element}{}", index + 1), &element, frac)?;
        if let Some(occupancy) = fields
            .iter()
            .skip(4)
            .find_map(|value| value.parse::<f64>().ok())
        {
            site = site.with_occupancy(occupancy.clamp(0.0, 1.0))?;
        }
        if !is_cartesian && fields.len() > 4 {
            let label = fields.get(4).copied().unwrap_or("");
            if label.parse::<f64>().is_err() {
                site = site.with_prop("castep_label", Value::String(label.to_owned()));
            }
        }
        builder = builder.site(site);
    }
    builder.build()
}

/// Writes a P1 structure as a CASTEP `.cell` file.
///
/// # Errors
///
/// Returns an error for an empty structure or a failed P1 expansion.
pub fn write_castep(structure: &Structure) -> Result<String, Error> {
    let p1 = structure.p1()?;
    if p1.atoms.is_empty() {
        return Err(Error::InvalidInput(
            "a CASTEP .cell cannot hold a structure with no atoms".to_owned(),
        ));
    }
    let mut species: Vec<String> = Vec::new();
    for atom in &p1.atoms {
        if !species.iter().any(|known| known == &atom.element) {
            species.push(atom.element.clone());
        }
    }
    let mut out = String::from("%BLOCK LATTICE_CART\n");
    for row in structure.lattice().matrix() {
        out.push_str(&format!("  {:.10} {:.10} {:.10}\n", row[0], row[1], row[2]));
    }
    out.push_str("%ENDBLOCK LATTICE_CART\n\n%BLOCK POSITIONS_FRAC\n");
    let mut used: std::collections::BTreeSet<String> = std::collections::BTreeSet::new();
    for atom in &p1.atoms {
        let label = crate::molecule::unique_label(&atom.element, &mut used);
        let occupancy = structure
            .sites()
            .get(atom.site_idx)
            .map_or(1.0, |site| site.occupancy());
        out.push_str(&format!(
            "  {} {:.10} {:.10} {:.10} {} {:.6}\n",
            atom.element, atom.frac[0], atom.frac[1], atom.frac[2], label, occupancy
        ));
    }
    out.push_str("%ENDBLOCK POSITIONS_FRAC\n\n%BLOCK SPECIES_POT\n");
    for element in &species {
        out.push_str(&format!("  {element} {element}.usp\n"));
    }
    out.push_str("%ENDBLOCK SPECIES_POT\n");
    Ok(out)
}

fn find_block(lines: &[&str], name: &str) -> Option<(usize, usize)> {
    let wanted = format!("%BLOCK {name}");
    let start = lines
        .iter()
        .position(|line| line.trim().eq_ignore_ascii_case(&wanted))?;
    let end = lines
        .iter()
        .enumerate()
        .skip(start + 1)
        .find(|(_, line)| {
            line.trim()
                .eq_ignore_ascii_case(&format!("%ENDBLOCK {name}"))
        })
        .map(|(index, _)| index)?;
    Some((start, end))
}

fn block_rows<'a>(lines: &'a [&'a str], start: usize, end: usize) -> Vec<&'a str> {
    lines
        .iter()
        .take(end)
        .skip(start + 1)
        .map(|line| line.trim())
        .filter(|line| !line.is_empty() && !line.starts_with('#'))
        .collect()
}
