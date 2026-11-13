//! Quantum ESPRESSO `pw.x` output reader (Task 13.2.1).
//!
//! Reads the final cell and positions, plus the last total energy. The cell
//! comes from the last `CELL_PARAMETERS` card of a variable-cell run, or from
//! `crystal axes` scaled by `lattice parameter (alat)`. Positions come from
//! the last `ATOMIC_POSITIONS` block, or from the `site n.` table of older
//! outputs. Forces and stress stay in the OUTCAR/vasprun readers.

#![allow(clippy::indexing_slicing)] // fixed-format parser, lengths checked

use serde_json::Value;
use xtalforge_core::elements::parse_symbol;
use xtalforge_core::{Error, Lattice, Site, Structure, StructureBuilder};

use super::{BOHR, RYDBERG_EV, card_unit, parenthesised, three_numbers, unit_scale};

pub(super) fn read(text: &str) -> Result<Structure, Error> {
    let lines: Vec<&str> = text.lines().collect();
    let alat = alat_angstrom(&lines);
    let lattice = output_lattice(&lines, alat)?;
    let nat = atoms_per_cell(&lines);
    let positions = output_positions(&lines, nat, &lattice, alat)?;
    if positions.is_empty() {
        return Err(Error::InvalidInput(
            "pw.x output has no atomic positions".to_owned(),
        ));
    }
    if let Some(count) = nat
        && count != positions.len()
    {
        return Err(Error::InvalidInput(format!(
            "pw.x output declares {count} atoms but prints {} positions",
            positions.len()
        )));
    }
    let mut builder = StructureBuilder::new("pwscf", lattice);
    for (index, (symbol, frac)) in positions.iter().enumerate() {
        let element = parse_symbol(symbol)?;
        builder = builder.site(Site::new(
            format!("{element}{}", index + 1),
            &element,
            *frac,
        )?);
    }
    builder = builder.metadata("pwscf_output", Value::Bool(true));
    if let Some(ry) = total_energy(&lines) {
        builder = builder.metadata("pwscf_energy_ry", Value::from(ry));
        builder = builder.metadata("pwscf_energy_ev", Value::from(ry * RYDBERG_EV));
    }
    builder.build()
}

fn alat_angstrom(lines: &[&str]) -> Option<f64> {
    for line in lines {
        let lower = line.to_ascii_lowercase();
        if (lower.contains("lattice parameter (alat)") || lower.contains("celldm(1)"))
            && let Some(value) = value_after_equals(line)
        {
            return Some(value * BOHR);
        }
    }
    None
}

fn value_after_equals(line: &str) -> Option<f64> {
    let (_, rest) = line.split_once('=')?;
    rest.split_whitespace()
        .find_map(|value| value.parse::<f64>().ok())
}

fn output_lattice(lines: &[&str], alat: Option<f64>) -> Result<Lattice, Error> {
    if let Some(index) = last_card(lines, "CELL_PARAMETERS") {
        let unit = card_unit(lines.get(index).copied().unwrap_or(""));
        let scale = unit_scale(&unit, alat)?;
        let mut matrix = [[0.0; 3]; 3];
        let mut found = 0usize;
        let mut cursor = index + 1;
        while found < 3 && cursor < lines.len() {
            let line = lines.get(cursor).copied().unwrap_or("").trim();
            cursor += 1;
            if line.is_empty() {
                continue;
            }
            let values = three_numbers(line)
                .ok_or_else(|| Error::InvalidInput("bad pw.x output cell row".to_owned()))?;
            if let Some(row) = matrix.get_mut(found) {
                *row = values.map(|value| value * scale);
            }
            found += 1;
        }
        if found == 3 {
            return Lattice::new(matrix);
        }
    }
    let alat = alat.ok_or_else(|| {
        Error::InvalidInput("pw.x output has no lattice parameter (alat)".to_owned())
    })?;
    let index = lines
        .iter()
        .rposition(|line| line.to_ascii_lowercase().contains("crystal axes:"))
        .ok_or_else(|| Error::InvalidInput("pw.x output has no crystal axes".to_owned()))?;
    let mut matrix = [[0.0; 3]; 3];
    let mut found = 0usize;
    let mut cursor = index + 1;
    while found < 3 && cursor < lines.len() {
        let line = lines.get(cursor).copied().unwrap_or("");
        cursor += 1;
        if let Some(values) = parenthesised(line) {
            if let Some(row) = matrix.get_mut(found) {
                *row = values.map(|value| value * alat);
            }
            found += 1;
        }
    }
    if found < 3 {
        return Err(Error::InvalidInput(
            "pw.x output crystal axes are incomplete".to_owned(),
        ));
    }
    Lattice::new(matrix)
}

fn last_card(lines: &[&str], keyword: &str) -> Option<usize> {
    lines.iter().rposition(|line| {
        line.trim()
            .to_ascii_uppercase()
            .starts_with(&keyword.to_ascii_uppercase())
    })
}

fn atoms_per_cell(lines: &[&str]) -> Option<usize> {
    lines
        .iter()
        .find(|line| line.to_ascii_lowercase().contains("number of atoms/cell"))
        .and_then(|line| value_after_equals(line).map(|value| value as usize))
}

fn output_positions(
    lines: &[&str],
    nat: Option<usize>,
    lattice: &Lattice,
    alat: Option<f64>,
) -> Result<Vec<(String, [f64; 3])>, Error> {
    if let Some(index) = last_card(lines, "ATOMIC_POSITIONS") {
        let unit = card_unit(lines.get(index).copied().unwrap_or(""));
        let crystal = matches!(unit.as_str(), "crystal" | "crystal_sg");
        let scale = if crystal {
            1.0
        } else {
            unit_scale(&unit, alat)?
        };
        let mut positions = Vec::new();
        let mut cursor = index + 1;
        while cursor < lines.len() {
            let line = lines.get(cursor).copied().unwrap_or("").trim();
            cursor += 1;
            if line.is_empty() {
                break;
            }
            let mut fields = line.split_whitespace();
            let Some(symbol) = fields.next() else {
                continue;
            };
            let values: Vec<f64> = fields
                .take(3)
                .filter_map(|value| value.parse::<f64>().ok())
                .collect();
            if values.len() < 3 {
                break;
            }
            let cart = [values[0] * scale, values[1] * scale, values[2] * scale];
            let frac = if crystal {
                [values[0], values[1], values[2]]
            } else {
                lattice.to_frac(cart)
            };
            positions.push((symbol.to_owned(), frac));
            if nat.is_some_and(|count| positions.len() >= count) {
                break;
            }
        }
        if positions.len() == nat.unwrap_or(positions.len()) {
            return Ok(positions);
        }
    }
    let index = lines
        .iter()
        .rposition(|line| line.to_ascii_lowercase().contains("positions (alat units)"))
        .or_else(|| {
            lines
                .iter()
                .rposition(|line| line.to_ascii_lowercase().contains("site n.     atom"))
        })
        .ok_or_else(|| {
            Error::InvalidInput("pw.x output has no atomic position table".to_owned())
        })?;
    let alat = alat.ok_or_else(|| {
        Error::InvalidInput("pw.x output positions in alat units need alat".to_owned())
    })?;
    let mut positions = Vec::new();
    let mut cursor = index + 1;
    while cursor < lines.len() {
        let line = lines.get(cursor).copied().unwrap_or("");
        cursor += 1;
        let lower = line.to_ascii_lowercase();
        if !lower.contains("tau(") {
            if positions.is_empty() && line.trim().is_empty() {
                continue;
            }
            break;
        }
        let mut fields = line.split_whitespace();
        let _index = fields.next();
        let Some(symbol) = fields.next() else {
            break;
        };
        let values = parenthesised(line).ok_or_else(|| {
            Error::InvalidInput("pw.x output site line has no Cartesian vector".to_owned())
        })?;
        positions.push((
            symbol.to_owned(),
            lattice.to_frac(values.map(|value| value * alat)),
        ));
        if nat.is_some_and(|count| positions.len() >= count) {
            break;
        }
    }
    Ok(positions)
}

fn total_energy(lines: &[&str]) -> Option<f64> {
    lines
        .iter()
        .rev()
        .filter(|line| line.contains('!') && line.to_ascii_lowercase().contains("total energy"))
        .find_map(|line| value_after_equals(line))
}
