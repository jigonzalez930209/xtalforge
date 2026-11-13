//! Quantum ESPRESSO `pw.x` input deck reader (Task 13.2.1).

#![allow(clippy::indexing_slicing)] // fixed-format parser, lengths checked

use serde_json::Value;
use xtalforge_core::elements::parse_symbol;
use xtalforge_core::{Error, Lattice, Site, Structure, StructureBuilder};

use super::{BOHR, card_unit, strip_comments, three_numbers, unit_scale, variable};

pub(super) fn read(text: &str) -> Result<Structure, Error> {
    let cleaned = strip_comments(text);
    let lines: Vec<&str> = cleaned.lines().collect();
    let alat = alat_angstrom(text);
    let nat = variable(text, "nat").and_then(|value| value.parse::<usize>().ok());
    let ntyp = variable(text, "ntyp").and_then(|value| value.parse::<usize>().ok());
    let lattice = read_lattice(&lines, text, alat)?;
    let species = read_species(&lines, ntyp)?;
    let (positions, count) = read_positions(&lines, nat, &lattice, alat)?;
    if positions.len() != count {
        return Err(Error::InvalidInput(format!(
            "pw.x declares {count} atoms but ATOMIC_POSITIONS has {}",
            positions.len()
        )));
    }
    let mut builder = StructureBuilder::new("pwscf", lattice);
    for (index, (symbol, frac)) in positions.iter().enumerate() {
        let element = species
            .iter()
            .find(|(label, _)| label.eq_ignore_ascii_case(symbol))
            .map(|(_, element)| element.clone())
            .ok_or_else(|| {
                Error::InvalidInput(format!(
                    "pw.x ATOMIC_POSITIONS uses species {symbol:?} that ATOMIC_SPECIES does not define"
                ))
            })?;
        builder = builder.site(Site::new(
            format!("{element}{}", index + 1),
            &element,
            *frac,
        )?);
    }
    if let Some(ibrav) = variable(text, "ibrav") {
        builder = builder.metadata("pwscf_ibrav", Value::String(ibrav));
    }
    if let Some(calculation) = variable(text, "calculation") {
        builder = builder.metadata("pwscf_calculation", Value::String(calculation));
    }
    if let Some(prefix) = variable(text, "prefix") {
        builder = builder.metadata("pwscf_prefix", Value::String(prefix));
    }
    builder.build()
}

fn alat_angstrom(text: &str) -> Option<f64> {
    if let Some(value) = variable(text, "A").and_then(|value| value.parse::<f64>().ok()) {
        return Some(value);
    }
    variable(text, "celldm(1)")
        .and_then(|value| value.parse::<f64>().ok())
        .map(|value| value * BOHR)
}

fn read_lattice(lines: &[&str], text: &str, alat: Option<f64>) -> Result<Lattice, Error> {
    if let Some(index) = card_index(lines, "CELL_PARAMETERS") {
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
                .ok_or_else(|| Error::InvalidInput("bad pw.x CELL_PARAMETERS row".to_owned()))?;
            if let Some(row) = matrix.get_mut(found) {
                *row = values.map(|value| value * scale);
            }
            found += 1;
        }
        if found < 3 {
            return Err(Error::InvalidInput(
                "pw.x CELL_PARAMETERS needs three rows".to_owned(),
            ));
        }
        return Lattice::new(matrix);
    }
    let ibrav = variable(text, "ibrav")
        .and_then(|value| value.parse::<i32>().ok())
        .unwrap_or(0);
    if ibrav == 0 {
        return Err(Error::InvalidInput(
            "pw.x input has no CELL_PARAMETERS and no supported ibrav".to_owned(),
        ));
    }
    let a = alat.ok_or_else(|| {
        Error::InvalidInput("pw.x ibrav needs celldm(1) or A in &system".to_owned())
    })?;
    let ratio = |name: &str| {
        variable(text, name)
            .and_then(|value| value.parse::<f64>().ok())
            .unwrap_or(1.0)
    };
    let matrix = match ibrav {
        1 => [[a, 0.0, 0.0], [0.0, a, 0.0], [0.0, 0.0, a]],
        2 => {
            let half = a / 2.0;
            [[-half, 0.0, half], [0.0, half, half], [-half, half, 0.0]]
        }
        3 => {
            let half = a / 2.0;
            [
                [half, half, half],
                [-half, half, half],
                [-half, -half, half],
            ]
        }
        6 => {
            let c = a * ratio("celldm(3)");
            [[a, 0.0, 0.0], [0.0, a, 0.0], [0.0, 0.0, c]]
        }
        8 => {
            let b = a * ratio("celldm(2)");
            let c = a * ratio("celldm(3)");
            [[a, 0.0, 0.0], [0.0, b, 0.0], [0.0, 0.0, c]]
        }
        other => {
            return Err(Error::InvalidInput(format!(
                "pw.x ibrav={other} is not supported; write CELL_PARAMETERS with ibrav=0"
            )));
        }
    };
    Lattice::new(matrix)
}

fn card_index(lines: &[&str], keyword: &str) -> Option<usize> {
    lines.iter().position(|line| {
        line.trim()
            .to_ascii_uppercase()
            .starts_with(&keyword.to_ascii_uppercase())
    })
}

fn read_species(lines: &[&str], ntyp: Option<usize>) -> Result<Vec<(String, String)>, Error> {
    let index = card_index(lines, "ATOMIC_SPECIES")
        .ok_or_else(|| Error::InvalidInput("pw.x input has no ATOMIC_SPECIES card".to_owned()))?;
    let mut species = Vec::new();
    let mut cursor = index + 1;
    while cursor < lines.len() {
        let line = lines.get(cursor).copied().unwrap_or("").trim();
        cursor += 1;
        if line.is_empty() {
            continue;
        }
        if line.starts_with('&') || line.starts_with('/') {
            break;
        }
        let mut fields = line.split_whitespace();
        let Some(label) = fields.next() else {
            continue;
        };
        if matches!(
            label.to_ascii_uppercase().as_str(),
            "CELL_PARAMETERS" | "ATOMIC_POSITIONS" | "K_POINTS" | "ATOMIC_SPECIES"
        ) {
            break;
        }
        if fields.next().is_none() {
            break;
        }
        let element = parse_symbol(label)?;
        species.push((label.to_owned(), element));
        if ntyp.is_some_and(|count| species.len() >= count) {
            break;
        }
    }
    if species.is_empty() {
        return Err(Error::InvalidInput(
            "pw.x ATOMIC_SPECIES defines no species".to_owned(),
        ));
    }
    Ok(species)
}

type Position = (String, [f64; 3]);

fn read_positions(
    lines: &[&str],
    nat: Option<usize>,
    lattice: &Lattice,
    alat: Option<f64>,
) -> Result<(Vec<Position>, usize), Error> {
    let index = card_index(lines, "ATOMIC_POSITIONS")
        .ok_or_else(|| Error::InvalidInput("pw.x input has no ATOMIC_POSITIONS card".to_owned()))?;
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
            continue;
        }
        if line.starts_with('&') || line.starts_with('/') || line.starts_with("K_POINTS") {
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
            return Err(Error::InvalidInput(
                "a pw.x ATOMIC_POSITIONS line needs three coordinates".to_owned(),
            ));
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
    let count = nat.unwrap_or(positions.len());
    Ok((positions, count))
}
