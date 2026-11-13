//! VASP `OUTCAR` reader (Task 13.2.4).
//!
//! Reads the ionic steps of a VASP run: Cartesian positions and forces
//! (`POSITION ... TOTAL-FORCE`), the total energy (`free energy TOTEN`), the
//! stress (`in kB`, converted to GPa) and the direct lattice vectors. Frames
//! use the existing [`crate::trajectory::Trajectory`] model; the structure
//! reader returns the last frame.

#![allow(clippy::indexing_slicing)] // fixed-format parser, lengths checked

use std::collections::BTreeMap;

use xtalforge_core::elements::parse_symbol;
use xtalforge_core::{Error, Lattice, Structure};

use super::super::trajectory::{Frame, Trajectory};

type RawStep = (Vec<[f64; 3]>, Vec<[f64; 3]>);

/// Reads an `OUTCAR` trajectory.
///
/// # Errors
///
/// Returns [`Error::InvalidInput`] when the file has no `TITEL`, no lattice or
/// no `POSITION ... TOTAL-FORCE` block, or when the atom counts disagree.
pub fn read_outcar(text: &str) -> Result<Trajectory, Error> {
    let lines: Vec<&str> = text.lines().collect();
    let mut titles: Vec<String> = Vec::new();
    let mut counts: Vec<usize> = Vec::new();
    let mut lattices: Vec<[[f64; 3]; 3]> = Vec::new();
    let mut raw_frames: Vec<RawStep> = Vec::new();
    let mut energies: Vec<f64> = Vec::new();
    let mut stresses: Vec<[[f64; 3]; 3]> = Vec::new();
    let mut index = 0usize;
    while index < lines.len() {
        let line = lines.get(index).copied().unwrap_or("");
        let trimmed = line.trim();
        if trimmed.starts_with("TITEL")
            && let Some(symbol) = title_element(trimmed)
        {
            titles.push(symbol);
            index += 1;
            continue;
        }
        if trimmed.starts_with("ions per type") {
            counts = integers_after_equals(trimmed);
            index += 1;
            continue;
        }
        if trimmed.starts_with("direct lattice vectors")
            && let Some((matrix, next)) = read_lattice(&lines, index + 1)
        {
            lattices.push(matrix);
            index = next;
            continue;
        }
        if trimmed.contains("POSITION") && trimmed.contains("TOTAL-FORCE") {
            let (frame, next) = read_force_block(&lines, index + 1)?;
            raw_frames.push(frame);
            index = next;
            continue;
        }
        let compact = trimmed.split_whitespace().collect::<Vec<_>>().join(" ");
        if compact.contains("free energy TOTEN")
            && let Some(value) = value_after_equals(trimmed)
        {
            energies.push(value);
        }
        if trimmed.starts_with("in kB")
            && let Some(values) = six_numbers(trimmed.get(5..).unwrap_or(""))
        {
            stresses.push(stress_matrix(&values));
        }
        index += 1;
    }
    if raw_frames.is_empty() {
        return Err(Error::InvalidInput(
            "OUTCAR has no POSITION / TOTAL-FORCE block".to_owned(),
        ));
    }
    if lattices.is_empty() {
        return Err(Error::InvalidInput(
            "OUTCAR has no direct lattice vectors".to_owned(),
        ));
    }
    let elements = elements_from_titles(&titles, &counts)?;
    let first_lattice =
        lattices
            .first()
            .copied()
            .unwrap_or([[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]]);
    let atom_count = raw_frames
        .first()
        .map_or(0, |(positions, _)| positions.len());
    if atom_count != elements.len() {
        return Err(Error::InvalidInput(format!(
            "OUTCAR has {atom_count} atoms per step but {} element entries",
            elements.len()
        )));
    }
    let variable = lattices.len() == raw_frames.len()
        && lattices.windows(2).any(|pair| pair.first() != pair.get(1));
    let fixed = if variable { None } else { Some(first_lattice) };
    let mut frames = Vec::with_capacity(raw_frames.len());
    for (step, (positions, forces)) in raw_frames.iter().enumerate() {
        let matrix = if variable {
            lattices.get(step).copied().unwrap_or(first_lattice)
        } else {
            first_lattice
        };
        let lattice = Lattice::new(matrix)?;
        if positions.len() != elements.len() || forces.len() != elements.len() {
            return Err(Error::InvalidInput(format!(
                "OUTCAR step {step} has {} positions and {} forces for {} atoms",
                positions.len(),
                forces.len(),
                elements.len()
            )));
        }
        let mut info = BTreeMap::new();
        info.insert("source".to_owned(), "outcar".to_owned());
        frames.push(Frame {
            frac: positions
                .iter()
                .map(|point| lattice.to_frac(*point))
                .collect(),
            lattice: variable.then_some(matrix),
            energy: energies.get(step).copied(),
            step: Some(step as u64),
            info,
            forces: Some(forces.clone()),
            stress: stresses.get(step).copied(),
        });
    }
    Ok(Trajectory {
        elements,
        lattice: fixed,
        frames,
    })
}

/// Reads the last ionic step of an `OUTCAR` as a structure.
///
/// # Errors
///
/// Returns the errors of [`read_outcar`] or of the frame-to-structure
/// conversion.
pub fn read_outcar_structure(text: &str) -> Result<Structure, Error> {
    let trajectory = read_outcar(text)?;
    let last = trajectory
        .frames
        .len()
        .checked_sub(1)
        .ok_or_else(|| Error::InvalidInput("OUTCAR has no frames".to_owned()))?;
    trajectory.structure(last)
}

fn title_element(line: &str) -> Option<String> {
    let (_, rest) = line.split_once('=')?;
    let fields: Vec<&str> = rest.split_whitespace().collect();
    fields
        .iter()
        .skip(1)
        .copied()
        .chain(fields.first().copied())
        .find_map(|token| parse_symbol(token).ok())
}

fn integers_after_equals(line: &str) -> Vec<usize> {
    let Some((_, rest)) = line.split_once('=') else {
        return Vec::new();
    };
    rest.split_whitespace()
        .filter_map(|value| value.parse().ok())
        .collect()
}

fn value_after_equals(line: &str) -> Option<f64> {
    let (_, rest) = line.split_once('=')?;
    rest.split_whitespace()
        .find_map(|value| value.parse::<f64>().ok())
}

fn read_lattice(lines: &[&str], start: usize) -> Option<([[f64; 3]; 3], usize)> {
    let mut matrix = [[0.0; 3]; 3];
    let mut found = 0usize;
    let mut index = start;
    while found < 3 && index < lines.len() {
        let line = lines.get(index).copied().unwrap_or("").trim();
        index += 1;
        if line.is_empty() {
            continue;
        }
        if let Some(values) = three_numbers(line) {
            if let Some(row) = matrix.get_mut(found) {
                *row = values;
            }
            found += 1;
        } else {
            break;
        }
    }
    (found == 3).then_some((matrix, index))
}

fn read_force_block(lines: &[&str], start: usize) -> Result<(RawStep, usize), Error> {
    let mut positions = Vec::new();
    let mut forces = Vec::new();
    let mut index = start;
    while index < lines.len() {
        let line = lines.get(index).copied().unwrap_or("").trim();
        index += 1;
        if line.is_empty() {
            continue;
        }
        if line.starts_with('-') {
            if positions.is_empty() {
                continue;
            }
            break;
        }
        let Some(values) = six_numbers(line) else {
            if positions.is_empty() {
                continue;
            }
            break;
        };
        positions.push([values[0], values[1], values[2]]);
        forces.push([values[3], values[4], values[5]]);
    }
    if positions.is_empty() {
        return Err(Error::InvalidInput(
            "OUTCAR force block has no atom rows".to_owned(),
        ));
    }
    Ok(((positions, forces), index))
}

fn three_numbers(line: &str) -> Option<[f64; 3]> {
    let mut values = line
        .split_whitespace()
        .filter_map(|value| value.parse::<f64>().ok());
    Some([values.next()?, values.next()?, values.next()?])
}

fn six_numbers(line: &str) -> Option<[f64; 6]> {
    let values: Vec<f64> = line
        .split_whitespace()
        .filter_map(|value| value.parse::<f64>().ok())
        .collect();
    if values.len() < 6 {
        return None;
    }
    Some([
        values[0], values[1], values[2], values[3], values[4], values[5],
    ])
}

fn stress_matrix(values: &[f64; 6]) -> [[f64; 3]; 3] {
    let kbar = 0.1;
    let xx = values[0] * kbar;
    let yy = values[1] * kbar;
    let zz = values[2] * kbar;
    let xy = values[3] * kbar;
    let yz = values[4] * kbar;
    let zx = values[5] * kbar;
    [[xx, xy, zx], [xy, yy, yz], [zx, yz, zz]]
}

fn elements_from_titles(titles: &[String], counts: &[usize]) -> Result<Vec<String>, Error> {
    if titles.is_empty() {
        return Err(Error::InvalidInput(
            "OUTCAR has no TITEL element lines".to_owned(),
        ));
    }
    if counts.is_empty() {
        return Ok(titles.to_vec());
    }
    if counts.len() != titles.len() {
        return Err(Error::InvalidInput(format!(
            "OUTCAR has {} TITEL lines but {} ions-per-type counts",
            titles.len(),
            counts.len()
        )));
    }
    let mut elements = Vec::new();
    for (title, count) in titles.iter().zip(counts.iter()) {
        for _ in 0..*count {
            elements.push(title.clone());
        }
    }
    Ok(elements)
}
