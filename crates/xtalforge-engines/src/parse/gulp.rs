//! GULP output parser (Task 18.2.4).
//!
//! Reads the lattice energy in eV, the final Cartesian coordinates and cell
//! parameters, forces (or gradients, negated) in eV/Angstrom, the stress
//! tensor in GPa, SCF/optimisation markers and CPU timings.

use xtalforge_core::Lattice;

use crate::Error;
use crate::result::{CalculationResult, Snapshot};

use super::common::{
    build_structure, fields, number, numbers, rposition, stress_from_rows, value_after, warn,
};

/// Parses a GULP output into the common result model.
///
/// # Errors
///
/// Returns [`Error::Parse`] when the output has no lattice energy or no final
/// coordinates and cell.
pub fn parse(text: &str) -> Result<CalculationResult, Error> {
    let lines: Vec<&str> = text.lines().collect();
    let mut result = CalculationResult::new("gulp");
    if let Some(energy) = energy(&lines) {
        result.energies.push(energy);
    } else {
        return Err(Error::Parse("GULP output has no lattice energy".to_owned()));
    }
    let cell = final_cell(&lines);
    let atoms = final_atoms(&lines);
    if let (Some(parameters), Some(atoms)) = (cell, atoms.as_ref()) {
        let [a, b, c, alpha, beta, gamma] = parameters;
        let lattice = Lattice::from_parameters(a, b, c, alpha, beta, gamma)
            .map_err(|error| Error::Parse(format!("GULP cell: {error}")))?;
        let matrix = lattice.matrix();
        let elements = atoms
            .iter()
            .map(|(element, _)| element.clone())
            .collect::<Vec<String>>();
        let frac = atoms
            .iter()
            .map(|(_, cart)| lattice.to_frac(*cart))
            .collect::<Vec<_>>();
        let structure = build_structure("gulp", matrix, &elements, &frac)?;
        result.structure = Some(structure);
    } else if atoms.is_some() {
        warn(
            &mut result,
            "gulp_cell_missing",
            "GULP prints coordinates without final cell parameters; no structure was built",
        );
    }
    result.forces = forces(&lines);
    result.stress = stress(&lines);
    let energy = result.energies.last().copied();
    if let Some(structure) = result.structure.clone() {
        result.trajectory.push(Snapshot {
            step: 0,
            structure,
            energy,
            forces: result.forces.clone(),
            stress: result.stress,
        });
    }
    convergence(&lines, &mut result);
    timings(&lines, &mut result);
    warnings(&lines, &mut result);
    Ok(result)
}

fn energy(lines: &[&str]) -> Option<f64> {
    lines
        .iter()
        .rev()
        .filter(|line| {
            let lower = line.to_ascii_lowercase();
            lower.contains("lattice energy") || lower.contains("total energy")
        })
        .find_map(|line| value_after(line, "="))
}

fn final_cell(lines: &[&str]) -> Option<[f64; 6]> {
    let start = rposition(lines, |line| {
        line.to_ascii_lowercase().contains("final cell parameters")
    })?;
    let mut collected = Vec::new();
    for line in lines.get(start + 1..)?.iter().take(8) {
        collected.extend(numbers(line));
        if collected.len() >= 6 {
            break;
        }
    }
    if collected.len() < 6 {
        return None;
    }
    Some([
        *collected.first()?,
        *collected.get(1)?,
        *collected.get(2)?,
        *collected.get(3)?,
        *collected.get(4)?,
        *collected.get(5)?,
    ])
}

fn final_atoms(lines: &[&str]) -> Option<Vec<(String, [f64; 3])>> {
    let start = rposition(lines, |line| {
        line.to_ascii_lowercase()
            .contains("final coordinates (angstroms)")
    })?;
    let mut atoms = Vec::new();
    for line in lines.get(start + 1..)? {
        let values = fields(line);
        let Some(first) = values.first() else {
            if atoms.is_empty() {
                continue;
            }
            break;
        };
        if first.parse::<usize>().is_err() {
            if atoms.is_empty() {
                continue;
            }
            break;
        }
        let Some(element) = values.get(1) else {
            break;
        };
        let coordinates = values
            .iter()
            .skip(2)
            .filter_map(|value| number(value))
            .collect::<Vec<f64>>();
        let (Some(x), Some(y), Some(z)) =
            (coordinates.first(), coordinates.get(1), coordinates.get(2))
        else {
            break;
        };
        atoms.push(((*element).to_owned(), [*x, *y, *z]));
    }
    (!atoms.is_empty()).then_some(atoms)
}

fn forces(lines: &[&str]) -> Option<Vec<[f64; 3]>> {
    let start = rposition(lines, |line| {
        let lower = line.to_ascii_lowercase();
        lower.contains("forces (ev/angstrom") || lower.contains("derivatives (ev/angstrom")
    })?;
    let negate = lines
        .get(start)
        .is_some_and(|line| line.to_ascii_lowercase().contains("derivatives"));
    let mut forces = Vec::new();
    for line in lines.get(start + 1..)? {
        let values = fields(line);
        let Some(first) = values.first() else {
            if forces.is_empty() {
                continue;
            }
            break;
        };
        if first.parse::<usize>().is_err() {
            if forces.is_empty() {
                continue;
            }
            break;
        }
        let components = values
            .iter()
            .skip(2)
            .filter_map(|value| number(value))
            .collect::<Vec<f64>>();
        let (Some(x), Some(y), Some(z)) =
            (components.first(), components.get(1), components.get(2))
        else {
            break;
        };
        let sign = if negate { -1.0 } else { 1.0 };
        forces.push([x * sign, y * sign, z * sign]);
    }
    (!forces.is_empty()).then_some(forces)
}

fn stress(lines: &[&str]) -> Option<[[f64; 3]; 3]> {
    let start = rposition(lines, |line| {
        let lower = line.to_ascii_lowercase();
        lower.contains("stress") && lower.contains("gpa")
    })?;
    let mut rows = Vec::new();
    for line in lines.get(start + 1..)? {
        let values = numbers(line);
        if values.is_empty() {
            if rows.is_empty() {
                continue;
            }
            break;
        }
        rows.push(values);
        if rows.len() == 3 {
            break;
        }
    }
    if let Some(row) = rows.first()
        && row.len() >= 6
    {
        let value = |index: usize| row.get(index).copied().unwrap_or(0.0);
        return Some([
            [value(0), value(3), value(4)],
            [value(3), value(1), value(5)],
            [value(4), value(5), value(2)],
        ]);
    }
    stress_from_rows(&rows, 1.0)
}

fn convergence(lines: &[&str], result: &mut CalculationResult) {
    let upper = lines
        .iter()
        .map(|line| line.to_ascii_uppercase())
        .collect::<Vec<String>>();
    let has = |needle: &str| upper.iter().any(|line| line.contains(needle));
    let scf = has("SCF CONVERGED");
    let ionic = has("OPTIMISATION COMPLETED") || has("OPTIMIZATION COMPLETED");
    result.convergence.scf = scf;
    result.convergence.scf_reported = has("SCF CONVERGED") || has("SCF NOT CONVERGED");
    result.convergence.ionic = ionic;
    result.convergence.ionic_reported = ionic;
    if !result.convergence.scf_reported {
        warn(
            result,
            "gulp_convergence_missing",
            "GULP output has no SCF convergence marker",
        );
    }
}

fn timings(lines: &[&str], result: &mut CalculationResult) {
    for line in lines {
        let lower = line.to_ascii_lowercase();
        let stage = if lower.contains("total cpu time") {
            Some("cpu")
        } else if lower.contains("total time") {
            Some("wall")
        } else {
            None
        };
        if let Some(stage) = stage
            && let Some(seconds) = value_after(line, "=")
        {
            result.push_timing(stage, seconds);
        }
    }
}

fn warnings(lines: &[&str], result: &mut CalculationResult) {
    for line in lines {
        let trimmed = line.trim();
        if trimmed.starts_with("WARNING") {
            result.push_warning("gulp_warning", trimmed);
        } else if trimmed.starts_with("ERROR") {
            result.push_warning("gulp_error", trimmed);
        }
    }
}
