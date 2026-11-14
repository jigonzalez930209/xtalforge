//! CP2K output parser (Task 18.2.4).
//!
//! Reads the `ENERGY|` lines (Hartree to eV), the last
//! `ATOMIC COORDINATES` block per optimisation step, `ATOMIC FORCES in
//! [a.u.]` (Hartree/bohr to eV/Angstrom), the analytical stress tensor (GPa),
//! SCF and geometry markers, and the final `CP2K ... seconds` timing.

use crate::Error;
use crate::result::{CalculationResult, ElectronicStep, Snapshot};

use super::common::{
    BOHR_ANGSTROM, HARTREE_BOHR_EV_ANGSTROM, HARTREE_EV, build_structure, fields, number, numbers,
    rposition, stress_from_rows, warn,
};

type AtomRow = (String, [f64; 3]);

/// Parses a CP2K output into the common result model.
///
/// # Errors
///
/// Returns [`Error::Parse`] when the output has no `ENERGY|` line, a
/// coordinate block without a readable cell, or an inconsistent atom count.
pub fn parse(text: &str) -> Result<CalculationResult, Error> {
    let lines: Vec<&str> = text.lines().collect();
    let mut result = CalculationResult::new("cp2k");
    for line in &lines {
        if let Some(energy) = energy_line(line) {
            result.energies.push(energy);
            result.steps.push(ElectronicStep {
                energy,
                fermi: None,
                max_force: None,
            });
        }
    }
    if result.energies.is_empty() {
        return Err(Error::Parse("CP2K output has no ENERGY| line".to_owned()));
    }
    let cell = cell_vectors(&lines);
    let blocks = coordinate_blocks(&lines);
    if let Some(matrix) = cell {
        let lattice = xtalforge_core::Lattice::new(matrix)
            .map_err(|error| Error::Parse(format!("CP2K cell: {error}")))?;
        for (index, atoms) in blocks.iter().enumerate() {
            let elements = atoms
                .iter()
                .map(|(element, _)| element.clone())
                .collect::<Vec<String>>();
            let frac = atoms
                .iter()
                .map(|(_, point)| lattice.to_frac(*point))
                .collect::<Vec<_>>();
            let structure = build_structure("cp2k", matrix, &elements, &frac)?;
            result.trajectory.push(Snapshot {
                step: index as u64,
                structure,
                energy: result.energies.get(index).copied(),
                forces: None,
                stress: None,
            });
        }
    } else if !blocks.is_empty() {
        warn(
            &mut result,
            "cp2k_cell_missing",
            "CP2K coordinate blocks were found without CELL| vectors; no structure was built",
        );
    }
    result.forces = atomic_forces(&lines);
    result.stress = stress_tensor(&lines);
    if let Some(last) = result.trajectory.last_mut() {
        last.forces = result.forces.clone();
        last.stress = result.stress;
    }
    result.close_trajectory();
    convergence(text, &mut result);
    timings(&lines, &mut result);
    warnings(&lines, &mut result);
    result
        .scalars
        .insert("ionic_steps".to_owned(), result.trajectory.len() as f64);
    Ok(result)
}

fn energy_line(line: &str) -> Option<f64> {
    let trimmed = line.trim();
    if !trimmed.starts_with("ENERGY|") {
        return None;
    }
    let (_, rest) = trimmed.rsplit_once(':')?;
    numbers(rest)
        .first()
        .copied()
        .map(|value| value * HARTREE_EV)
}

fn cell_vectors(lines: &[&str]) -> Option<[[f64; 3]; 3]> {
    let mut rows: [Option<[f64; 3]>; 3] = [None; 3];
    for line in lines {
        let trimmed = line.trim();
        if !trimmed.starts_with("CELL| Vector") {
            continue;
        }
        let slot = if trimmed.contains("Vector a") {
            0
        } else if trimmed.contains("Vector b") {
            1
        } else if trimmed.contains("Vector c") {
            2
        } else {
            continue;
        };
        let scale = if trimmed.contains("bohr") {
            BOHR_ANGSTROM
        } else {
            1.0
        };
        let Some((_, rest)) = trimmed.split_once(':') else {
            continue;
        };
        let values = numbers(rest);
        let (Some(x), Some(y), Some(z)) = (values.first(), values.get(1), values.get(2)) else {
            continue;
        };
        if let Some(target) = rows.get_mut(slot) {
            *target = Some([x * scale, y * scale, z * scale]);
        }
    }
    Some([
        rows.first().copied().flatten()?,
        rows.get(1).copied().flatten()?,
        rows.get(2).copied().flatten()?,
    ])
}

fn coordinate_blocks(lines: &[&str]) -> Vec<Vec<AtomRow>> {
    let mut blocks = Vec::new();
    let mut index = 0usize;
    while index < lines.len() {
        let header = lines.get(index).copied().unwrap_or("").trim();
        if !header.to_ascii_uppercase().contains("ATOMIC COORDINATES") {
            index += 1;
            continue;
        }
        let scale = if header.to_ascii_lowercase().contains("a.u.") {
            BOHR_ANGSTROM
        } else {
            1.0
        };
        index += 1;
        let mut atoms = Vec::new();
        while index < lines.len() {
            let row = lines.get(index).copied().unwrap_or("").trim();
            if row.is_empty() {
                index += 1;
                continue;
            }
            let values = fields(row);
            let numbered = values
                .first()
                .and_then(|value| value.parse::<usize>().ok())
                .is_some();
            if !numbered {
                if atoms.is_empty() {
                    index += 1;
                    continue;
                }
                break;
            }
            let Some(element) = values.get(2) else {
                break;
            };
            let point = (
                values.get(3).and_then(|value| number(value)),
                values.get(4).and_then(|value| number(value)),
                values.get(5).and_then(|value| number(value)),
            );
            if let (Some(x), Some(y), Some(z)) = point {
                atoms.push(((*element).to_owned(), [x * scale, y * scale, z * scale]));
            }
            index += 1;
        }
        if !atoms.is_empty() {
            blocks.push(atoms);
        }
    }
    blocks
}

fn atomic_forces(lines: &[&str]) -> Option<Vec<[f64; 3]>> {
    let start = rposition(lines, |line| line.contains("ATOMIC FORCES in"))?;
    let scale = if lines.get(start).is_some_and(|line| line.contains("a.u.")) {
        HARTREE_BOHR_EV_ANGSTROM
    } else {
        1.0
    };
    let mut forces = Vec::new();
    for line in lines.get(start + 1..)? {
        let trimmed = line.trim();
        if trimmed.is_empty() {
            continue;
        }
        let values = fields(trimmed);
        if values
            .first()
            .and_then(|value| value.parse::<usize>().ok())
            .is_none()
        {
            if forces.is_empty() {
                continue;
            }
            break;
        }
        let point = (
            values.get(3).and_then(|value| number(value)),
            values.get(4).and_then(|value| number(value)),
            values.get(5).and_then(|value| number(value)),
        );
        if let (Some(x), Some(y), Some(z)) = point {
            forces.push([x * scale, y * scale, z * scale]);
        }
    }
    (!forces.is_empty()).then_some(forces)
}

fn stress_tensor(lines: &[&str]) -> Option<[[f64; 3]; 3]> {
    let start = rposition(lines, |line| {
        line.contains("STRESS|")
            && line
                .to_ascii_lowercase()
                .contains("analytical stress tensor")
    })?;
    let mut rows = Vec::new();
    for line in lines.get(start + 1..)? {
        if !line.contains("STRESS|") {
            if rows.is_empty() {
                continue;
            }
            break;
        }
        let (_, rest) = line.split_once('|')?;
        let values = numbers(rest);
        if values.len() < 3 {
            continue;
        }
        let Some(row) = values.get(..3) else {
            continue;
        };
        rows.push(row.to_vec());
        if rows.len() == 3 {
            break;
        }
    }
    stress_from_rows(&rows, 1.0)
}

fn convergence(text: &str, result: &mut CalculationResult) {
    let upper = text.to_ascii_uppercase();
    let converged = upper.contains("SCF RUN CONVERGED");
    let failed = upper.contains("SCF RUN NOT CONVERGED");
    if converged || failed {
        result.convergence.scf = converged;
        result.convergence.scf_reported = true;
    } else {
        warn(
            result,
            "cp2k_convergence_missing",
            "CP2K output has no SCF run convergence marker",
        );
    }
    if failed {
        warn(
            result,
            "cp2k_scf_not_converged",
            "CP2K reports that the SCF run did not converge",
        );
    }
    if upper.contains("GEOMETRY OPTIMIZATION COMPLETED") {
        result.convergence.ionic = true;
        result.convergence.ionic_reported = true;
    } else if upper.contains("MAXIMUM NUMBER OF OPTIMIZATION STEPS REACHED") {
        result.convergence.ionic = false;
        result.convergence.ionic_reported = true;
        warn(
            result,
            "cp2k_ionic_limit",
            "CP2K reached the maximum number of optimisation steps",
        );
    } else if upper.contains("GEO_OPT") || upper.contains("RUN TYPE: GEO_OPT") {
        warn(
            result,
            "cp2k_ionic_convergence_missing",
            "CP2K geometry optimisation has no completion marker",
        );
    }
}

fn timings(lines: &[&str], result: &mut CalculationResult) {
    for line in lines {
        let trimmed = line.trim();
        if trimmed.starts_with("CP2K")
            && trimmed.ends_with("seconds")
            && let Some(seconds) = numbers(trimmed).last().copied()
        {
            result.push_timing("cp2k", seconds);
        }
    }
}

fn warnings(lines: &[&str], result: &mut CalculationResult) {
    for line in lines {
        let trimmed = line.trim();
        if trimmed.contains("*** WARNING") {
            result.push_warning("cp2k_warning", trimmed);
        }
    }
}
