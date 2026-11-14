//! Quantum ESPRESSO `pw.x` text output parser (Task 18.2.3).
//!
//! The final structure comes from the Phase 13 `pwscf` reader. This module
//! adds the energies (Rydberg to eV), SCF steps, forces (Ry/bohr to
//! eV/Angstrom), stress (Ry/bohr^3 to GPa), Fermi level, convergence markers
//! and timings.

use xtalforge_io::pwscf::read_pwscf;

use crate::Error;
use crate::result::{CalculationResult, ElectronicStep, Snapshot};

use super::common::{
    RYDBERG_BOHR_EV_ANGSTROM, RYDBERG_BOHR3_GPA, RYDBERG_EV, number, numbers, rposition,
    value_after, warn,
};

/// Parses a `pw.x` text output into the common result model.
///
/// # Errors
///
/// Returns [`Error::Parse`] when the cell or positions cannot be read, or when
/// the output ends before the final total energy.
pub fn parse(text: &str) -> Result<CalculationResult, Error> {
    let structure =
        read_pwscf(text).map_err(|error| Error::Parse(format!("pw.x output: {error}")))?;
    let mut result = CalculationResult::new("quantum_espresso");
    let lines: Vec<&str> = text.lines().collect();
    for line in &lines {
        let trimmed = line.trim();
        if trimmed.starts_with('!')
            && trimmed.contains("total energy")
            && let Some(ry) = value_after(trimmed, "=")
        {
            result.energies.push(ry * RYDBERG_EV);
            continue;
        }
        if trimmed.contains("total energy")
            && trimmed.contains('=')
            && let Some(ry) = value_after(trimmed, "=")
        {
            result.steps.push(ElectronicStep {
                energy: ry * RYDBERG_EV,
                fermi: None,
                max_force: None,
            });
        }
    }
    if result.energies.is_empty() {
        return Err(Error::Parse(
            "pw.x output has no final total energy".to_owned(),
        ));
    }
    result.fermi = lines
        .iter()
        .rev()
        .filter(|line| line.contains("highest occupied, lowest unoccupied level"))
        .find_map(|line| value_after(line, ":"))
        .or_else(|| {
            lines
                .iter()
                .rev()
                .filter(|line| line.contains("the Fermi energy is"))
                .find_map(|line| value_after(line, "is"))
        });
    result.forces = last_forces(&lines);
    result.stress = last_stress(&lines);
    result.magnetic_moments = Vec::new();
    for line in &lines {
        if line.contains("total magnetization")
            && let Some(value) = value_after(line, "=")
        {
            result
                .scalars
                .insert("total_magnetization".to_owned(), value);
        } else if line.contains("absolute magnetization")
            && let Some(value) = value_after(line, "=")
        {
            result
                .scalars
                .insert("absolute_magnetization".to_owned(), value);
        }
    }
    convergence(&lines, &mut result);
    timings(&lines, &mut result);
    warnings(&lines, &mut result);
    let final_energy = result.energies.last().copied();
    result.trajectory.push(Snapshot {
        step: 0,
        structure: structure.clone(),
        energy: final_energy,
        forces: result.forces.clone(),
        stress: result.stress,
    });
    result.structure = Some(structure);
    Ok(result)
}

fn last_forces(lines: &[&str]) -> Option<Vec<[f64; 3]>> {
    let start = rposition(lines, |line| line.contains("Forces acting on atoms"))?;
    let mut forces = Vec::new();
    for line in lines.get(start + 1..)? {
        if line.contains("force =") {
            let (_, rest) = line.split_once("force =")?;
            let values = numbers(rest);
            if let (Some(x), Some(y), Some(z)) = (values.first(), values.get(1), values.get(2)) {
                forces.push([
                    x * RYDBERG_BOHR_EV_ANGSTROM,
                    y * RYDBERG_BOHR_EV_ANGSTROM,
                    z * RYDBERG_BOHR_EV_ANGSTROM,
                ]);
            }
            continue;
        }
        if !forces.is_empty() {
            break;
        }
    }
    (!forces.is_empty()).then_some(forces)
}

fn last_stress(lines: &[&str]) -> Option<[[f64; 3]; 3]> {
    let start = rposition(lines, |line| {
        line.contains("total   stress") || line.contains("total  stress")
    })?;
    let mut rows = Vec::new();
    for line in lines.get(start + 1..)? {
        let values = numbers(line);
        if values.len() < 3 {
            if rows.is_empty() {
                continue;
            }
            break;
        }
        let Some(row) = values.get(..3) else {
            continue;
        };
        rows.push(row.to_vec());
        if rows.len() == 3 {
            break;
        }
    }
    super::common::stress_from_rows(&rows, RYDBERG_BOHR3_GPA)
}

fn convergence(lines: &[&str], result: &mut CalculationResult) {
    let scf = lines
        .iter()
        .any(|line| line.contains("convergence has been achieved"));
    let ionic = lines.iter().any(|line| {
        line.contains("bfgs converged")
            || line.contains("geometry optimization completed")
            || line.contains("Optimization completed")
    });
    result.convergence.scf = scf;
    result.convergence.scf_reported = lines.iter().any(|line| {
        line.contains("convergence has been achieved") || line.contains("not converged")
    });
    result.convergence.ionic = ionic;
    result.convergence.ionic_reported = ionic;
    if !result.convergence.scf_reported {
        warn(
            result,
            "qe_convergence_missing",
            "pw.x output has no SCF convergence marker",
        );
    }
}

fn timings(lines: &[&str], result: &mut CalculationResult) {
    for line in lines {
        if line.contains("total cpu time spent up to now is") {
            if let Some(seconds) = value_after(line, "is") {
                result.push_timing("total_cpu", seconds);
            }
            continue;
        }
        if !(line.contains("s CPU") || line.contains("s WALL")) {
            continue;
        }
        let Some((stage, rest)) = line.split_once(':') else {
            continue;
        };
        let stage = stage.trim();
        if stage.is_empty() || stage.len() > 40 {
            continue;
        }
        let seconds = rest
            .trim_start()
            .split_once('s')
            .and_then(|(value, _)| number(value))
            .or_else(|| numbers(rest).first().copied());
        if let Some(seconds) = seconds {
            result.push_timing(stage, seconds);
        }
    }
}

fn warnings(lines: &[&str], result: &mut CalculationResult) {
    for line in lines {
        let trimmed = line.trim();
        if trimmed.starts_with("Warning:") {
            result.push_warning("qe_warning", trimmed);
        } else if trimmed.contains("Error in routine") {
            result.push_warning("qe_error", trimmed);
        }
    }
}
