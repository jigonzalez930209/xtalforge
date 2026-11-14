//! VASP `OUTCAR` result parser (Task 18.2.2).
//!
//! Reuses the Phase 13 `OUTCAR` reader for the trajectory (positions, forces,
//! stress and energies) and adds magnetic moments, the Fermi level, total
//! charge, convergence markers, timings and warnings.

use xtalforge_io::trajectory::Trajectory;
use xtalforge_io::vasp::read_outcar;

use crate::Error;
use crate::result::{CalculationResult, ElectronicStep, Snapshot};

use super::common::{max_force, number, value_after, warn};

/// Parses an `OUTCAR` into the common result model.
///
/// # Errors
///
/// Returns [`Error::Parse`] when the file has no `TITEL`, lattice or
/// `POSITION`/`TOTAL-FORCE` block, or when atom counts disagree.
pub fn parse(text: &str) -> Result<CalculationResult, Error> {
    let trajectory = read_outcar(text).map_err(|error| Error::Parse(format!("OUTCAR: {error}")))?;
    let mut result = CalculationResult::new("vasp");
    snapshots(&trajectory, &mut result)?;
    let lines: Vec<&str> = text.lines().collect();
    result.fermi = lines
        .iter()
        .rev()
        .filter(|line| line.contains("E-fermi"))
        .find_map(|line| value_after(line, "E-fermi"))
        .or_else(|| {
            lines
                .iter()
                .rev()
                .filter(|line| line.contains("E-fermi"))
                .find_map(|line| value_after(line, ":"))
        });
    result.magnetic_moments = last_magnetization(&lines).unwrap_or_default();
    let sites = result
        .structure
        .as_ref()
        .map_or(0, |structure| structure.sites().len());
    if sites != 0 && !result.magnetic_moments.is_empty() && result.magnetic_moments.len() != sites {
        warn(
            &mut result,
            "vasp_magnetization_mismatch",
            "OUTCAR magnetization rows do not match the number of sites",
        );
        result.magnetic_moments.clear();
    }
    if let Some(charge) = lines
        .iter()
        .rev()
        .filter(|line| line.contains("total charge"))
        .find_map(|line| value_after(line, "="))
    {
        result.scalars.insert("total_charge".to_owned(), charge);
    }
    convergence(&lines, &mut result);
    timings(&lines, &mut result);
    warnings(&lines, &mut result);
    result.close_trajectory();
    Ok(result)
}

fn snapshots(trajectory: &Trajectory, result: &mut CalculationResult) -> Result<(), Error> {
    for (index, frame) in trajectory.frames.iter().enumerate() {
        let structure = trajectory
            .structure(index)
            .map_err(|error| Error::Parse(format!("OUTCAR step {index}: {error}")))?;
        result.trajectory.push(Snapshot {
            step: frame.step.unwrap_or(index as u64),
            structure,
            energy: frame.energy,
            forces: frame.forces.clone(),
            stress: frame.stress,
        });
        if let Some(energy) = frame.energy {
            result.energies.push(energy);
            result.steps.push(ElectronicStep {
                energy,
                fermi: None,
                max_force: frame.forces.as_deref().map(max_force),
            });
        }
    }
    Ok(())
}

fn last_magnetization(lines: &[&str]) -> Option<Vec<f64>> {
    let start = lines
        .iter()
        .rposition(|line| line.to_ascii_lowercase().contains("magnetization (x)"))?;
    let mut moments = Vec::new();
    for line in lines.get(start + 1..)? {
        let fields: Vec<&str> = line.split_whitespace().collect();
        let Some(first) = fields.first() else {
            if moments.is_empty() {
                continue;
            }
            break;
        };
        if first.parse::<usize>().is_err() {
            if moments.is_empty() {
                continue;
            }
            break;
        }
        let Some(value) = fields.get(1).and_then(|value| number(value)) else {
            if moments.is_empty() {
                continue;
            }
            break;
        };
        moments.push(value);
    }
    (!moments.is_empty()).then_some(moments)
}

fn convergence(lines: &[&str], result: &mut CalculationResult) {
    let reached = lines
        .iter()
        .any(|line| line.contains("reached required accuracy"));
    let ediff = lines.iter().any(|line| {
        line.contains("EDIFF is reached") || line.contains("aborting loop because EDIFF")
    });
    result.convergence.scf = reached || ediff;
    result.convergence.scf_reported = reached || ediff;
    result.convergence.ionic = reached;
    result.convergence.ionic_reported = reached;
    if !result.convergence.scf_reported {
        warn(
            result,
            "vasp_convergence_missing",
            "OUTCAR has no EDIFF or accuracy convergence marker",
        );
    }
}

fn timings(lines: &[&str], result: &mut CalculationResult) {
    for line in lines {
        let stage = if line.contains("Total CPU time used") {
            Some("cpu")
        } else if line.contains("Elapsed time") {
            Some("wall")
        } else {
            None
        };
        if let Some(stage) = stage
            && let Some(seconds) = value_after(line, ":")
        {
            result.push_timing(stage, seconds);
        }
    }
}

fn warnings(lines: &[&str], result: &mut CalculationResult) {
    for line in lines {
        let trimmed = line.trim();
        if trimmed.starts_with("WARNING") {
            result.push_warning("vasp_warning", trimmed);
        } else if trimmed.contains("ZBRENT") {
            result.push_warning("vasp_zbrent", trimmed);
        }
    }
}
