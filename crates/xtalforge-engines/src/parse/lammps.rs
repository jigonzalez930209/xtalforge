//! LAMMPS log and dump parsers (Task 18.2.4).
//!
//! The log reader follows the thermo table headers, converts energies with
//! the active `units` style and records pressure, timings and error lines. The
//! dump reader reuses the Phase 3 `read_lammps_dump_with_types` reader for the
//! cell and coordinates but returns the common result model.

use xtalforge_io::trajectory::{Trajectory, read_lammps_dump_with_types};

use crate::Error;
use crate::result::{CalculationResult, ElectronicStep, Snapshot};

use super::common::{ATM_GPA, BAR_GPA, KCAL_MOL_EV, fields, value_after, warn};

/// Parses a LAMMPS log or dump, detected from the first `ITEM: TIMESTEP`.
///
/// # Errors
///
/// Returns the errors of [`parse_log`] or [`parse_dump`].
pub fn parse(text: &str) -> Result<CalculationResult, Error> {
    if text.contains("ITEM: TIMESTEP") {
        parse_dump(text)
    } else {
        parse_log(text)
    }
}

/// Parses a LAMMPS log into the common result model.
///
/// # Errors
///
/// Returns [`Error::Parse`] when the log has no thermo table with an energy
/// column.
pub fn parse_log(text: &str) -> Result<CalculationResult, Error> {
    let lines: Vec<&str> = text.lines().collect();
    let mut result = CalculationResult::new("lammps");
    let units = units(&lines);
    if units.is_none() {
        warn(
            &mut result,
            "lammps_units_assumed",
            "LAMMPS log has no echoed units command; metal units were assumed",
        );
    }
    let style = units.as_deref().unwrap_or("metal");
    let energy_factor = match style {
        "real" => KCAL_MOL_EV,
        "si" => 6.241_509_074e18,
        "metal" | "electron" | "micro" => 1.0,
        "lj" => 1.0,
        _ => 1.0,
    };
    if style == "lj" {
        warn(
            &mut result,
            "lammps_lj_units",
            "LAMMPS log uses reduced (lj) units; energies are dimensionless",
        );
    }
    let pressure_factor = match style {
        "real" => ATM_GPA,
        "si" => 1e-9,
        "metal" | "electron" | "micro" => BAR_GPA,
        "lj" => 1.0,
        _ => BAR_GPA,
    };
    let mut columns: Option<Vec<&str>> = None;
    let mut rows = 0usize;
    for line in &lines {
        let values = fields(line);
        if values.first() == Some(&"Step")
            && values
                .iter()
                .any(|value| matches!(*value, "TotEng" | "PotEng" | "E_pair" | "E_mol"))
        {
            columns = Some(values);
            continue;
        }
        let Some(header) = columns.as_ref() else {
            continue;
        };
        if values.is_empty() {
            continue;
        }
        let Some(step) = values.first().and_then(|value| value.parse::<i64>().ok()) else {
            columns = None;
            continue;
        };
        let index = |name: &str| header.iter().position(|column| *column == name);
        let energy = ["TotEng", "PotEng", "E_pair"]
            .iter()
            .find_map(|name| index(name).and_then(|column| values.get(column)))
            .and_then(|value| value.parse::<f64>().ok());
        let Some(energy) = energy else {
            continue;
        };
        result.energies.push(energy * energy_factor);
        let temperature = index("Temp")
            .and_then(|column| values.get(column))
            .and_then(|value| value.parse::<f64>().ok());
        result.steps.push(ElectronicStep {
            energy: energy * energy_factor,
            fermi: None,
            max_force: None,
        });
        if let Some(temperature) = temperature {
            result.scalars.insert("temperature".to_owned(), temperature);
        }
        if let Some(pressure) = index("Press")
            .and_then(|column| values.get(column))
            .and_then(|value| value.parse::<f64>().ok())
        {
            result
                .scalars
                .insert("pressure_gpa".to_owned(), pressure * pressure_factor);
        }
        result.scalars.insert("last_step".to_owned(), step as f64);
        rows += 1;
    }
    if rows == 0 {
        return Err(Error::Parse(
            "LAMMPS log has no thermo table with an energy column".to_owned(),
        ));
    }
    if text.contains("Minimization stats:") {
        result.convergence.ionic = true;
        result.convergence.ionic_reported = true;
    }
    for line in &lines {
        let trimmed = line.trim();
        let stage = if trimmed.contains("Loop time of") {
            Some("loop")
        } else if trimmed.starts_with("Total wall time") {
            Some("wall")
        } else {
            None
        };
        if let Some(stage) = stage
            && let Some(seconds) = value_after(trimmed, "of").or_else(|| value_after(trimmed, ":"))
        {
            result.push_timing(stage, seconds);
        }
        if trimmed.starts_with("ERROR") {
            result.push_warning("lammps_error", trimmed);
        } else if trimmed.starts_with("WARNING") {
            result.push_warning("lammps_warning", trimmed);
        }
    }
    warn(
        &mut result,
        "lammps_log_has_no_structure",
        "LAMMPS logs do not carry coordinates; read a dump for the final structure",
    );
    result.scalars.insert("thermo_rows".to_owned(), rows as f64);
    Ok(result)
}

/// Parses a LAMMPS text dump with an `element` column.
///
/// # Errors
///
/// Returns [`Error::Parse`] for a truncated dump or a missing element column.
pub fn parse_dump(text: &str) -> Result<CalculationResult, Error> {
    parse_dump_with_types(text, &[])
}

/// Parses a LAMMPS text dump, mapping numeric types through `types`.
///
/// # Errors
///
/// Returns [`Error::Parse`] for a truncated dump, an unknown type or an
/// inconsistent element sequence.
pub fn parse_dump_with_types(text: &str, types: &[&str]) -> Result<CalculationResult, Error> {
    let trajectory = read_lammps_dump_with_types(text, types)
        .map_err(|error| Error::Parse(format!("LAMMPS dump: {error}")))?;
    let mut result = CalculationResult::new("lammps");
    snapshots(&trajectory, &mut result)?;
    result.close_trajectory();
    warn(
        &mut result,
        "lammps_dump_has_no_thermo",
        "LAMMPS dumps carry coordinates only; read the log for energies and timings",
    );
    Ok(result)
}

fn snapshots(trajectory: &Trajectory, result: &mut CalculationResult) -> Result<(), Error> {
    for (index, frame) in trajectory.frames.iter().enumerate() {
        let structure = trajectory
            .structure(index)
            .map_err(|error| Error::Parse(format!("LAMMPS dump frame {index}: {error}")))?;
        result.trajectory.push(Snapshot {
            step: frame.step.unwrap_or(index as u64),
            structure,
            energy: None,
            forces: None,
            stress: None,
        });
    }
    Ok(())
}

fn units(lines: &[&str]) -> Option<String> {
    lines.iter().find_map(|line| {
        let values = fields(line);
        if values.first() == Some(&"units") {
            values.get(1).map(|value| (*value).to_owned())
        } else {
            None
        }
    })
}
