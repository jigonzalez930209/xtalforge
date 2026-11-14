//! VASP `OSZICAR` parser (Task 18.2.2).
//!
//! `OSZICAR` is the compact energy log: one line per SCF iteration and one
//! `F=` line per ionic step. The parser keeps both series and detects the
//! structural convergence marker.

use crate::Error;
use crate::result::{CalculationResult, ElectronicStep};

use super::common::{number, value_after, warn};

/// Parses an `OSZICAR` into the common result model.
///
/// # Errors
///
/// Returns [`Error::Parse`] when the file has no ionic `F=` energy line, which
/// is the signature of a truncated log.
pub fn parse(text: &str) -> Result<CalculationResult, Error> {
    let mut result = CalculationResult::new("vasp");
    let mut electronic_steps = 0usize;
    let mut ionic_steps = 0usize;
    let mut reached = false;
    for line in text.lines() {
        let trimmed = line.trim();
        if trimmed.is_empty() {
            continue;
        }
        if trimmed.contains("reached required accuracy") {
            reached = true;
            continue;
        }
        if trimmed.to_ascii_lowercase().contains("magnetization") {
            if let Some(value) = value_after(trimmed, "magnetization") {
                result
                    .scalars
                    .insert("total_magnetization".to_owned(), value);
            }
            continue;
        }
        let fields: Vec<&str> = trimmed.split_whitespace().collect();
        let label = fields.first().copied().unwrap_or("");
        if is_electronic_label(label) {
            if let Some(energy) = fields.get(2).and_then(|value| number(value)) {
                result.steps.push(ElectronicStep {
                    energy,
                    fermi: None,
                    max_force: None,
                });
                electronic_steps += 1;
            }
            continue;
        }
        if let Some(energy) = value_after(trimmed, "F=") {
            result.energies.push(energy);
            ionic_steps += 1;
        }
    }
    if ionic_steps == 0 {
        return Err(Error::Parse(
            "OSZICAR has no ionic step (F=) energy line".to_owned(),
        ));
    }
    result.convergence.ionic = reached;
    result.convergence.ionic_reported = reached;
    if reached {
        result.convergence.scf = true;
        result.convergence.scf_reported = false;
    } else {
        warn(
            &mut result,
            "vasp_convergence_missing",
            "OSZICAR has no 'reached required accuracy' marker",
        );
    }
    result
        .scalars
        .insert("electronic_steps".to_owned(), electronic_steps as f64);
    result
        .scalars
        .insert("ionic_steps".to_owned(), ionic_steps as f64);
    Ok(result)
}

fn is_electronic_label(label: &str) -> bool {
    matches!(
        label,
        "DAV:" | "RMM:" | "CG:" | "SDA:" | "DIA:" | "EDDAV:" | "EDRMM:"
    )
}
