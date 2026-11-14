//! Quantum ESPRESSO `projwfc.x` output parser (Task 18.2.3).
//!
//! Reads the `bands (ev)` blocks (k-point, band energies and occupations) and
//! the Löwdin charges printed after `Lowdin Charges:`. The Löwdin charges are
//! mapped onto `CalculationResult::charges` in the atom order printed by QE.

use crate::Error;
use crate::result::{BandData, CalculationResult};

use super::common::{numbers, value_after, warn};

/// Parses a `projwfc.x` output into the common result model.
///
/// # Errors
///
/// Returns [`Error::Parse`] when the output has neither band blocks nor Löwdin
/// charges.
pub fn parse(text: &str) -> Result<CalculationResult, Error> {
    let lines: Vec<&str> = text.lines().collect();
    let mut result = CalculationResult::new("quantum_espresso");
    let mut kpoints = Vec::new();
    let mut energies: Vec<Vec<f64>> = Vec::new();
    let mut occupations: Vec<Vec<f64>> = Vec::new();
    let mut lowdin: Vec<(usize, f64)> = Vec::new();
    let mut index = 0usize;
    while index < lines.len() {
        let line = lines.get(index).copied().unwrap_or("").trim();
        if line.starts_with("k =") && line.contains("bands") {
            let point = numbers(line.split_once('=').map_or("", |(_, rest)| rest));
            let Some(values) = point.get(..3) else {
                index += 1;
                continue;
            };
            kpoints.push([
                values.first().copied().unwrap_or(0.0),
                values.get(1).copied().unwrap_or(0.0),
                values.get(2).copied().unwrap_or(0.0),
            ]);
            energies.push(Vec::new());
            occupations.push(Vec::new());
            index += 1;
            while index < lines.len() {
                let row = lines.get(index).copied().unwrap_or("");
                let trimmed = row.trim();
                if trimmed.is_empty() {
                    index += 1;
                    continue;
                }
                if trimmed.starts_with("k =")
                    || trimmed.starts_with("Lowdin")
                    || trimmed.starts_with("Program")
                    || trimmed.starts_with("psi")
                    || trimmed.starts_with("sum")
                {
                    break;
                }
                let Some(first) = trimmed.split_whitespace().next() else {
                    break;
                };
                let indexed = first.parse::<usize>().is_ok();
                let row_values = numbers(trimmed);
                let (energy, occupation) = if indexed {
                    if row_values.len() < 3 {
                        break;
                    }
                    (row_values.get(1).copied(), row_values.get(2).copied())
                } else {
                    (row_values.first().copied(), row_values.get(1).copied())
                };
                let (Some(energy), Some(occupation)) = (energy, occupation) else {
                    break;
                };
                if let (Some(row_energy), Some(row_occupation)) =
                    (energies.last_mut(), occupations.last_mut())
                {
                    row_energy.push(energy);
                    row_occupation.push(occupation);
                }
                index += 1;
            }
            continue;
        }
        if line.starts_with("Lowdin Charges") {
            index += 1;
            while index < lines.len() {
                let row = lines.get(index).copied().unwrap_or("").trim();
                if row.is_empty() {
                    index += 1;
                    continue;
                }
                if !row.starts_with("Atom") {
                    break;
                }
                let site = row
                    .split_once(':')
                    .and_then(|(left, _)| left.split_whitespace().last())
                    .and_then(|value| value.parse::<usize>().ok());
                let charge = value_after(row, "total charge");
                if let (Some(site), Some(charge)) = (site, charge) {
                    lowdin.push((site, charge));
                }
                index += 1;
            }
            continue;
        }
        index += 1;
    }
    if energies.is_empty() && lowdin.is_empty() {
        return Err(Error::Parse(
            "projwfc output has no bands or Löwdin charges".to_owned(),
        ));
    }
    if !kpoints.is_empty() {
        let band_count = energies.iter().map(Vec::len).max().unwrap_or(0);
        let kpoint_count = kpoints.len();
        result.bands = Some(BandData {
            kpoints,
            weights: Vec::new(),
            energies,
            occupations,
            fermi: None,
        });
        result.scalars.insert("nbnd".to_owned(), band_count as f64);
        result.scalars.insert("nks".to_owned(), kpoint_count as f64);
        warn(
            &mut result,
            "qe_projwfc_no_fermi",
            "projwfc output does not carry a Fermi level",
        );
    }
    if !lowdin.is_empty() {
        lowdin.sort_by_key(|(site, _)| *site);
        result.charges = lowdin.into_iter().map(|(_, charge)| charge).collect();
    }
    Ok(result)
}
