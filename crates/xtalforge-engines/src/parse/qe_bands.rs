//! Quantum ESPRESSO `bands.dat` / `*.bands` parser (Task 18.2.3).
//!
//! The file starts with `&plot nbnd=..., nks=... /` and then alternates one
//! k-point line with one line of band energies in eV.

use crate::Error;
use crate::result::{BandData, CalculationResult};

use super::common::{numbers, warn};

/// Parses a QE bands file into the common result model.
///
/// # Errors
///
/// Returns [`Error::Parse`] for a missing `&plot` header or a truncated band
/// list.
pub fn parse(text: &str) -> Result<CalculationResult, Error> {
    let lines: Vec<&str> = text.lines().collect();
    let header = lines
        .iter()
        .position(|line| line.trim_start().starts_with("&plot"))
        .ok_or_else(|| Error::Parse("QE bands file has no &plot header".to_owned()))?;
    let line = lines.get(header).copied().unwrap_or("");
    let nbnd = header_count(line, "nbnd")?;
    let nks = header_count(line, "nks")?;
    if nbnd == 0 || nks == 0 {
        return Err(Error::Parse(
            "QE bands header declares no bands or no k-points".to_owned(),
        ));
    }
    let mut kpoints = Vec::new();
    let mut energies: Vec<Vec<f64>> = Vec::new();
    let mut index = header + 1;
    while index < lines.len() {
        let line = lines.get(index).copied().unwrap_or("").trim();
        if line.is_empty() {
            index += 1;
            continue;
        }
        if line.starts_with('&') || line == "/" {
            break;
        }
        let point = numbers(line);
        let (Some(x), Some(y), Some(z)) = (point.first(), point.get(1), point.get(2)) else {
            return Err(Error::Parse(format!(
                "QE bands k-point line {} is incomplete",
                kpoints.len() + 1
            )));
        };
        kpoints.push([*x, *y, *z]);
        index += 1;
        let mut row = Vec::new();
        while row.len() < nbnd && index < lines.len() {
            let next = lines.get(index).copied().unwrap_or("").trim();
            if next.starts_with('&') {
                break;
            }
            let values = numbers(next);
            if values.is_empty() {
                index += 1;
                continue;
            }
            if values.len() > nbnd - row.len() {
                return Err(Error::Parse(format!(
                    "QE bands file has too many values for k-point {}",
                    kpoints.len()
                )));
            }
            row.extend(values);
            index += 1;
        }
        if row.len() < nbnd {
            return Err(Error::Parse(format!(
                "QE bands file ends inside k-point {} ({nbnd} bands expected)",
                kpoints.len()
            )));
        }
        row.truncate(nbnd);
        energies.push(row);
        if energies.len() > nks {
            return Err(Error::Parse(format!(
                "QE bands file has more than the declared {nks} k-points"
            )));
        }
    }
    if energies.len() != nks {
        return Err(Error::Parse(format!(
            "QE bands header declares {nks} k-points but {} were read",
            energies.len()
        )));
    }
    let mut result = CalculationResult::new("quantum_espresso");
    result.bands = Some(BandData {
        kpoints,
        weights: Vec::new(),
        energies,
        occupations: Vec::new(),
        fermi: None,
    });
    result.scalars.insert("nbnd".to_owned(), nbnd as f64);
    result.scalars.insert("nks".to_owned(), nks as f64);
    warn(
        &mut result,
        "qe_bands_no_fermi",
        "QE bands files do not carry a Fermi level; supply it from the SCF output",
    );
    Ok(result)
}

fn header_count(line: &str, key: &str) -> Result<usize, Error> {
    let start = line
        .find(key)
        .ok_or_else(|| Error::Parse(format!("QE bands header has no {key}=")))?;
    let rest = line
        .get(start + key.len()..)
        .ok_or_else(|| Error::Parse(format!("QE bands header has no {key}=")))?;
    let rest = rest.trim_start().strip_prefix('=').unwrap_or(rest);
    let digits: String = rest
        .chars()
        .skip_while(|character| !character.is_ascii_digit())
        .take_while(|character| character.is_ascii_digit())
        .collect();
    digits
        .parse()
        .map_err(|_| Error::Parse(format!("QE bands header has a malformed {key} value")))
}
