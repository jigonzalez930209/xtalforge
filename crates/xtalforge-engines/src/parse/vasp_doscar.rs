//! VASP `DOSCAR` parser (Task 18.2.2).
//!
//! Reads the total DOS block and its integrated values. Partial (site and
//! orbital projected) blocks are counted but not expanded into channels; the
//! parser records a warning so callers know they were present.

use crate::Error;
use crate::result::{CalculationResult, DosChannel, DosData};

use super::common::{number, warn};

/// Parses a `DOSCAR` into the common result model.
///
/// # Errors
///
/// Returns [`Error::Parse`] when the DOS header or the declared number of DOS
/// rows is missing or truncated.
pub fn parse(text: &str) -> Result<CalculationResult, Error> {
    let lines: Vec<&str> = text.lines().collect();
    let first = lines
        .first()
        .ok_or_else(|| Error::Parse("DOSCAR is empty".to_owned()))?;
    let header: Vec<usize> = first
        .split_whitespace()
        .filter_map(|value| value.parse().ok())
        .collect();
    let natoms = header
        .first()
        .copied()
        .ok_or_else(|| Error::Parse("DOSCAR header has no atom count".to_owned()))?;
    let nspin = header.get(1).copied().unwrap_or(1);
    let (header_index, nedos, fermi) = dos_header(&lines)
        .ok_or_else(|| Error::Parse("DOSCAR has no Emax/Emin/NEDOS header".to_owned()))?;
    let mut result = CalculationResult::new("vasp");
    let mut energies = Vec::new();
    let mut up = Vec::new();
    let mut down = Vec::new();
    let mut integrated_up = Vec::new();
    let mut integrated_down = Vec::new();
    let mut consumed = 0usize;
    for line in lines.iter().skip(header_index + 1) {
        if consumed >= nedos {
            break;
        }
        let values = line
            .split_whitespace()
            .filter_map(number)
            .collect::<Vec<f64>>();
        if values.is_empty() {
            continue;
        }
        let (Some(energy), Some(dos)) = (values.first(), values.get(1)) else {
            continue;
        };
        energies.push(*energy);
        if values.len() >= 5 {
            let Some(spin_down) = values.get(2) else {
                break;
            };
            let Some(integral_up) = values.get(3) else {
                break;
            };
            let Some(integral_down) = values.get(4) else {
                break;
            };
            up.push(*dos);
            down.push(*spin_down);
            integrated_up.push(*integral_up);
            integrated_down.push(*integral_down);
        } else {
            up.push(*dos);
            integrated_up.push(values.get(2).copied().unwrap_or(0.0));
        }
        consumed += 1;
    }
    if consumed < nedos {
        return Err(Error::Parse(format!(
            "DOSCAR declares {nedos} DOS rows but only {consumed} are present"
        )));
    }
    let spin = !down.is_empty();
    let total = if spin {
        up.iter()
            .zip(down.iter())
            .map(|(first, second)| first + second)
            .collect::<Vec<f64>>()
    } else {
        up.clone()
    };
    let integrated = if spin {
        integrated_up
            .iter()
            .zip(integrated_down.iter())
            .map(|(first, second)| first + second)
            .collect::<Vec<f64>>()
    } else {
        integrated_up
    };
    let mut channels = Vec::new();
    if spin {
        channels.push(DosChannel {
            label: "spin_up".to_owned(),
            values: up,
        });
        channels.push(DosChannel {
            label: "spin_down".to_owned(),
            values: down,
        });
    }
    result.dos = Some(DosData {
        energies,
        total,
        integrated,
        channels,
        fermi: Some(fermi),
    });
    result.fermi = Some(fermi);
    result.scalars.insert("natoms".to_owned(), natoms as f64);
    result.scalars.insert("nspin".to_owned(), nspin as f64);
    result.scalars.insert("nedos".to_owned(), nedos as f64);
    if lines.len() > header_index + 1 + nedos {
        warn(
            &mut result,
            "doscar_partial_states_ignored",
            "DOSCAR contains partial DOS blocks; only the total DOS was read",
        );
    }
    warn(
        &mut result,
        "doscar_has_no_convergence",
        "DOSCAR does not carry convergence flags or total energies",
    );
    Ok(result)
}

fn dos_header(lines: &[&str]) -> Option<(usize, usize, f64)> {
    for (index, line) in lines.iter().enumerate().take(16) {
        let values = line
            .split_whitespace()
            .filter_map(number)
            .collect::<Vec<f64>>();
        if values.len() != 5 {
            continue;
        }
        let Some(count) = values.get(2) else {
            continue;
        };
        if count.fract() != 0.0 || !(1.0..=100_000.0).contains(count) {
            continue;
        }
        let Some(fermi) = values.get(3) else {
            continue;
        };
        let nedos = *count as usize;
        return Some((index, nedos, *fermi));
    }
    None
}
