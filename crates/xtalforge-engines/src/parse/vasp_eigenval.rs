//! VASP `EIGENVAL` parser (Task 18.2.2).
//!
//! The file lists every k-point followed by one row per band with its energy
//! in eV and occupation. K-points without band rows (the extra header line in
//! the format) are skipped, and truncated files are rejected.

use crate::Error;
use crate::result::{BandData, CalculationResult};

use super::common::{number, warn};

/// Parses an `EIGENVAL` into the common result model.
///
/// # Errors
///
/// Returns [`Error::Parse`] for a truncated header, missing `CAR` marker or a
/// k-point with an incomplete band list.
pub fn parse(text: &str) -> Result<CalculationResult, Error> {
    let mut lines = text.lines().filter(|line| !line.trim().is_empty());
    let header = lines
        .next()
        .ok_or_else(|| Error::Parse("EIGENVAL has no header".to_owned()))?;
    let counts: Vec<usize> = header
        .split_whitespace()
        .filter_map(|value| value.parse().ok())
        .collect();
    let natoms = counts
        .first()
        .copied()
        .ok_or_else(|| Error::Parse("EIGENVAL header has no atom count".to_owned()))?;
    let nbands = counts
        .get(1)
        .copied()
        .ok_or_else(|| Error::Parse("EIGENVAL header has no band count".to_owned()))?;
    let mut found_car = false;
    for line in lines.by_ref() {
        if line.trim() == "CAR" {
            found_car = true;
            break;
        }
    }
    if !found_car {
        return Err(Error::Parse(
            "EIGENVAL is missing the CAR species marker".to_owned(),
        ));
    }
    let mut result = CalculationResult::new("vasp");
    let species: Vec<String> = lines
        .next()
        .ok_or_else(|| Error::Parse("EIGENVAL ends before the species line".to_owned()))?
        .split_whitespace()
        .map(str::to_owned)
        .collect();
    let species_counts: Vec<usize> = lines
        .next()
        .ok_or_else(|| Error::Parse("EIGENVAL ends before the species counts".to_owned()))?
        .split_whitespace()
        .filter_map(|value| value.parse().ok())
        .collect();
    if !species_counts.is_empty() && species_counts.iter().sum::<usize>() != natoms {
        warn(
            &mut result,
            "vasp_atom_count_mismatch",
            "EIGENVAL species counts do not add up to the header atom count",
        );
    }
    let mut series = Series::default();
    for line in lines {
        let fields: Vec<&str> = line.split_whitespace().collect();
        if fields.len() == 4 {
            let point = (
                fields.first().and_then(|value| number(value)),
                fields.get(1).and_then(|value| number(value)),
                fields.get(2).and_then(|value| number(value)),
                fields.get(3).and_then(|value| number(value)),
            );
            if let (Some(x), Some(y), Some(z), Some(weight)) = point {
                let inside = [x, y, z].iter().all(|value| value.abs() <= 1.0 + 1e-6);
                if inside && weight > 0.0 {
                    series.finish();
                    series.current = Some(([x, y, z], weight));
                    continue;
                }
            }
        }
        if fields.len() == 3 {
            let band = fields.first().and_then(|value| value.parse::<usize>().ok());
            let energy = fields.get(1).and_then(|value| number(value));
            let occupation = fields.get(2).and_then(|value| number(value));
            if band.is_some()
                && let (Some(energy), Some(occupation)) = (energy, occupation)
            {
                series.energies_current.push(energy);
                series.occupations_current.push(occupation);
            }
        }
    }
    series.finish();
    if series.energies.is_empty() {
        return Err(Error::Parse("EIGENVAL has no k-point band rows".to_owned()));
    }
    if let Some((index, row)) = series
        .energies
        .iter()
        .enumerate()
        .find(|(_, row)| row.len() != nbands)
    {
        return Err(Error::Parse(format!(
            "EIGENVAL k-point {} has {} bands but the header declares {nbands}",
            index + 1,
            row.len()
        )));
    }
    result.bands = Some(BandData {
        kpoints: series.kpoints,
        weights: series.weights,
        energies: series.energies,
        occupations: series.occupations,
        fermi: None,
    });
    result.scalars.insert("natoms".to_owned(), natoms as f64);
    result.scalars.insert("nbands".to_owned(), nbands as f64);
    result
        .scalars
        .insert("species".to_owned(), species.len() as f64);
    warn(
        &mut result,
        "eigenval_has_no_convergence",
        "EIGENVAL does not carry convergence flags or total energies",
    );
    Ok(result)
}

#[derive(Default)]
struct Series {
    kpoints: Vec<[f64; 3]>,
    weights: Vec<f64>,
    energies: Vec<Vec<f64>>,
    occupations: Vec<Vec<f64>>,
    current: Option<([f64; 3], f64)>,
    energies_current: Vec<f64>,
    occupations_current: Vec<f64>,
}

impl Series {
    fn finish(&mut self) {
        if let Some((point, weight)) = self.current.take()
            && !self.energies_current.is_empty()
        {
            self.kpoints.push(point);
            self.weights.push(weight);
            self.energies
                .push(std::mem::take(&mut self.energies_current));
            self.occupations
                .push(std::mem::take(&mut self.occupations_current));
        }
    }
}
