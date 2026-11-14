//! VASP `PROCAR` parser (Task 18.2.2).
//!
//! Reads k-points, band energies and occupations, and the site/orbital
//! projection weights printed for every band. Atom and band indices in the
//! result are zero-based; orbital labels are kept exactly as printed, except
//! the `tot` summary column which duplicates their sum.

use crate::Error;
use crate::result::{BandData, BandProjection, CalculationResult, OrbitalWeight, ProjectionData};

use super::common::{number, value_after, warn};

/// Parses a `PROCAR` into the common result model.
///
/// # Errors
///
/// Returns [`Error::Parse`] when the file has no k-point or band blocks.
pub fn parse(text: &str) -> Result<CalculationResult, Error> {
    let mut result = CalculationResult::new("vasp");
    let mut kpoints = Vec::new();
    let mut weights = Vec::new();
    let mut energies: Vec<Vec<f64>> = Vec::new();
    let mut occupations: Vec<Vec<f64>> = Vec::new();
    let mut current_kpoint: Option<usize> = None;
    let mut current_point = [0.0; 3];
    let mut labels: Vec<String> = Vec::new();
    let mut projections = ProjectionData::default();
    let mut band_count = 0usize;
    let mut declared = [0usize; 3];
    for line in text.lines() {
        let trimmed = line.trim();
        if trimmed.is_empty() {
            continue;
        }
        if trimmed.contains("# of k-points") {
            declared = [
                count_after(trimmed, "# of k-points:"),
                count_after(trimmed, "# of bands:"),
                count_after(trimmed, "# of ions:"),
            ];
            continue;
        }
        if trimmed.starts_with("k-point") {
            let Some((_, rest)) = trimmed.split_once(':') else {
                continue;
            };
            let (point_text, weight_text) = rest.split_once("weight").unwrap_or((rest, ""));
            let point = crate::parse::common::numbers(point_text);
            current_point = [
                point.first().copied().unwrap_or(0.0),
                point.get(1).copied().unwrap_or(0.0),
                point.get(2).copied().unwrap_or(0.0),
            ];
            kpoints.push(current_point);
            weights.push(value_after(weight_text, "=").unwrap_or(0.0));
            energies.push(Vec::new());
            occupations.push(Vec::new());
            current_kpoint = kpoints.len().checked_sub(1);
            continue;
        }
        if trimmed.starts_with("band ") {
            let Some((_, rest)) = trimmed.split_once("energy") else {
                continue;
            };
            let Some(energy) = crate::parse::common::numbers(rest).first().copied() else {
                continue;
            };
            let occupation = value_after(trimmed, "# occ.").unwrap_or(0.0);
            let Some(index) = current_kpoint else {
                return Err(Error::Parse(
                    "PROCAR has a band block before the first k-point".to_owned(),
                ));
            };
            if let (Some(row), Some(occ)) = (energies.get_mut(index), occupations.get_mut(index)) {
                row.push(energy);
                occ.push(occupation);
            }
            projections.bands.push(BandProjection {
                kpoint: current_point,
                band: band_count,
                energy,
                occupation,
                weights: Vec::new(),
            });
            band_count += 1;
            continue;
        }
        if trimmed.starts_with("ion") {
            let fields: Vec<&str> = trimmed.split_whitespace().collect();
            labels = fields
                .iter()
                .skip(1)
                .filter(|label| !label.eq_ignore_ascii_case("tot"))
                .map(|label| (*label).to_owned())
                .collect();
            continue;
        }
        let fields: Vec<&str> = trimmed.split_whitespace().collect();
        let Some(atom) = fields.first().and_then(|value| value.parse::<usize>().ok()) else {
            continue;
        };
        let Some(projection) = projections.bands.last_mut() else {
            continue;
        };
        for (label, value) in labels.iter().zip(fields.iter().skip(1)) {
            if let Some(weight) = number(value) {
                projection.weights.push(OrbitalWeight {
                    atom: atom.saturating_sub(1),
                    orbital: label.clone(),
                    weight,
                });
            }
        }
    }
    if projections.bands.is_empty() {
        return Err(Error::Parse("PROCAR has no band blocks".to_owned()));
    }
    if energies.is_empty() {
        return Err(Error::Parse("PROCAR has no k-point blocks".to_owned()));
    }
    result.bands = Some(BandData {
        kpoints,
        weights,
        energies,
        occupations,
        fermi: None,
    });
    result.projections = Some(projections);
    result
        .scalars
        .insert("kpoints".to_owned(), declared[0] as f64);
    result
        .scalars
        .insert("bands".to_owned(), declared[1] as f64);
    result.scalars.insert("ions".to_owned(), declared[2] as f64);
    warn(
        &mut result,
        "procar_has_no_convergence",
        "PROCAR does not carry convergence flags or total energies",
    );
    Ok(result)
}

fn count_after(line: &str, marker: &str) -> usize {
    value_after(line, marker).map_or(0, |value| value as usize)
}
