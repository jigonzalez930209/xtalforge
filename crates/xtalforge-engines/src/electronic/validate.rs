//! Shape and value validation shared by the electronic-structure types.

use std::collections::BTreeSet;

use super::bands::{Band, KPointLabel};
use super::dos::Projection;
use xtalforge_core::Error;

/// Builds an input error with a message.
pub(crate) fn invalid(message: impl Into<String>) -> Error {
    Error::InvalidInput(message.into())
}

/// Checks that every k-point coordinate is finite.
pub(crate) fn kpoints(kpoints: &[[f64; 3]]) -> Result<(), Error> {
    for (index, point) in kpoints.iter().enumerate() {
        if point.iter().any(|value| !value.is_finite()) {
            return Err(invalid(format!("k-point {index} is not finite")));
        }
    }
    Ok(())
}

/// Checks a `[band][kpoint]` matrix for shape and finiteness.
pub(crate) fn energies(energies: &[Vec<f64>], kpoints: usize, field: &str) -> Result<(), Error> {
    if energies.is_empty() {
        return Err(invalid(format!("{field} must contain at least one band")));
    }
    for (band, row) in energies.iter().enumerate() {
        if row.len() != kpoints {
            return Err(invalid(format!(
                "{field} band {band} has {} values but the path has {kpoints} k-points",
                row.len()
            )));
        }
        if row.iter().any(|value| !value.is_finite()) {
            return Err(invalid(format!("{field} band {band} is not finite")));
        }
    }
    Ok(())
}

/// Checks one occupation row against its expected length.
fn occupation_row(row: &[f64], expected: usize, band: usize, field: &str) -> Result<(), Error> {
    if row.len() != expected {
        return Err(invalid(format!(
            "{field} band {band} has {} values but the energies have {expected}",
            row.len()
        )));
    }
    for value in row {
        if !value.is_finite() || *value < 0.0 || *value > 2.0 {
            return Err(invalid(format!(
                "{field} band {band} contains {value}, expected a value in [0, 2]"
            )));
        }
    }
    Ok(())
}

/// Checks an occupation matrix against its energies.
pub(crate) fn occupations(
    occupations: &[Vec<f64>],
    energies: &[Vec<f64>],
    field: &str,
) -> Result<(), Error> {
    if occupations.len() != energies.len() {
        return Err(invalid(format!(
            "{field} has {} bands but the energies have {}",
            occupations.len(),
            energies.len()
        )));
    }
    for (band, (row, reference)) in occupations.iter().zip(energies).enumerate() {
        occupation_row(row, reference.len(), band, field)?;
    }
    Ok(())
}

/// Checks band rows for shape and finiteness.
pub(crate) fn bands(bands: &[Band], kpoints: usize, field: &str) -> Result<(), Error> {
    if bands.is_empty() {
        return Err(invalid(format!("{field} must contain at least one band")));
    }
    for (index, band) in bands.iter().enumerate() {
        if band.energies.len() != kpoints {
            return Err(invalid(format!(
                "{field} band {index} has {} values but the path has {kpoints} k-points",
                band.energies.len()
            )));
        }
        if band.energies.iter().any(|value| !value.is_finite()) {
            return Err(invalid(format!("{field} band {index} is not finite")));
        }
    }
    Ok(())
}

/// Checks an occupation matrix against band rows.
pub(crate) fn occupations_against_bands(
    occupations: &[Vec<f64>],
    bands: &[Band],
    field: &str,
) -> Result<(), Error> {
    if occupations.len() != bands.len() {
        return Err(invalid(format!(
            "{field} has {} bands but the energies have {}",
            occupations.len(),
            bands.len()
        )));
    }
    for (band, (row, reference)) in occupations.iter().zip(bands).enumerate() {
        occupation_row(row, reference.energies.len(), band, field)?;
    }
    Ok(())
}

/// Checks labels for non-empty text and in-range, unique indices.
pub(crate) fn labels(labels: &[KPointLabel], kpoints: usize) -> Result<(), Error> {
    let mut seen = BTreeSet::new();
    for label in labels {
        if label.index >= kpoints {
            return Err(invalid(format!(
                "label {:?} points at k-point {} but the path has {kpoints} k-points",
                label.label, label.index
            )));
        }
        if label.label.trim().is_empty() {
            return Err(invalid("a k-point label must not be empty"));
        }
        if !seen.insert(label.index) {
            return Err(invalid(format!(
                "k-point {} carries more than one label",
                label.index
            )));
        }
    }
    Ok(())
}

/// Checks a path coordinate for shape, finiteness and monotonicity.
pub(crate) fn distances(distances: &[f64], kpoints: usize, field: &str) -> Result<(), Error> {
    if distances.len() != kpoints {
        return Err(invalid(format!(
            "{field} has {} values but the path has {kpoints} k-points",
            distances.len()
        )));
    }
    if distances.iter().any(|value| !value.is_finite()) {
        return Err(invalid(format!("{field} contains a non-finite value")));
    }
    for window in distances.windows(2) {
        if let [left, right] = window
            && right < left
        {
            return Err(invalid(format!("{field} must be non-decreasing")));
        }
    }
    Ok(())
}

/// Checks projection indices and weights.
pub(crate) fn projections(
    projections: &[Projection],
    bands: usize,
    kpoints: usize,
) -> Result<(), Error> {
    let mut states = BTreeSet::new();
    for (index, projection) in projections.iter().enumerate() {
        if !states.insert((projection.band, projection.kpoint)) {
            return Err(invalid(format!(
                "projection {index} repeats band {} at k-point {}",
                projection.band, projection.kpoint
            )));
        }
        if projection.band >= bands {
            return Err(invalid(format!(
                "projection {index} refers to band {} but there are {bands} bands",
                projection.band
            )));
        }
        if projection.kpoint >= kpoints {
            return Err(invalid(format!(
                "projection {index} refers to k-point {} but the path has {kpoints} k-points",
                projection.kpoint
            )));
        }
        for weight in &projection.weights {
            if weight.element.trim().is_empty() {
                return Err(invalid(format!(
                    "projection {index} has a weight with an empty element"
                )));
            }
            if weight.orbital.as_str().trim().is_empty() {
                return Err(invalid(format!(
                    "projection {index} has a weight with an empty orbital"
                )));
            }
            if !weight.weight.is_finite() || weight.weight < 0.0 {
                return Err(invalid(format!(
                    "projection {index} has weight {} for {}/{}, expected a finite non-negative value",
                    weight.weight,
                    weight.element,
                    weight.orbital.as_str()
                )));
            }
        }
    }
    Ok(())
}

/// Checks an optional Fermi energy.
pub(crate) fn fermi(fermi: Option<f64>) -> Result<(), Error> {
    if let Some(value) = fermi
        && !value.is_finite()
    {
        return Err(invalid("the Fermi energy is not finite"));
    }
    Ok(())
}
