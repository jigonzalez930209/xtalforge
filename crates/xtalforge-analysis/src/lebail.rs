//! Le Bail intensity extraction (Task 7.6.1).

#![allow(clippy::indexing_slicing)] // reflection indices are bounded by the loops

use serde::{Deserialize, Serialize};
use xtalforge_core::Error;
use xtalforge_io::xy::PatternPoint;

use crate::profiles::{Profile, ProfileKind};
use crate::pxrd::{PxrdOptions, Radiation, calculate_pxrd};
use xtalforge_core::Structure;

/// Options for the Le Bail iteration.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct LebailOptions {
    /// Profile FWHM in degrees.
    pub fwhm: f64,
    /// Integration window half-width in FWHM units.
    pub window: f64,
    /// Number of iterations.
    pub iterations: usize,
    /// Profile shape.
    pub kind: ProfileKind,
    /// Pseudo-Voigt mixing.
    pub eta: f64,
}

impl Default for LebailOptions {
    fn default() -> Self {
        Self {
            fwhm: 0.1,
            window: 3.0,
            iterations: 12,
            kind: ProfileKind::PseudoVoigt,
            eta: 0.5,
        }
    }
}

/// Result of the Le Bail extraction.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct LebailResult {
    /// Extracted intensity per reflection, in reflection order.
    pub intensities: Vec<f64>,
    /// Reflection positions in degrees.
    pub positions: Vec<f64>,
    /// Calculated pattern at the observed points.
    pub pattern: Vec<PatternPoint>,
    /// Profile R factor on the fitted points.
    pub r_factor: f64,
}

/// Extracts reflection intensities from a measured pattern.
///
/// # Errors
///
/// Returns [`Error::InvalidInput`] for no reflections, too few points or bad
/// options.
pub fn le_bail(
    structure: &Structure,
    observed: &[PatternPoint],
    options: &LebailOptions,
) -> Result<LebailResult, Error> {
    if observed.len() < 3 {
        return Err(Error::InvalidInput(
            "Le Bail needs at least three observed points".to_owned(),
        ));
    }
    if !options.fwhm.is_finite() || options.fwhm <= 0.0 || options.iterations == 0 {
        return Err(Error::InvalidInput(
            "Le Bail options must be positive".to_owned(),
        ));
    }
    let pxrd = calculate_pxrd(
        structure,
        &PxrdOptions {
            radiation: Radiation::CuKa,
            two_theta_range: [
                observed.first().map_or(5.0, |point| point.two_theta),
                observed.last().map_or(80.0, |point| point.two_theta),
            ],
            step: 0.02,
            fwhm: options.fwhm,
        },
    )?;
    if pxrd.reflections.is_empty() {
        return Err(Error::InvalidInput(
            "the structure has no reflections in the observed range".to_owned(),
        ));
    }
    let positions: Vec<f64> = pxrd
        .reflections
        .iter()
        .map(|reflection| reflection.two_theta)
        .collect();
    let mut intensities: Vec<f64> = pxrd
        .reflections
        .iter()
        .map(|reflection| reflection.intensity.max(1e-9))
        .collect();
    for _ in 0..options.iterations {
        let calculated = profile_sum(&positions, &intensities, observed, options);
        let mut next = vec![0.0; positions.len()];
        for (index, position) in positions.iter().enumerate() {
            let half_width = options.window * options.fwhm;
            let mut observed_sum = 0.0;
            let mut calculated_sum = 0.0;
            for (point, value) in observed.iter().zip(calculated.iter()) {
                if (point.two_theta - position).abs() <= half_width {
                    observed_sum += point.intensity;
                    calculated_sum += value;
                }
            }
            if calculated_sum > 1e-12 {
                next[index] = intensities[index] * observed_sum / calculated_sum;
            } else {
                next[index] = intensities[index];
            }
        }
        intensities = next;
    }
    let pattern_values = profile_sum(&positions, &intensities, observed, options);
    let pattern: Vec<PatternPoint> = observed
        .iter()
        .zip(pattern_values.iter())
        .map(|(point, value)| PatternPoint {
            two_theta: point.two_theta,
            intensity: *value,
            uncertainty: None,
        })
        .collect();
    let numerator: f64 = observed
        .iter()
        .zip(pattern.iter())
        .map(|(point, fit)| (point.intensity - fit.intensity).abs())
        .sum();
    let denominator: f64 = observed.iter().map(|point| point.intensity.abs()).sum();
    Ok(LebailResult {
        r_factor: if denominator > 0.0 {
            numerator / denominator
        } else {
            0.0
        },
        intensities,
        positions,
        pattern,
    })
}

fn profile_sum(
    positions: &[f64],
    intensities: &[f64],
    observed: &[PatternPoint],
    options: &LebailOptions,
) -> Vec<f64> {
    let mut values = vec![0.0; observed.len()];
    for (position, intensity) in positions.iter().zip(intensities.iter()) {
        let profile = Profile {
            kind: options.kind,
            center: *position,
            fwhm: options.fwhm,
            eta: options.eta,
            height: *intensity,
        };
        for (value, point) in values.iter_mut().zip(observed.iter()) {
            *value += profile.value(point.two_theta);
        }
    }
    values
}

/// Pawley intensities from a linear least-squares fit.
///
/// # Errors
///
/// Returns [`Error::InvalidInput`] for too few points or a singular normal
/// matrix.
pub fn pawley(
    structure: &Structure,
    observed: &[PatternPoint],
    options: &LebailOptions,
) -> Result<LebailResult, Error> {
    if observed.len() < 3 {
        return Err(Error::InvalidInput(
            "Pawley needs at least three observed points".to_owned(),
        ));
    }
    let pxrd = calculate_pxrd(
        structure,
        &PxrdOptions {
            radiation: Radiation::CuKa,
            two_theta_range: [
                observed.first().map_or(5.0, |point| point.two_theta),
                observed.last().map_or(80.0, |point| point.two_theta),
            ],
            step: 0.02,
            fwhm: options.fwhm,
        },
    )?;
    let first = observed.first().map_or(0.0, |point| point.two_theta);
    let last = observed.last().map_or(0.0, |point| point.two_theta);
    let margin = options.window * options.fwhm;
    let positions: Vec<f64> = pxrd
        .reflections
        .iter()
        .filter(|reflection| {
            reflection.intensity > 1e-8
                && reflection.two_theta >= first - margin
                && reflection.two_theta <= last + margin
        })
        .map(|reflection| reflection.two_theta)
        .collect();
    if positions.is_empty() {
        return Err(Error::InvalidInput(
            "the structure has no reflections in the observed range".to_owned(),
        ));
    }
    let count = positions.len();
    let mut normal = vec![vec![0.0; count]; count];
    let mut rhs = vec![0.0; count];
    for point in observed {
        let row: Vec<f64> = positions
            .iter()
            .map(|position| {
                Profile {
                    kind: options.kind,
                    center: *position,
                    fwhm: options.fwhm,
                    eta: options.eta,
                    height: 1.0,
                }
                .value(point.two_theta)
            })
            .collect();
        for first in 0..count {
            rhs[first] += row[first] * point.intensity;
            for second in 0..count {
                normal[first][second] += row[first] * row[second];
            }
        }
    }
    let scale = normal
        .iter()
        .flat_map(|row| row.iter())
        .fold(0.0_f64, |value, entry| value.max(entry.abs()))
        .max(1.0);
    for (index, row) in normal.iter_mut().enumerate() {
        if let Some(value) = row.get_mut(index) {
            *value += 1e-10 * scale;
        }
    }
    let Some(intensities) = solve(&normal, &rhs) else {
        return Err(Error::InvalidInput(
            "the Pawley normal matrix is singular".to_owned(),
        ));
    };
    let pattern_values = profile_sum(&positions, &intensities, observed, options);
    let pattern: Vec<PatternPoint> = observed
        .iter()
        .zip(pattern_values.iter())
        .map(|(point, value)| PatternPoint {
            two_theta: point.two_theta,
            intensity: *value,
            uncertainty: None,
        })
        .collect();
    let numerator: f64 = observed
        .iter()
        .zip(pattern.iter())
        .map(|(point, fit)| (point.intensity - fit.intensity).abs())
        .sum();
    let denominator: f64 = observed.iter().map(|point| point.intensity.abs()).sum();
    Ok(LebailResult {
        intensities,
        positions,
        pattern,
        r_factor: if denominator > 0.0 {
            numerator / denominator
        } else {
            0.0
        },
    })
}

fn solve(matrix: &[Vec<f64>], rhs: &[f64]) -> Option<Vec<f64>> {
    let n = matrix.len();
    let mut augmented: Vec<Vec<f64>> = matrix
        .iter()
        .zip(rhs.iter())
        .map(|(row, value)| {
            let mut line = row.clone();
            line.push(*value);
            line
        })
        .collect();
    for column in 0..n {
        let pivot = (column..n).max_by(|a, b| {
            augmented[*a][column]
                .abs()
                .partial_cmp(&augmented[*b][column].abs())
                .unwrap_or(std::cmp::Ordering::Equal)
        })?;
        if augmented[pivot][column].abs() < 1e-14 {
            return None;
        }
        augmented.swap(column, pivot);
        let divisor = augmented[column][column];
        for value in augmented[column].iter_mut() {
            *value /= divisor;
        }
        let pivot_row = augmented[column].clone();
        for (row, line) in augmented.iter_mut().enumerate() {
            if row == column {
                continue;
            }
            let factor = line[column];
            for (slot, pivot_value) in line.iter_mut().zip(pivot_row.iter()) {
                *slot -= factor * pivot_value;
            }
        }
    }
    Some(augmented.iter().map(|row| row[n]).collect())
}
