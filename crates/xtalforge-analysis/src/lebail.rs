//! Le Bail intensity extraction (Task 7.6.1).

#![allow(clippy::indexing_slicing)] // reflection indices are bounded by the loops

use serde::{Deserialize, Serialize};
use xtalforge_core::Error;
use xtalforge_io::xy::PatternPoint;

use crate::profiles::{Profile, ProfileKind};
use crate::pxrd::{PxrdOptions, Radiation, calculate_pxrd};
use xtalforge_core::Structure;

mod esd;

use self::esd::{intensity_esds, solve, statistics, weights};

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
    /// Weighted profile R factor.
    pub rwp: f64,
    /// Reduced chi-squared of the extracted intensities.
    pub chi_squared: f64,
    /// Weighted residual sum of squares of the profile fit.
    pub weighted_sum: f64,
    /// Standard uncertainty of every intensity, in reflection order.
    pub intensity_esds: Vec<f64>,
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
    let (background, net) = split_background(observed, options);
    for _ in 0..options.iterations {
        let calculated = profile_sum(&positions, &intensities, &net, options);
        let mut next = vec![0.0; positions.len()];
        // Le Bail partition: every observed point is shared among the
        // reflections in proportion to their calculated contribution,
        // H_k <- H_k * sum_i P_ik (y_obs,i / y_calc,i) / sum_i P_ik.
        for (index, position) in positions.iter().enumerate() {
            let half_width = options.window * options.fwhm;
            let shape = Profile {
                kind: options.kind,
                center: *position,
                fwhm: options.fwhm,
                eta: options.eta,
                height: 1.0,
            };
            let mut numerator = 0.0;
            let mut denominator = 0.0;
            for (point, value) in net.iter().zip(calculated.iter()) {
                if (point.two_theta - position).abs() > half_width || *value <= 1e-12 {
                    continue;
                }
                let weight = shape.value(point.two_theta);
                numerator += weight * point.intensity.max(0.0) / value;
                denominator += weight;
            }
            next[index] = if denominator > 1e-12 {
                intensities[index] * numerator / denominator
            } else {
                intensities[index]
            };
        }
        intensities = next;
    }
    let pattern_values: Vec<f64> = profile_sum(&positions, &intensities, observed, options)
        .into_iter()
        .zip(&background)
        .map(|(value, base)| value + base)
        .collect();
    let pattern: Vec<PatternPoint> = observed
        .iter()
        .zip(pattern_values.iter())
        .map(|(point, value)| PatternPoint {
            two_theta: point.two_theta,
            intensity: *value,
            uncertainty: None,
        })
        .collect();
    let weights = weights(observed);
    let (r_factor, rwp, chi_squared, weighted_sum) =
        statistics(observed, &pattern_values, &weights);
    let intensity_esds = intensity_esds(&positions, observed, options, &weights, chi_squared);
    Ok(LebailResult {
        r_factor,
        rwp,
        chi_squared,
        weighted_sum,
        intensity_esds,
        intensities,
        positions,
        pattern,
    })
}

/// SNIP background under the observed pattern and the net observed points.
///
/// The clipping window spans about ten peak widths, so a Bragg line is
/// clipped while a slowly varying baseline survives.
fn split_background(
    observed: &[PatternPoint],
    options: &LebailOptions,
) -> (Vec<f64>, Vec<PatternPoint>) {
    let step = match (observed.first(), observed.last()) {
        (Some(first), Some(last)) if observed.len() > 1 => {
            ((last.two_theta - first.two_theta) / (observed.len() - 1) as f64).abs()
        }
        _ => 0.0,
    };
    let window = if step > 0.0 {
        ((10.0 * options.fwhm / step).ceil() as usize).max(1)
    } else {
        0
    };
    let background = if window == 0 {
        vec![0.0; observed.len()]
    } else {
        crate::peaks::background_envelope(observed, window)
    };
    let net = observed
        .iter()
        .zip(&background)
        .map(|(point, base)| PatternPoint {
            intensity: point.intensity - base,
            ..*point
        })
        .collect();
    (background, net)
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
    let weights = weights(observed);
    let (background, net) = split_background(observed, options);
    let mut normal = vec![vec![0.0; count]; count];
    let mut rhs = vec![0.0; count];
    for (row_index, point) in net.iter().enumerate() {
        let weight = weights.get(row_index).copied().unwrap_or(1.0);
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
            rhs[first] += weight * row[first] * point.intensity;
            for second in 0..count {
                normal[first][second] += weight * row[first] * row[second];
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
    let pattern_values: Vec<f64> = profile_sum(&positions, &intensities, observed, options)
        .into_iter()
        .zip(&background)
        .map(|(value, base)| value + base)
        .collect();
    let pattern: Vec<PatternPoint> = observed
        .iter()
        .zip(pattern_values.iter())
        .map(|(point, value)| PatternPoint {
            two_theta: point.two_theta,
            intensity: *value,
            uncertainty: None,
        })
        .collect();
    let (r_factor, rwp, chi_squared, weighted_sum) =
        statistics(observed, &pattern_values, &weights);
    let intensity_esds = intensity_esds(&positions, observed, options, &weights, chi_squared);
    Ok(LebailResult {
        intensities,
        positions,
        pattern,
        r_factor,
        rwp,
        chi_squared,
        weighted_sum,
        intensity_esds,
    })
}
