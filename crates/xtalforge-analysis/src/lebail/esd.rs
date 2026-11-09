//! Weights, R-factors and intensity covariances for the Le Bail fit (Task
//! 7.6.3).

use xtalforge_io::xy::PatternPoint;

use super::LebailOptions;
use crate::profiles::Profile;

/// Reflection count above which the intensity esds use the covariance
/// diagonal instead of a full inverse.
const MAX_EXACT_ESD: usize = 128;

/// Confidence weight of one observed point.
pub(super) fn weights(observed: &[PatternPoint]) -> Vec<f64> {
    observed
        .iter()
        .map(|point| 1.0 / point.intensity.abs().max(1.0))
        .collect()
}

/// `Rp`, `Rwp` and the reduced chi-squared on the fitted points.
pub(super) fn statistics(
    observed: &[PatternPoint],
    calculated: &[f64],
    weights: &[f64],
) -> (f64, f64, f64, f64) {
    let mut numerator = 0.0;
    let mut denominator = 0.0;
    let mut weighted = 0.0;
    let mut weighted_observed = 0.0;
    for (index, (point, value)) in observed.iter().zip(calculated.iter()).enumerate() {
        let weight = weights.get(index).copied().unwrap_or(1.0);
        numerator += (point.intensity - value).abs();
        denominator += point.intensity.abs();
        weighted += weight * (point.intensity - value).powi(2);
        weighted_observed += weight * point.intensity.powi(2);
    }
    let r_factor = if denominator > 0.0 {
        numerator / denominator
    } else {
        0.0
    };
    let rwp = if weighted_observed > 0.0 {
        (weighted / weighted_observed).sqrt()
    } else {
        0.0
    };
    let degrees = observed
        .len()
        .saturating_sub(calculated.len().max(1))
        .max(1) as f64;
    (r_factor, rwp, weighted / degrees, weighted)
}

/// Covariance diagonal of the intensities at fixed positions and profile.
///
/// The normal matrix is assembled over the observed points and inverted when
/// the reflection list is small enough; above [`MAX_EXACT_ESD`] only the
/// diagonal is used, which ignores overlap correlations.
pub(super) fn intensity_esds(
    positions: &[f64],
    observed: &[PatternPoint],
    options: &LebailOptions,
    weights: &[f64],
    variance: f64,
) -> Vec<f64> {
    let count = positions.len();
    if count == 0 || !variance.is_finite() || variance < 0.0 {
        return vec![0.0; count];
    }
    let mut normal = vec![vec![0.0; count]; count];
    for (row_index, point) in observed.iter().enumerate() {
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
            normal[first][first] += weight * row[first] * row[first];
            for second in (first + 1)..count {
                let value = weight * row[first] * row[second];
                normal[first][second] += value;
                normal[second][first] += value;
            }
        }
    }
    if count <= MAX_EXACT_ESD {
        let scale = normal
            .iter()
            .flat_map(|row| row.iter())
            .fold(0.0_f64, |value, entry| value.max(entry.abs()))
            .max(1.0);
        for (index, row) in normal.iter_mut().enumerate() {
            if let Some(value) = row.get_mut(index) {
                *value += 1e-12 * scale;
            }
        }
        if let Some(covariance) = invert(&normal) {
            return (0..count)
                .map(|index| {
                    covariance
                        .get(index)
                        .and_then(|row| row.get(index))
                        .map_or(0.0, |value| (value * variance).max(0.0).sqrt())
                })
                .collect();
        }
    }
    (0..count)
        .map(|index| {
            let diagonal = normal
                .get(index)
                .and_then(|row| row.get(index))
                .copied()
                .unwrap_or(0.0);
            if diagonal > 0.0 {
                (variance / diagonal).max(0.0).sqrt()
            } else {
                0.0
            }
        })
        .collect()
}

/// Inverts a square matrix, or returns `None` when it is singular.
fn invert(matrix: &[Vec<f64>]) -> Option<Vec<Vec<f64>>> {
    let n = matrix.len();
    if n == 0 || matrix.iter().any(|row| row.len() != n) {
        return None;
    }
    let mut augmented: Vec<Vec<f64>> = matrix
        .iter()
        .enumerate()
        .map(|(row, values)| {
            let mut line = values.clone();
            line.extend((0..n).map(|column| if column == row { 1.0 } else { 0.0 }));
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
            for (slot, value) in line.iter_mut().zip(pivot_row.iter()) {
                *slot -= factor * value;
            }
        }
    }
    Some(
        augmented
            .iter()
            .map(|row| row[n..].to_vec())
            .collect::<Vec<_>>(),
    )
}

pub(super) fn solve(matrix: &[Vec<f64>], rhs: &[f64]) -> Option<Vec<f64>> {
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
