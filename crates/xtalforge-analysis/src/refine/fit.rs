//! Linear algebra and analytic scale/background fit for [`super::refine_phase`].

use xtalforge_core::Error;
use xtalforge_io::xy::PatternPoint;

/// Fits `scale * curve + background(t)` and returns the coefficients and the
/// model values. Background terms are `1`, `t`, `t^2` with `t` centred and
/// scaled to `[-1, 1]` over the observed range for conditioning.
pub(super) fn fit_scale_background(
    curve: &[f64],
    observed: &[PatternPoint],
    terms: usize,
) -> Result<(f64, Vec<f64>, Vec<f64>), Error> {
    let first = observed.first().map_or(0.0, |point| point.two_theta);
    let last = observed.last().map_or(first, |point| point.two_theta);
    let center = 0.5 * (first + last);
    let half = 0.5 * (last - first);
    let span = if half > 0.0 { half } else { 1.0 };
    let columns = terms.clamp(1, 3) + 1;
    let mut normal = vec![vec![0.0; columns]; columns];
    let mut rhs = vec![0.0; columns];
    for (unit, point) in curve.iter().zip(observed.iter()) {
        let t = (point.two_theta - center) / span;
        for row in 0..columns {
            let value = basis(*unit, t, row);
            if let Some(slot) = rhs.get_mut(row) {
                *slot += value * point.intensity;
            }
            for column in 0..columns {
                let other = basis(*unit, t, column);
                if let Some(slot) = normal.get_mut(row).and_then(|line| line.get_mut(column)) {
                    *slot += value * other;
                }
            }
        }
    }
    for index in 0..columns {
        if let Some(slot) = normal.get_mut(index).and_then(|row| row.get_mut(index)) {
            *slot += 1e-12 * slot.abs().max(1.0);
        }
    }
    let Some(coefficients) = solve_linear(&normal, &rhs) else {
        return Err(Error::InvalidInput(
            "the scale and background normal matrix is singular".to_owned(),
        ));
    };
    let scale = coefficients.first().copied().unwrap_or(0.0);
    let background: Vec<f64> = (1..columns)
        .map(|index| coefficients.get(index).copied().unwrap_or(0.0))
        .collect();
    if !scale.is_finite() || background.iter().any(|value| !value.is_finite()) {
        return Err(Error::InvalidInput(
            "the scale and background fit produced non-finite values".to_owned(),
        ));
    }
    let mut model = Vec::with_capacity(observed.len());
    for (unit, point) in curve.iter().zip(observed.iter()) {
        let t = (point.two_theta - center) / span;
        model.push(scale * unit + background_value(&background, t));
    }
    Ok((scale, background, model))
}

fn basis(unit: f64, t: f64, index: usize) -> f64 {
    match index {
        0 => unit,
        1 => 1.0,
        2 => t,
        3 => t * t,
        _ => 0.0,
    }
}

fn background_value(background: &[f64], t: f64) -> f64 {
    let mut value = 0.0;
    let mut power = 1.0;
    for coefficient in background {
        value += coefficient * power;
        power *= t;
    }
    value
}

pub(super) fn weighted_rwp(observed: &[PatternPoint], model: &[f64]) -> Result<f64, Error> {
    let mut numerator = 0.0;
    let mut denominator = 0.0;
    for (point, value) in observed.iter().zip(model.iter()) {
        let weight = 1.0 / point.intensity.abs().max(1.0);
        numerator += weight * (point.intensity - value).powi(2);
        denominator += weight * point.intensity.powi(2);
    }
    if !denominator.is_finite() || denominator <= 0.0 {
        return Err(Error::InvalidInput(
            "the observed pattern has no intensity to fit".to_owned(),
        ));
    }
    Ok((numerator / denominator).sqrt())
}

pub(super) fn weighted_residual(observed: &[PatternPoint], model: &[f64]) -> Vec<f64> {
    observed
        .iter()
        .zip(model.iter())
        .map(|(point, value)| {
            let root = 1.0 / point.intensity.abs().max(1.0).sqrt();
            root * (value - point.intensity)
        })
        .collect()
}

pub(super) fn normal_matrix(jacobian: &[Vec<f64>]) -> Vec<Vec<f64>> {
    let count = jacobian.first().map_or(0, Vec::len);
    let mut normal = vec![vec![0.0; count]; count];
    for row in jacobian {
        for first in 0..count {
            let Some(a) = row.get(first).copied() else {
                continue;
            };
            for second in 0..count {
                let Some(b) = row.get(second).copied() else {
                    continue;
                };
                if let Some(slot) = normal.get_mut(first).and_then(|line| line.get_mut(second)) {
                    *slot += a * b;
                }
            }
        }
    }
    normal
}

pub(super) fn gradient(jacobian: &[Vec<f64>], residual: &[f64]) -> Vec<f64> {
    let count = jacobian.first().map_or(0, Vec::len);
    let mut gradient = vec![0.0; count];
    for (row, value) in jacobian.iter().zip(residual.iter()) {
        for column in 0..count {
            let Some(entry) = row.get(column).copied() else {
                continue;
            };
            if let Some(slot) = gradient.get_mut(column) {
                *slot += entry * value;
            }
        }
    }
    gradient
}

pub(super) fn solve_linear(matrix: &[Vec<f64>], rhs: &[f64]) -> Option<Vec<f64>> {
    let n = matrix.len();
    if n == 0 || rhs.len() != n {
        return None;
    }
    let mut augmented: Vec<Vec<f64>> = matrix
        .iter()
        .zip(rhs.iter())
        .map(|(row, value)| {
            let mut line = row.clone();
            line.push(*value);
            line
        })
        .collect();
    let entry = |matrix: &[Vec<f64>], row: usize, column: usize| {
        matrix
            .get(row)
            .and_then(|line| line.get(column))
            .copied()
            .unwrap_or(0.0)
    };
    for column in 0..n {
        let pivot = (column..n).max_by(|a, b| {
            entry(&augmented, *a, column)
                .abs()
                .partial_cmp(&entry(&augmented, *b, column).abs())
                .unwrap_or(std::cmp::Ordering::Equal)
        })?;
        if entry(&augmented, pivot, column).abs() < 1e-14 {
            return None;
        }
        augmented.swap(column, pivot);
        let divisor = entry(&augmented, column, column);
        if let Some(row) = augmented.get_mut(column) {
            for value in row.iter_mut() {
                *value /= divisor;
            }
        }
        let pivot_row = augmented.get(column)?.clone();
        for (row, line) in augmented.iter_mut().enumerate() {
            if row == column {
                continue;
            }
            let factor = line.get(column).copied().unwrap_or(0.0);
            for (slot, value) in line.iter_mut().zip(pivot_row.iter()) {
                *slot -= factor * value;
            }
        }
    }
    Some(
        augmented
            .iter()
            .map(|row| row.get(n).copied().unwrap_or(0.0))
            .collect(),
    )
}
