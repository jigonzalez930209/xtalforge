#![allow(clippy::indexing_slicing, clippy::needless_range_loop)] // fixed-size numeric kernels

//! Equation-of-state fits (Task 9.6.2).
//!
//! Birch-Murnaghan (third order), Vinet and Murnaghan energy-volume fits with
//! a Levenberg-Marquardt solver and parameter uncertainties. Energies are eV,
//! volumes cubic Angstrom and the reported bulk modulus GPa.

use serde::{Deserialize, Serialize};

use crate::Error;

/// One energy-volume point.
#[derive(Debug, Clone, Copy, PartialEq, Serialize, Deserialize)]
pub struct EosPoint {
    /// Volume in cubic Angstrom.
    pub volume: f64,
    /// Energy in eV.
    pub energy: f64,
}

/// Equation-of-state model.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
#[serde(rename_all = "snake_case")]
pub enum EosKind {
    /// Third-order Birch-Murnaghan.
    BirchMurnaghan,
    /// Vinet (universal) equation.
    Vinet,
    /// Murnaghan equation.
    Murnaghan,
}

/// A fitted equation of state.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct EosFit {
    /// Model that was fitted.
    pub kind: EosKind,
    /// Equilibrium volume in cubic Angstrom.
    pub v0: f64,
    /// Equilibrium energy in eV.
    pub e0: f64,
    /// Bulk modulus in GPa.
    pub b0: f64,
    /// Pressure derivative of the bulk modulus.
    pub b0_prime: f64,
    /// Root mean square residual in eV.
    pub rms: f64,
    /// Standard errors of `[v0, e0, b0, b0']` in the reported units.
    pub errors: [f64; 4],
    /// Number of points used.
    pub n_points: usize,
    /// Predicted energy at each input volume.
    pub predicted: Vec<f64>,
}

const GPA_PER_EV_A3: f64 = 160.217_662_08;

/// Fits an equation of state to energy-volume points.
///
/// # Errors
///
/// Returns [`Error::InvalidInput`] for fewer than four points, non-finite or
/// non-positive volumes, duplicate volumes or a singular fit.
pub fn fit_eos(points: &[EosPoint], kind: EosKind) -> Result<EosFit, Error> {
    if points.len() < 4 {
        return Err(Error::InvalidInput(
            "an equation-of-state fit needs at least four points".to_owned(),
        ));
    }
    for point in points {
        if !point.volume.is_finite() || point.volume <= 0.0 || !point.energy.is_finite() {
            return Err(Error::InvalidInput(
                "EOS points need positive finite volumes and finite energies".to_owned(),
            ));
        }
    }
    let mut sorted = points.to_vec();
    sorted.sort_by(|first, second| {
        first
            .volume
            .partial_cmp(&second.volume)
            .unwrap_or(std::cmp::Ordering::Equal)
    });
    for window in sorted.windows(2) {
        if let [first, second] = window
            && (first.volume - second.volume).abs() < 1e-9
        {
            return Err(Error::InvalidInput(
                "EOS points contain duplicate volumes".to_owned(),
            ));
        }
    }
    let minimum = sorted
        .iter()
        .enumerate()
        .min_by(|first, second| {
            first
                .1
                .energy
                .partial_cmp(&second.1.energy)
                .unwrap_or(std::cmp::Ordering::Equal)
        })
        .map(|(index, _)| index)
        .unwrap_or(0);
    let v0 = sorted.get(minimum).map_or(1.0, |point| point.volume);
    let e0 = sorted.get(minimum).map_or(0.0, |point| point.energy);
    let curvature = second_derivative(&sorted, minimum).unwrap_or(1.0);
    let mut parameters = [v0, e0, (v0 * curvature).max(1e-4), 4.0];
    let mut lambda = 1e-3;
    let mut previous = residual_sum(&sorted, kind, &parameters);
    for _ in 0..400 {
        let jacobian = jacobian(&sorted, kind, &parameters);
        let residual = residuals(&sorted, kind, &parameters);
        let normal = normal_matrix(&jacobian, lambda);
        let gradient = {
            let mut vector = [0.0; 4];
            for (row, values) in jacobian.iter().enumerate() {
                let Some(value) = residual.get(row) else {
                    continue;
                };
                for column in 0..4 {
                    vector[column] -= values[column] * value;
                }
            }
            vector
        };
        let Some(delta) = solve4(normal, gradient) else {
            lambda *= 10.0;
            if lambda > 1e12 {
                break;
            }
            continue;
        };
        let candidate = [
            parameters[0] + delta[0],
            parameters[1] + delta[1],
            parameters[2] + delta[2],
            parameters[3] + delta[3],
        ];
        if candidate[0] <= 0.0 || candidate[2] <= 0.0 || !candidate.iter().all(|v| v.is_finite()) {
            lambda *= 10.0;
            if lambda > 1e12 {
                break;
            }
            continue;
        }
        let sum = residual_sum(&sorted, kind, &candidate);
        if sum < previous {
            parameters = candidate;
            if (previous - sum).abs() <= 1e-12 * previous.max(1.0) {
                break;
            }
            previous = sum;
            lambda = (lambda / 3.0).max(1e-12);
        } else {
            lambda *= 10.0;
            if lambda > 1e12 {
                break;
            }
        }
    }
    if !parameters.iter().all(|value| value.is_finite()) || parameters[0] <= 0.0 {
        return Err(Error::InvalidInput(
            "the equation-of-state fit did not converge".to_owned(),
        ));
    }
    let residual = residuals(&sorted, kind, &parameters);
    let rms =
        (residual.iter().map(|value| value * value).sum::<f64>() / residual.len() as f64).sqrt();
    let jacobian = jacobian(&sorted, kind, &parameters);
    let normal = normal_matrix(&jacobian, 0.0);
    let errors = covariance_errors(normal, rms, residual.len()).unwrap_or([0.0; 4]);
    let predicted = sorted
        .iter()
        .map(|point| eos_energy(kind, &parameters, point.volume))
        .collect();
    Ok(EosFit {
        kind,
        v0: parameters[0],
        e0: parameters[1],
        b0: parameters[2] * GPA_PER_EV_A3,
        b0_prime: parameters[3],
        rms,
        errors: [errors[0], errors[1], errors[2] * GPA_PER_EV_A3, errors[3]],
        n_points: sorted.len(),
        predicted,
    })
}

/// Energy of one volume for the given parameters `[v0, e0, b0, b0']`, with
/// `b0` in eV/Angstrom^3.
#[must_use]
pub fn eos_energy(kind: EosKind, parameters: &[f64; 4], volume: f64) -> f64 {
    let [v0, e0, b0, b0_prime] = *parameters;
    if v0 <= 0.0 || volume <= 0.0 {
        return f64::NAN;
    }
    match kind {
        EosKind::BirchMurnaghan => {
            let eta = (v0 / volume).powf(1.0 / 3.0);
            let strain = eta * eta - 1.0;
            e0 + 9.0 * v0 * b0 / 16.0
                * (strain.powi(3) * b0_prime + strain.powi(2) * (6.0 - 4.0 * eta * eta))
        }
        EosKind::Vinet => {
            let x = (volume / v0).powf(1.0 / 3.0);
            let factor = 3.0 * (b0_prime - 1.0) * (x - 1.0) / 2.0;
            e0 + 2.0 * b0 * v0 / (b0_prime - 1.0).powi(2)
                * (2.0 - (5.0 + 3.0 * b0_prime * (x - 1.0) - 3.0 * x) * (-factor).exp())
        }
        EosKind::Murnaghan => {
            e0 + b0 * volume / b0_prime * ((v0 / volume).powf(b0_prime) / (b0_prime - 1.0) + 1.0)
                - b0 * v0 / (b0_prime - 1.0)
        }
    }
}

fn residuals(points: &[EosPoint], kind: EosKind, parameters: &[f64; 4]) -> Vec<f64> {
    points
        .iter()
        .map(|point| eos_energy(kind, parameters, point.volume) - point.energy)
        .collect()
}

fn residual_sum(points: &[EosPoint], kind: EosKind, parameters: &[f64; 4]) -> f64 {
    residuals(points, kind, parameters)
        .iter()
        .map(|value| value * value)
        .sum()
}

fn jacobian(points: &[EosPoint], kind: EosKind, parameters: &[f64; 4]) -> Vec<[f64; 4]> {
    let mut rows = Vec::with_capacity(points.len());
    for point in points {
        let mut row = [0.0; 4];
        for column in 0..4 {
            let step = (parameters[column].abs() * 1e-6).max(1e-6);
            let mut plus = *parameters;
            let mut minus = *parameters;
            plus[column] += step;
            minus[column] -= step;
            row[column] = (eos_energy(kind, &plus, point.volume)
                - eos_energy(kind, &minus, point.volume))
                / (2.0 * step);
        }
        rows.push(row);
    }
    rows
}

fn normal_matrix(jacobian: &[[f64; 4]], lambda: f64) -> [[f64; 4]; 4] {
    let mut normal = [[0.0; 4]; 4];
    for row in jacobian {
        for i in 0..4 {
            for j in 0..4 {
                normal[i][j] += row[i] * row[j];
            }
        }
    }
    for (index, diagonal) in normal.iter_mut().enumerate() {
        if let Some(value) = diagonal.get_mut(index) {
            *value += lambda * value.abs().max(1.0);
        }
    }
    normal
}

fn covariance_errors(normal: [[f64; 4]; 4], rms: f64, count: usize) -> Option<[f64; 4]> {
    if count <= 4 {
        return None;
    }
    let inverse = inverse4(normal)?;
    let variance = rms * rms * count as f64 / (count as f64 - 4.0);
    Some([
        (variance * inverse[0][0].abs()).sqrt(),
        (variance * inverse[1][1].abs()).sqrt(),
        (variance * inverse[2][2].abs()).sqrt(),
        (variance * inverse[3][3].abs()).sqrt(),
    ])
}

fn second_derivative(points: &[EosPoint], index: usize) -> Option<f64> {
    let previous = points.get(index.checked_sub(1)?)?;
    let current = points.get(index)?;
    let next = points.get(index + 1)?;
    let h1 = current.volume - previous.volume;
    let h2 = next.volume - current.volume;
    if h1 <= 0.0 || h2 <= 0.0 {
        return None;
    }
    let f1 = (current.energy - previous.energy) / h1;
    let f2 = (next.energy - current.energy) / h2;
    Some(2.0 * (f2 - f1) / (h1 + h2))
}

fn solve4(matrix: [[f64; 4]; 4], vector: [f64; 4]) -> Option<[f64; 4]> {
    let mut augmented = [[0.0; 5]; 4];
    for row in 0..4 {
        for column in 0..4 {
            augmented[row][column] = matrix[row][column];
        }
        augmented[row][4] = vector[row];
    }
    for pivot in 0..4 {
        let mut best = pivot;
        for row in (pivot + 1)..4 {
            if augmented[row][pivot].abs() > augmented[best][pivot].abs() {
                best = row;
            }
        }
        if augmented[best][pivot].abs() < 1e-14 {
            return None;
        }
        augmented.swap(pivot, best);
        let divisor = augmented[pivot][pivot];
        for column in pivot..5 {
            augmented[pivot][column] /= divisor;
        }
        for row in 0..4 {
            if row == pivot {
                continue;
            }
            let factor = augmented[row][pivot];
            for column in pivot..5 {
                augmented[row][column] -= factor * augmented[pivot][column];
            }
        }
    }
    Some([
        augmented[0][4],
        augmented[1][4],
        augmented[2][4],
        augmented[3][4],
    ])
}

fn inverse4(matrix: [[f64; 4]; 4]) -> Option<[[f64; 4]; 4]> {
    let mut inverse = [[0.0; 4]; 4];
    for column in 0..4 {
        let mut unit = [0.0; 4];
        unit[column] = 1.0;
        let solution = solve4(matrix, unit)?;
        for row in 0..4 {
            inverse[row][column] = solution[row];
        }
    }
    Some(inverse)
}
