//! Small non-linear least-squares engine (Task 7.4.1).

#![allow(clippy::indexing_slicing)] // local matrix indices are bounded

use serde::{Deserialize, Serialize};
use xtalforge_core::Error;

/// Result of a least-squares fit.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct FitResult {
    /// Best parameters.
    pub parameters: Vec<f64>,
    /// Residual vector at the solution.
    pub residuals: Vec<f64>,
    /// Sum of squared residuals.
    pub chi_squared: f64,
    /// Iterations used.
    pub iterations: usize,
    /// Whether the parameter step fell below the tolerance.
    pub converged: bool,
}

/// Options for the Levenberg-Marquardt fit.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct FitOptions {
    /// Maximum iterations.
    pub max_iterations: usize,
    /// Relative step tolerance.
    pub tolerance: f64,
    /// Initial damping.
    pub lambda: f64,
}

impl Default for FitOptions {
    fn default() -> Self {
        Self {
            max_iterations: 100,
            tolerance: 1e-8,
            lambda: 1e-3,
        }
    }
}

/// Fits parameters to a residual function with finite-difference Jacobians.
///
/// # Errors
///
/// Returns [`Error::InvalidInput`] for empty parameters, an empty residual or
/// non-finite values.
pub fn least_squares<F>(
    initial: &[f64],
    options: &FitOptions,
    residuals: F,
) -> Result<FitResult, Error>
where
    F: Fn(&[f64]) -> Result<Vec<f64>, Error>,
{
    if initial.is_empty() {
        return Err(Error::InvalidInput(
            "a fit needs at least one parameter".to_owned(),
        ));
    }
    if !options.tolerance.is_finite() || options.tolerance <= 0.0 {
        return Err(Error::InvalidInput(
            "fit tolerance must be positive".to_owned(),
        ));
    }
    let mut parameters = initial.to_vec();
    let mut current = residuals(&parameters)?;
    if current.is_empty() {
        return Err(Error::InvalidInput(
            "the residual function returned no values".to_owned(),
        ));
    }
    let mut chi_squared = sum_squares(&current);
    let mut lambda = options.lambda.max(1e-12);
    let mut converged = false;
    let mut iterations = 0;
    for step in 0..options.max_iterations {
        iterations = step + 1;
        let count = parameters.len();
        let mut jacobian = vec![vec![0.0; count]; current.len()];
        for column in 0..count {
            let mut shifted = parameters.clone();
            let delta = 1e-6 * (1.0 + parameters[column].abs());
            shifted[column] += delta;
            let shifted_residuals = residuals(&shifted)?;
            for (row, value) in shifted_residuals.iter().enumerate() {
                if let Some(slot) = jacobian.get_mut(row).and_then(|line| line.get_mut(column)) {
                    *slot = (value - current[row]) / delta;
                }
            }
        }
        let normal = normal_matrix(&jacobian);
        let gradient = gradient(&jacobian, &current);
        let mut damped = normal;
        for (index, row) in damped.iter_mut().enumerate() {
            if let Some(value) = row.get_mut(index) {
                *value += lambda * (1.0 + value.abs());
            }
        }
        let Some(delta) = solve(&damped, &gradient) else {
            lambda *= 10.0;
            continue;
        };
        let candidate: Vec<f64> = parameters
            .iter()
            .zip(delta.iter())
            .map(|(value, step)| value - step)
            .collect();
        let next = residuals(&candidate)?;
        let next_chi = sum_squares(&next);
        if next_chi < chi_squared {
            let change = delta.iter().map(|value| value * value).sum::<f64>().sqrt();
            parameters = candidate;
            current = next;
            chi_squared = next_chi;
            lambda = (lambda / 3.0).max(1e-12);
            if change < options.tolerance * (1.0 + norm(&parameters)) {
                converged = true;
                break;
            }
        } else {
            lambda *= 5.0;
            if lambda > 1e12 {
                break;
            }
        }
    }
    Ok(FitResult {
        parameters,
        residuals: current,
        chi_squared,
        iterations,
        converged,
    })
}

fn sum_squares(values: &[f64]) -> f64 {
    values.iter().map(|value| value * value).sum()
}

fn norm(values: &[f64]) -> f64 {
    values.iter().map(|value| value * value).sum::<f64>().sqrt()
}

fn normal_matrix(jacobian: &[Vec<f64>]) -> Vec<Vec<f64>> {
    let count = jacobian.first().map_or(0, Vec::len);
    let mut normal = vec![vec![0.0; count]; count];
    for row in jacobian {
        for (first, line) in normal.iter_mut().enumerate() {
            let Some(a) = row.get(first) else {
                continue;
            };
            for (second, slot) in line.iter_mut().enumerate() {
                if let Some(b) = row.get(second) {
                    *slot += a * b;
                }
            }
        }
    }
    normal
}

fn gradient(jacobian: &[Vec<f64>], residuals: &[f64]) -> Vec<f64> {
    let count = jacobian.first().map_or(0, Vec::len);
    let mut gradient = vec![0.0; count];
    for (row, residual) in jacobian.iter().zip(residuals.iter()) {
        for (column, value) in row.iter().enumerate() {
            if let Some(slot) = gradient.get_mut(column) {
                *slot += value * residual;
            }
        }
    }
    gradient
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
            matrix[*a][column]
                .abs()
                .partial_cmp(&matrix[*b][column].abs())
                .unwrap_or(std::cmp::Ordering::Equal)
        })?;
        if matrix[pivot][column].abs() < 1e-15 {
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
