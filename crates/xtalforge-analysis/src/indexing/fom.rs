//! Figures of merit for a candidate cell (Task 7.5.2).
//!
//! de Wolff's `M20 = Q20 / (2 * <|dQ|> * N20)` measures how closely a cell
//! predicts the first twenty observed lines: `Q20` is the reciprocal squared
//! of the twentieth observed reflection, `N20` the number of calculated
//! reflections in one hemisphere with `Q <= Q20`, and `<|dQ|>` the mean
//! absolute difference over the indexed lines. Smith and Snyder's
//! `F_N = N / (<|d2theta|> * N_possible)` divides the same idea by the number
//! of calculated lines that could appear before the last observed one.
//!
//! Both figures count reflections in one hemisphere (the `2` in `M20` fills
//! the sphere), and both use the mean over the matched lines only. A cell that
//! indexes few lines cannot buy a high figure with one very good match,
//! because `N20` and `N_possible` grow with the cell volume.

#![allow(clippy::indexing_slicing)] // bounded local windows

use serde::{Deserialize, Serialize};
use xtalforge_core::Error;

use crate::peaks::Peak;

use super::{CellMetric, IndexedPeak};

/// Highest index enumerated when predicting lines for the figures of merit.
const MAX_FOM_INDEX: i32 = 8;
/// Floor that keeps an exactly predicted pattern finite.
const EPSILON_FLOOR: f64 = 1e-12;

/// Figures of merit of one cell against one peak list.
#[derive(Debug, Clone, Copy, PartialEq, Serialize, Deserialize)]
pub struct FiguresOfMerit {
    /// de Wolff M20 (higher is better).
    pub m20: f64,
    /// Smith-Snyder F_N (higher is better).
    pub f_n: f64,
    /// Number of calculated reflections in one hemisphere up to `Q20`.
    pub n20: usize,
    /// Number of calculated reflections up to the last observed line.
    pub n_possible: usize,
    /// Mean absolute `dQ` over the matched lines.
    pub mean_delta_q: f64,
    /// Mean absolute two-theta difference over the matched lines.
    pub mean_delta_two_theta: f64,
    /// Fraction of the predicted lines below `Q20` that the observed peaks
    /// explain. A cell that predicts many unobserved lines scores below one.
    pub coverage: f64,
}

/// A candidate cell matched against a peak list.
#[derive(Debug, Clone, PartialEq)]
pub(crate) struct Evaluation {
    /// Matched peaks in observed order.
    pub indexed: Vec<IndexedPeak>,
    /// Fraction of observed peaks that matched.
    pub completeness: f64,
    /// Mean absolute observed minus calculated two-theta.
    pub mean_error: f64,
    /// Figures of merit.
    pub figures: FiguresOfMerit,
}

/// Reciprocal squared `Q = 1/d^2` of a peak position.
///
/// # Errors
///
/// Returns [`Error::InvalidInput`] for a non-positive two-theta or wavelength.
pub fn peak_q(peak: &Peak, wavelength: f64) -> Result<f64, Error> {
    if !peak.two_theta.is_finite() || peak.two_theta <= 0.0 {
        return Err(Error::InvalidInput(format!(
            "peak position must be positive, got {}",
            peak.two_theta
        )));
    }
    if !wavelength.is_finite() || wavelength <= 0.0 {
        return Err(Error::InvalidInput(
            "wavelength must be positive and finite".to_owned(),
        ));
    }
    let theta = peak.two_theta.to_radians() / 2.0;
    let sine = theta.sin();
    Ok((2.0 * sine / wavelength).powi(2))
}

/// Two-theta of a reciprocal squared value.
#[must_use]
pub fn two_theta_of_q(q: f64, wavelength: f64) -> f64 {
    let sine = (q.sqrt() * wavelength / 2.0).clamp(0.0, 1.0);
    2.0 * sine.asin().to_degrees()
}

/// Largest |h|, |k| or |l| that can reach `max_q` in this cell.
///
/// The bound uses the reciprocal axis lengths, so a large cell with a small
/// `max_q` is searched over a handful of indices instead of a fixed cube.
#[must_use]
pub(crate) fn index_limit(cell: &CellMetric, max_q: f64) -> i32 {
    let alpha = cell.alpha.to_radians();
    let beta = cell.beta.to_radians();
    let gamma = cell.gamma.to_radians();
    let volume = cell.volume();
    if volume <= 0.0 || max_q <= 0.0 || !max_q.is_finite() {
        return 0;
    }
    let d100 = volume / (cell.b * cell.c * alpha.sin());
    let d010 = volume / (cell.a * cell.c * beta.sin());
    let d001 = volume / (cell.a * cell.b * gamma.sin());
    let reach = max_q.sqrt();
    let limit = [d100, d010, d001]
        .iter()
        .map(|spacing| (reach * spacing).ceil() as i32)
        .max()
        .unwrap_or(1);
    limit.clamp(1, MAX_FOM_INDEX)
}

/// Predicted reflections with `Q <= max_q`, sorted by `Q`.
///
/// Only one of each `+/-hkl` pair is generated, which is why `M20` carries the
/// factor two.
#[must_use]
pub fn predict_reflections(cell: &CellMetric, max_q: f64) -> Vec<(f64, [i32; 3])> {
    let limit = index_limit(cell, max_q);
    let mut lines = Vec::new();
    for h in -limit..=limit {
        for k in -limit..=limit {
            for l in -limit..=limit {
                if h < 0 || (h == 0 && k < 0) || (h == 0 && k == 0 && l <= 0) {
                    continue;
                }
                let hkl = [h, k, l];
                let Some(q) = cell.reciprocal_squared(hkl) else {
                    continue;
                };
                if q <= max_q * (1.0 + 1e-9) {
                    lines.push((q, hkl));
                }
            }
        }
    }
    lines.sort_by(|first, second| first.0.total_cmp(&second.0));
    lines
}

fn nearest_line(lines: &[(f64, [i32; 3])], q: f64) -> Option<(f64, [i32; 3])> {
    if lines.is_empty() {
        return None;
    }
    let index = lines.partition_point(|line| line.0 < q);
    let before = index.checked_sub(1).and_then(|slot| lines.get(slot));
    let after = lines.get(index);
    match (before, after) {
        (Some(first), Some(second)) => {
            if (first.0 - q).abs() <= (second.0 - q).abs() {
                Some(*first)
            } else {
                Some(*second)
            }
        }
        (Some(first), None) => Some(*first),
        (None, Some(second)) => Some(*second),
        (None, None) => None,
    }
}

/// Matches every peak to the nearest predicted reflection within `tolerance`.
#[must_use]
pub fn match_peaks(
    cell: &CellMetric,
    peaks: &[Peak],
    wavelength: f64,
    tolerance: f64,
) -> Vec<IndexedPeak> {
    evaluate_cell(cell, peaks, wavelength, tolerance)
        .map(|evaluation| evaluation.indexed)
        .unwrap_or_default()
}

/// Matches a peak list and computes both figures of merit from one prediction.
pub(crate) fn evaluate_cell(
    cell: &CellMetric,
    peaks: &[Peak],
    wavelength: f64,
    tolerance: f64,
) -> Result<Evaluation, Error> {
    if peaks.is_empty() {
        return Err(Error::InvalidInput(
            "figures of merit need at least one peak".to_owned(),
        ));
    }
    if !tolerance.is_finite() || tolerance <= 0.0 {
        return Err(Error::InvalidInput(
            "figure-of-merit tolerance must be positive".to_owned(),
        ));
    }
    let observed = peaks
        .iter()
        .map(|peak| peak_q(peak, wavelength))
        .collect::<Result<Vec<_>, _>>()?;
    let last_q = observed.iter().copied().fold(0.0, f64::max);
    let predicted = predict_reflections(cell, last_q);
    if predicted.is_empty() {
        return Err(Error::InvalidInput(
            "the cell predicts no reflection inside the observed range".to_owned(),
        ));
    }
    let mut indexed = Vec::new();
    let mut matched: Vec<(f64, f64, f64)> = Vec::new();
    let mut error_sum = 0.0;
    for (position, peak) in peaks.iter().enumerate() {
        let q = observed.get(position).copied().unwrap_or(0.0);
        let Some((line_q, hkl)) = nearest_line(&predicted, q) else {
            continue;
        };
        let calculated = two_theta_of_q(line_q, wavelength);
        let difference = (calculated - peak.two_theta).abs();
        if difference > tolerance {
            continue;
        }
        indexed.push(IndexedPeak {
            hkl,
            two_theta: peak.two_theta,
        });
        error_sum += difference;
        matched.push((q, (line_q - q).abs(), difference));
    }
    let completeness = indexed.len() as f64 / peaks.len() as f64;
    let mean_error = if indexed.is_empty() {
        f64::INFINITY
    } else {
        error_sum / indexed.len() as f64
    };
    matched.sort_by(|first, second| first.0.total_cmp(&second.0));
    let figures = figures_from(&observed, &matched, &predicted)?;
    Ok(Evaluation {
        indexed,
        completeness,
        mean_error,
        figures,
    })
}

/// Computes M20 and F_N of `cell` against `peaks`.
///
/// # Errors
///
/// Returns [`Error::InvalidInput`] for an empty peak list, a non-positive
/// wavelength or a cell that predicts no reflection in range.
pub fn figures_of_merit(
    cell: &CellMetric,
    peaks: &[Peak],
    wavelength: f64,
    tolerance: f64,
) -> Result<FiguresOfMerit, Error> {
    Ok(evaluate_cell(cell, peaks, wavelength, tolerance)?.figures)
}

fn figures_from(
    observed: &[f64],
    matched: &[(f64, f64, f64)],
    predicted: &[(f64, [i32; 3])],
) -> Result<FiguresOfMerit, Error> {
    if matched.is_empty() {
        return Err(Error::InvalidInput(
            "the cell indexes none of the observed peaks".to_owned(),
        ));
    }
    let window = matched.len().min(20);
    let q20 = matched
        .get(window - 1)
        .map(|entry| entry.0)
        .unwrap_or_else(|| observed.last().copied().unwrap_or(0.0));
    let n20 = predicted.iter().filter(|line| line.0 <= q20).count();
    if n20 == 0 {
        return Err(Error::InvalidInput(
            "the cell predicts no reflection below Q20".to_owned(),
        ));
    }
    let count = window as f64;
    let mean_delta_q = (matched
        .iter()
        .take(window)
        .map(|entry| entry.1)
        .sum::<f64>()
        / count)
        .max(EPSILON_FLOOR);
    let mean_delta_two_theta = (matched
        .iter()
        .take(window)
        .map(|entry| entry.2)
        .sum::<f64>()
        / count)
        .max(EPSILON_FLOOR);
    let coverage = (window as f64 / n20 as f64).clamp(0.0, 1.0);
    Ok(FiguresOfMerit {
        m20: q20 / (2.0 * mean_delta_q * n20 as f64),
        f_n: observed.len() as f64 / (mean_delta_two_theta * predicted.len() as f64),
        n20,
        n_possible: predicted.len(),
        mean_delta_q,
        mean_delta_two_theta,
        coverage,
    })
}
