//! Cubic powder indexing from peak positions (Task 7.5).

use serde::{Deserialize, Serialize};
use xtalforge_core::Error;

use crate::peaks::Peak;

use super::IndexedPeak;

/// A cubic indexing candidate.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct CubicIndex {
    /// Cubic lattice parameter in Angstrom.
    pub a: f64,
    /// Indexed peaks.
    pub peaks: Vec<IndexedPeak>,
    /// Fraction of peaks that found an integer hkl.
    pub completeness: f64,
    /// Mean absolute difference between observed and calculated 2-theta.
    pub mean_error: f64,
}

/// Indexes a peak list assuming a cubic cell.
///
/// Tries every hkl up to `max_index` for the first peak, derives `a`, and
/// matches the remaining peaks within `tolerance` degrees. The best candidate
/// is the one with the most matches; among near-equal errors the smallest cell
/// wins, so a face-centred pattern indexes on its own cubic cell and not on a
/// larger multiple of it.
///
/// # Errors
///
/// Returns [`Error::InvalidInput`] for fewer than two peaks or a non-positive
/// wavelength.
pub fn index_cubic(
    peaks: &[Peak],
    wavelength: f64,
    max_index: i32,
    tolerance: f64,
) -> Result<Option<CubicIndex>, Error> {
    if peaks.len() < 2 {
        return Err(Error::InvalidInput(
            "indexing needs at least two peaks".to_owned(),
        ));
    }
    if !wavelength.is_finite() || wavelength <= 0.0 {
        return Err(Error::InvalidInput(
            "wavelength must be positive".to_owned(),
        ));
    }
    let mut candidates: Vec<CubicIndex> = Vec::new();
    for h in 1..=max_index {
        for k in 0..=max_index {
            for l in 0..=max_index {
                let sum = h * h + k * k + l * l;
                if sum == 0 {
                    continue;
                }
                let Some(first) = peaks.first() else {
                    continue;
                };
                let theta = first.two_theta.to_radians() / 2.0;
                let d = wavelength / (2.0 * theta.sin());
                let a = d * (sum as f64).sqrt();
                let mut indexed = Vec::new();
                let mut error = 0.0;
                for peak in peaks {
                    if let Some(hkl) =
                        match_cubic(a, peak.two_theta, wavelength, max_index, tolerance)
                    {
                        let calculated = two_theta_cubic(a, hkl, wavelength);
                        error += (calculated - peak.two_theta).abs();
                        indexed.push(IndexedPeak {
                            hkl,
                            two_theta: peak.two_theta,
                        });
                    }
                }
                let completeness = indexed.len() as f64 / peaks.len() as f64;
                let mean_error = if indexed.is_empty() {
                    f64::INFINITY
                } else {
                    error / indexed.len() as f64
                };
                let candidate = CubicIndex {
                    a,
                    peaks: indexed,
                    completeness,
                    mean_error,
                };
                candidates.push(candidate);
            }
        }
    }
    // Several cells of the same family index the same lines (a, a sqrt 2,
    // a sqrt 3, ... when the first line is assigned a larger h^2 + k^2 + l^2).
    // Keep the best completeness, and among the candidates whose error is
    // within a tenth of the tolerance of the smallest, the smallest cell.
    let completeness = candidates
        .iter()
        .map(|candidate| candidate.completeness)
        .fold(0.0, f64::max);
    let least_error = candidates
        .iter()
        .filter(|candidate| candidate.completeness >= completeness)
        .map(|candidate| candidate.mean_error)
        .fold(f64::INFINITY, f64::min);
    let best = candidates
        .into_iter()
        .filter(|candidate| {
            candidate.completeness >= completeness
                && candidate.mean_error <= least_error + 0.1 * tolerance
        })
        .min_by(|first, second| first.a.total_cmp(&second.a));
    Ok(best)
}

fn match_cubic(
    a: f64,
    two_theta: f64,
    wavelength: f64,
    max_index: i32,
    tolerance: f64,
) -> Option<[i32; 3]> {
    let mut best: Option<([i32; 3], f64)> = None;
    for h in 0..=max_index {
        for k in 0..=max_index {
            for l in 0..=max_index {
                if h == 0 && k == 0 && l == 0 {
                    continue;
                }
                let hkl = [h, k, l];
                let difference = (two_theta_cubic(a, hkl, wavelength) - two_theta).abs();
                if difference <= tolerance && best.is_none_or(|(_, value)| difference < value) {
                    best = Some((hkl, difference));
                }
            }
        }
    }
    best.map(|(hkl, _)| hkl)
}

fn two_theta_cubic(a: f64, hkl: [i32; 3], wavelength: f64) -> f64 {
    let sum = (hkl[0] * hkl[0] + hkl[1] * hkl[1] + hkl[2] * hkl[2]) as f64;
    let d = a / sum.sqrt();
    let sine = wavelength / (2.0 * d);
    if sine.abs() >= 1.0 {
        f64::NAN
    } else {
        2.0 * sine.asin().to_degrees()
    }
}
