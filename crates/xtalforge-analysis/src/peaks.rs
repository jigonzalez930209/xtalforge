//! Peak finding on measured XY patterns (Tasks 7.2 and 7.3).

#![allow(clippy::indexing_slicing)] // local window indices are bounded

use serde::{Deserialize, Serialize};
use xtalforge_core::Error;
use xtalforge_io::xy::PatternPoint;

/// One detected peak.
#[derive(Debug, Clone, Copy, PartialEq, Serialize, Deserialize)]
pub struct Peak {
    /// Peak position in degrees two-theta.
    pub two_theta: f64,
    /// Peak height.
    pub intensity: f64,
    /// Full width at half maximum in degrees, estimated from the neighbours.
    pub fwhm: f64,
}

/// Options for peak finding.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct PeakOptions {
    /// Minimum height.
    pub min_height: f64,
    /// Minimum separation in degrees.
    pub min_distance: f64,
    /// Number of smoothing passes (0 disables smoothing).
    pub smoothing: usize,
}

impl Default for PeakOptions {
    fn default() -> Self {
        Self {
            min_height: 0.0,
            min_distance: 0.2,
            smoothing: 0,
        }
    }
}

/// Finds local maxima in a pattern.
///
/// # Errors
///
/// Returns [`Error::InvalidInput`] for fewer than three points or invalid
/// options.
pub fn find_peaks(points: &[PatternPoint], options: &PeakOptions) -> Result<Vec<Peak>, Error> {
    if points.len() < 3 {
        return Err(Error::InvalidInput(
            "peak finding needs at least three points".to_owned(),
        ));
    }
    if !options.min_height.is_finite()
        || !options.min_distance.is_finite()
        || options.min_distance <= 0.0
    {
        return Err(Error::InvalidInput(
            "peak options must be finite and positive".to_owned(),
        ));
    }
    let mut values: Vec<f64> = points.iter().map(|point| point.intensity).collect();
    for _ in 0..options.smoothing {
        values = values
            .windows(3)
            .map(|window| (window[0] + window[1] + window[2]) / 3.0)
            .collect();
    }
    let offset = options.smoothing;
    let mut peaks: Vec<Peak> = Vec::new();
    for index in 1..values.len().saturating_sub(1) {
        let value = values[index];
        if value < options.min_height || value <= values[index - 1] || value < values[index + 1] {
            continue;
        }
        let source = 1 + offset + index;
        let Some(point) = points.get(source) else {
            continue;
        };
        if let Some(last) = peaks.last()
            && point.two_theta - last.two_theta < options.min_distance
        {
            if value > last.intensity {
                peaks.pop();
            } else {
                continue;
            }
        }
        peaks.push(Peak {
            two_theta: point.two_theta,
            intensity: value,
            fwhm: estimate_fwhm(&values, index, points, source),
        });
    }
    Ok(peaks)
}

fn estimate_fwhm(values: &[f64], index: usize, points: &[PatternPoint], source: usize) -> f64 {
    let half = values[index] / 2.0;
    let mut left = index;
    while left > 0 && values[left] > half {
        left -= 1;
    }
    let mut right = index;
    while right + 1 < values.len() && values[right] > half {
        right += 1;
    }
    let left_theta = points
        .get(source - (index - left))
        .map_or(points[source].two_theta, |point| point.two_theta);
    let right_theta = points
        .get(source + (right - index))
        .map_or(points[source].two_theta, |point| point.two_theta);
    (right_theta - left_theta).abs()
}
