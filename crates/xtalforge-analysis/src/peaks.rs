//! Peak finding on measured XY patterns (Tasks 7.2 and 7.3).
//!
//! Detection is a smoothed local maximum above a height threshold; when
//! `background` is on, a SNIP envelope is subtracted first. `refine` fits
//! every detected line with a Gaussian plus a linear background in its own
//! window, which recovers overlapping lines and asymmetric profiles far
//! better than the raw maximum. `fit_peaks` exposes the same fit at named
//! positions, which is what an indexing run wants to refine.

#![allow(clippy::indexing_slicing)] // local window indices are bounded

mod background;
mod whole;

pub use background::background_envelope;
pub use whole::{
    WholeRangeFit, WholeRangeOptions, refine_whole_range, refine_whole_range_with_control,
};

use serde::{Deserialize, Serialize};
use xtalforge_core::Error;
use xtalforge_io::xy::PatternPoint;

use crate::fit::{FitOptions, least_squares};

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
    /// Minimum height, applied to the background-subtracted signal.
    pub min_height: f64,
    /// Minimum separation in degrees.
    pub min_distance: f64,
    /// Number of smoothing passes (0 disables smoothing).
    pub smoothing: usize,
    /// Subtract a SNIP background envelope before detection.
    pub background: bool,
    /// SNIP clipping window in points.
    pub background_window: usize,
    /// Refine every detected line with a Gaussian and a linear background.
    pub refine: bool,
    /// Half-width in degrees of the per-line fit window.
    pub fit_window: f64,
}

impl Default for PeakOptions {
    fn default() -> Self {
        Self {
            min_height: 0.0,
            min_distance: 0.2,
            smoothing: 0,
            background: false,
            background_window: 40,
            refine: false,
            fit_window: 0.5,
        }
    }
}

impl PeakOptions {
    fn validate(&self) -> Result<(), Error> {
        if !self.min_height.is_finite()
            || !self.min_distance.is_finite()
            || self.min_distance <= 0.0
        {
            return Err(Error::InvalidInput(
                "peak options must be finite and positive".to_owned(),
            ));
        }
        if self.background && self.background_window == 0 {
            return Err(Error::InvalidInput(
                "the background window must be positive".to_owned(),
            ));
        }
        if self.refine && (!self.fit_window.is_finite() || self.fit_window <= 0.0) {
            return Err(Error::InvalidInput(
                "the fit window must be positive and finite".to_owned(),
            ));
        }
        Ok(())
    }
}

impl PeakOptions {
    /// Options for a pattern with a background and noise: a SNIP background is
    /// subtracted first, then the lines kept are above 1 % of the strongest one
    /// and above five times the noise (the median absolute deviation of the
    /// background-free signal), after one smoothing pass. `run_automatic` and
    /// `powder index --background` use these.
    #[must_use]
    pub fn noise_limited(points: &[PatternPoint]) -> Self {
        let window = NOISE_LIMITED_BACKGROUND_WINDOW;
        let envelope = background_envelope(points, window);
        let mut signal: Vec<f64> = points
            .iter()
            .zip(envelope.iter())
            .map(|(point, background)| (point.intensity - background).max(0.0))
            .collect();
        let maximum = signal.iter().copied().fold(0.0_f64, f64::max);
        signal.sort_by(f64::total_cmp);
        let median = signal.get(signal.len() / 2).copied().unwrap_or(0.0);
        let mut deviations: Vec<f64> = signal.iter().map(|value| (value - median).abs()).collect();
        deviations.sort_by(f64::total_cmp);
        let noise = 1.4826 * deviations.get(deviations.len() / 2).copied().unwrap_or(0.0);
        Self {
            min_height: (0.01 * maximum).max(5.0 * noise),
            min_distance: (2.0 * median_step(points)).max(0.01),
            smoothing: 1,
            background: true,
            background_window: window,
            ..Self::default()
        }
    }
}

/// SNIP window, in points, of [`PeakOptions::noise_limited`].
const NOISE_LIMITED_BACKGROUND_WINDOW: usize = 40;

fn median_step(points: &[PatternPoint]) -> f64 {
    let mut steps: Vec<f64> = points
        .windows(2)
        .filter_map(|pair| match pair {
            [first, second] if second.two_theta > first.two_theta => {
                Some(second.two_theta - first.two_theta)
            }
            _ => None,
        })
        .collect();
    steps.sort_by(f64::total_cmp);
    steps
        .get(steps.len() / 2)
        .copied()
        .unwrap_or(0.02)
        .max(0.001)
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
    options.validate()?;
    let baseline = if options.background {
        background_envelope(points, options.background_window)
    } else {
        vec![0.0; points.len()]
    };
    let mut values: Vec<f64> = points
        .iter()
        .zip(baseline.iter())
        .map(|(point, background)| (point.intensity - background).max(0.0))
        .collect();
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
    if options.refine && !peaks.is_empty() {
        let positions: Vec<f64> = peaks.iter().map(|peak| peak.two_theta).collect();
        peaks = fit_peaks(points, &positions, options)?;
    }
    Ok(peaks)
}

/// Fits one line at each named position.
///
/// Each fit is a Gaussian of free centre, height and FWHM over a linear
/// background in a window of `options.fit_window` degrees around the position.
/// A fit that leaves the window is replaced by the raw maximum, so one bad
/// line does not discard the whole peak list.
///
/// # Errors
///
/// Returns [`Error::InvalidInput`] for no positions, fewer than five points in
/// a window or invalid options.
pub fn fit_peaks(
    points: &[PatternPoint],
    positions: &[f64],
    options: &PeakOptions,
) -> Result<Vec<Peak>, Error> {
    options.validate()?;
    if !options.fit_window.is_finite() || options.fit_window <= 0.0 {
        return Err(Error::InvalidInput(
            "the fit window must be positive and finite".to_owned(),
        ));
    }
    if positions.is_empty() {
        return Err(Error::InvalidInput(
            "line fitting needs at least one position".to_owned(),
        ));
    }
    if points.len() < 5 {
        return Err(Error::InvalidInput(
            "line fitting needs at least five pattern points".to_owned(),
        ));
    }
    let mut peaks = Vec::with_capacity(positions.len());
    for position in positions {
        peaks.push(fit_one(points, *position, options));
    }
    Ok(peaks)
}

fn fit_one(points: &[PatternPoint], position: f64, options: &PeakOptions) -> Peak {
    let window = options.fit_window;
    let selected: Vec<&PatternPoint> = points
        .iter()
        .filter(|point| (point.two_theta - position).abs() <= window)
        .collect();
    let fallback = raw_peak(&selected, position);
    if selected.len() < 5 {
        return fallback;
    }
    let origin = selected.first().map_or(position, |point| point.two_theta);
    let low = selected
        .iter()
        .map(|point| point.intensity)
        .fold(f64::INFINITY, f64::min);
    let high = selected
        .iter()
        .map(|point| point.intensity)
        .fold(f64::NEG_INFINITY, f64::max);
    let height = (high - low).max(1e-6);
    let initial = [height, position, window / 3.0, low, 0.0];
    let residual = |parameters: &[f64]| -> Result<Vec<f64>, Error> {
        let [amplitude, center, fwhm, offset, slope] = [
            parameters[0],
            parameters[1],
            parameters[2],
            parameters[3],
            parameters[4],
        ];
        let width = fwhm.max(1e-6);
        Ok(selected
            .iter()
            .map(|point| {
                let delta = (point.two_theta - center) / width;
                let gaussian = amplitude * (-4.0 * std::f64::consts::LN_2 * delta * delta).exp();
                let background = offset + slope * (point.two_theta - origin);
                gaussian + background - point.intensity
            })
            .collect())
    };
    let Ok(fit) = least_squares(
        &initial,
        &FitOptions {
            max_iterations: 60,
            tolerance: 1e-10,
            lambda: 1e-3,
        },
        residual,
    ) else {
        return fallback;
    };
    let [amplitude, center, fwhm, ..] = [
        fit.parameters[0],
        fit.parameters[1],
        fit.parameters[2],
        fit.parameters[3],
        fit.parameters[4],
    ];
    if !amplitude.is_finite()
        || amplitude <= 0.0
        || !fwhm.is_finite()
        || fwhm <= 1e-4
        || fwhm > 4.0 * window
        || (center - position).abs() > window
    {
        return fallback;
    }
    Peak {
        two_theta: center,
        intensity: amplitude,
        fwhm,
    }
}

fn raw_peak(points: &[&PatternPoint], position: f64) -> Peak {
    let mut best = position;
    let mut height = 0.0;
    for point in points {
        if point.intensity > height {
            height = point.intensity;
            best = point.two_theta;
        }
    }
    Peak {
        two_theta: best,
        intensity: height,
        fwhm: 0.2,
    }
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
