//! Whole-range peak refinement (Task 14.1.5).
//!
//! [`fit_peaks`](super::fit_peaks) fits each line in its own window, which is
//! fast but ignores the overlap between neighbouring lines. This module fits
//! every detected line at once: one shared FWHM, one shared mixing, one linear
//! background and a free centre and height per line, minimised by the
//! Levenberg-Marquardt engine of [`crate::fit`].
//!
//! The model at `2 theta = x` is
//!
//! ```text
//! y(x) = b0 + b1 t + sum_k H_k P(x; c_k, FWHM, eta),
//! t = (x - first) / (last - first),
//! ```
//!
//! with `P` a Gaussian, Lorentzian or pseudo-Voigt of unit height. Heights are
//! profile heights, not integrated areas, the same convention as
//! [`Peak`](super::Peak). The fit is deterministic for a given input and
//! options.

#![allow(clippy::indexing_slicing)] // parameter indices are bounded by construction

use std::cell::Cell;

use xtalforge_core::Error;
use xtalforge_io::xy::PatternPoint;

use crate::fit::{FitOptions, least_squares};
use crate::instrument::LineShape;
use crate::progress::RunControl;

use super::{Peak, PeakOptions, find_peaks};

/// Smallest shared FWHM in degrees.
const MIN_FWHM: f64 = 1e-3;
/// Largest shared FWHM in degrees.
const MAX_FWHM: f64 = 5.0;
/// Residual evaluations between progress events.
const PROGRESS_EVERY: usize = 25;

/// Options for [`refine_whole_range`].
#[derive(Debug, Clone, PartialEq)]
pub struct WholeRangeOptions {
    /// Peak detection options used to seed the fit.
    pub detection: PeakOptions,
    /// Line shape; TCH is not supported because it needs a full width model.
    pub shape: LineShape,
    /// Initial and fixed mixing when [`refine_eta`](Self::refine_eta) is off.
    pub eta: f64,
    /// Refine the shared mixing.
    pub refine_eta: bool,
    /// Largest number of lines the fit keeps (the strongest are kept).
    pub max_peaks: usize,
    /// Maximum Levenberg-Marquardt iterations.
    pub max_iterations: usize,
    /// Relative step tolerance.
    pub tolerance: f64,
}

impl Default for WholeRangeOptions {
    fn default() -> Self {
        Self {
            detection: PeakOptions::default(),
            shape: LineShape::PseudoVoigt,
            eta: 0.5,
            refine_eta: true,
            max_peaks: 60,
            max_iterations: 60,
            tolerance: 1e-8,
        }
    }
}

impl WholeRangeOptions {
    /// Validates the options.
    ///
    /// # Errors
    ///
    /// Returns [`Error::InvalidInput`] for a TCH shape, an out-of-range mixing,
    /// no peak slot, a non-positive iteration budget or tolerance, or a
    /// cancelled run.
    pub fn validate(&self) -> Result<(), Error> {
        if self.shape == LineShape::Tch {
            return Err(Error::InvalidInput(
                "whole-range peak fitting does not accept the TCH shape; use a fixed FWHM"
                    .to_owned(),
            ));
        }
        if !self.eta.is_finite() || !(0.0..=1.0).contains(&self.eta) {
            return Err(Error::InvalidInput(format!(
                "the mixing must be in [0, 1], got {}",
                self.eta
            )));
        }
        if self.max_peaks == 0 {
            return Err(Error::InvalidInput(
                "whole-range fitting needs at least one peak slot".to_owned(),
            ));
        }
        if self.max_iterations == 0 {
            return Err(Error::InvalidInput(
                "whole-range fitting needs at least one iteration".to_owned(),
            ));
        }
        if !self.tolerance.is_finite() || self.tolerance <= 0.0 {
            return Err(Error::InvalidInput(
                "the fit tolerance must be positive and finite".to_owned(),
            ));
        }
        Ok(())
    }
}

/// Result of a whole-range peak fit.
#[derive(Debug, Clone, PartialEq, serde::Serialize, serde::Deserialize)]
pub struct WholeRangeFit {
    /// Fitted lines, ordered by position.
    pub peaks: Vec<Peak>,
    /// Shared FWHM in degrees.
    pub fwhm: f64,
    /// Shared mixing used by the model.
    pub eta: f64,
    /// Background `b0 + b1 t`, with `t` normalised over the observed range.
    pub background: [f64; 2],
    /// Weighted profile R factor.
    pub rwp: f64,
    /// Weighted sum of squared residuals.
    pub chi_squared: f64,
    /// Whether the Levenberg-Marquardt step converged.
    pub converged: bool,
}

/// Refines every detected line with a shared width and a linear background.
///
/// # Errors
///
/// Returns [`Error::InvalidInput`] for fewer than five points, invalid
/// options, no detected line, or a non-finite fit.
pub fn refine_whole_range(
    points: &[PatternPoint],
    options: &WholeRangeOptions,
) -> Result<WholeRangeFit, Error> {
    refine_whole_range_with_control(points, options, &RunControl::default())
}

/// [`refine_whole_range`] with cancellation and progress controls.
///
/// The residual function checks the cancellation flag before every evaluation,
/// so a cancelled run returns quickly and reports no partial result.
///
/// # Errors
///
/// Returns [`Error::InvalidInput`] for fewer than five points, invalid
/// options, no detected line, a non-finite fit, or a cancellation.
pub fn refine_whole_range_with_control(
    points: &[PatternPoint],
    options: &WholeRangeOptions,
    control: &RunControl,
) -> Result<WholeRangeFit, Error> {
    options.validate()?;
    control.check()?;
    if points.len() < 5 {
        return Err(Error::InvalidInput(
            "whole-range fitting needs at least five pattern points".to_owned(),
        ));
    }
    let detected = find_peaks(points, &options.detection)?;
    if detected.is_empty() {
        return Err(Error::InvalidInput(
            "whole-range fitting found no peak to refine".to_owned(),
        ));
    }
    let mut peaks = detected;
    if peaks.len() > options.max_peaks {
        peaks.sort_by(|first, second| second.intensity.total_cmp(&first.intensity));
        peaks.truncate(options.max_peaks);
    }
    peaks.sort_by(|first, second| first.two_theta.total_cmp(&second.two_theta));
    if peaks.is_empty() {
        return Err(Error::InvalidInput(
            "whole-range fitting found no peak to refine".to_owned(),
        ));
    }
    let first = points.first().map_or(0.0, |point| point.two_theta);
    let last = points.last().map_or(1.0, |point| point.two_theta);
    let span = (last - first).max(1e-9);
    let fwhm = median_fwhm(&peaks).clamp(MIN_FWHM, MAX_FWHM);
    let mut initial = vec![fwhm];
    let eta_index = if options.refine_eta
        && options.shape != LineShape::Gaussian
        && options.shape != LineShape::Lorentzian
    {
        initial.push(options.eta);
        Some(1)
    } else {
        None
    };
    let background_index = initial.len();
    initial.push(
        points
            .iter()
            .map(|point| point.intensity)
            .fold(f64::INFINITY, f64::min),
    );
    initial.push(0.0);
    let peak_index = initial.len();
    for peak in &peaks {
        initial.push(peak.two_theta);
        initial.push(peak.intensity.max(1e-6));
    }
    let fixed_eta = match options.shape {
        LineShape::Gaussian => 0.0,
        LineShape::Lorentzian => 1.0,
        LineShape::PseudoVoigt => options.eta,
        LineShape::Tch => options.eta,
    };
    let evaluations = Cell::new(0usize);
    let control = control.clone();
    let weights: Vec<f64> = points.iter().map(weight).collect();
    let residual = |parameters: &[f64]| -> Result<Vec<f64>, Error> {
        control.check()?;
        let width = parameters[0].abs().clamp(MIN_FWHM, MAX_FWHM);
        let eta = eta_index.map_or(fixed_eta, |index| parameters[index].clamp(0.0, 1.0));
        let base = parameters[background_index];
        let slope = parameters[background_index + 1];
        let mut values = Vec::with_capacity(points.len());
        let mut chi_squared = 0.0;
        for (row, point) in points.iter().enumerate() {
            let t = (point.two_theta - first) / span;
            let mut model = base + slope * t;
            let (pairs, _) = parameters[peak_index..].as_chunks::<2>();
            for pair in pairs {
                let centre = pair[0];
                let height = pair[1];
                let delta = (point.two_theta - centre) / width;
                model += height * shape(options.shape, eta, delta);
            }
            let value = weights.get(row).copied().unwrap_or(1.0).sqrt() * (model - point.intensity);
            chi_squared += value * value;
            values.push(value);
        }
        let done = evaluations.get() + 1;
        evaluations.set(done);
        if done == 1 || done.is_multiple_of(PROGRESS_EVERY) {
            control.report("whole_range", done, Some(chi_squared));
        }
        Ok(values)
    };
    let fit = least_squares(
        &initial,
        &FitOptions {
            max_iterations: options.max_iterations,
            tolerance: options.tolerance,
            lambda: 1e-3,
        },
        residual,
    )?;
    if !fit.parameters.iter().all(|value| value.is_finite()) {
        return Err(Error::InvalidInput(
            "the whole-range fit produced a non-finite parameter".to_owned(),
        ));
    }
    let width = fit.parameters[0].abs().clamp(MIN_FWHM, MAX_FWHM);
    let eta = eta_index.map_or(fixed_eta, |index| fit.parameters[index].clamp(0.0, 1.0));
    let base = fit.parameters[background_index];
    let slope = fit.parameters[background_index + 1];
    let mut fitted = Vec::with_capacity(peaks.len());
    let (pairs, _) = fit.parameters[peak_index..].as_chunks::<2>();
    for pair in pairs {
        let centre = pair[0];
        let height = pair[1].max(0.0);
        fitted.push(Peak {
            two_theta: centre,
            intensity: height,
            fwhm: width,
        });
    }
    fitted.sort_by(|first, second| first.two_theta.total_cmp(&second.two_theta));
    let weighted_observed: f64 = points
        .iter()
        .zip(weights.iter())
        .map(|(point, weight)| weight * point.intensity * point.intensity)
        .sum();
    let rwp = if weighted_observed > 0.0 {
        (fit.chi_squared / weighted_observed).sqrt()
    } else {
        0.0
    };
    Ok(WholeRangeFit {
        peaks: fitted,
        fwhm: width,
        eta,
        background: [base, slope],
        rwp,
        chi_squared: fit.chi_squared,
        converged: fit.converged,
    })
}

/// Unit-height line shape of one profile.
fn shape(shape: LineShape, eta: f64, delta: f64) -> f64 {
    let gaussian = (-4.0 * std::f64::consts::LN_2 * delta * delta).exp();
    let lorentzian = 1.0 / (1.0 + 4.0 * delta * delta);
    match shape {
        LineShape::Gaussian => gaussian,
        LineShape::Lorentzian => lorentzian,
        LineShape::PseudoVoigt | LineShape::Tch => {
            let eta = eta.clamp(0.0, 1.0);
            eta * lorentzian + (1.0 - eta) * gaussian
        }
    }
}

/// Median FWHM of the detected lines, clamped to a sane window.
fn median_fwhm(peaks: &[Peak]) -> f64 {
    let mut widths: Vec<f64> = peaks
        .iter()
        .map(|peak| peak.fwhm)
        .filter(|value| value.is_finite() && *value > 0.0)
        .collect();
    widths.sort_by(f64::total_cmp);
    widths
        .get(widths.len() / 2)
        .copied()
        .unwrap_or(0.1)
        .clamp(0.01, 1.0)
}

fn weight(point: &PatternPoint) -> f64 {
    match point.uncertainty {
        Some(value) if value.is_finite() && value > 0.0 => 1.0 / (value * value),
        _ => 1.0 / point.intensity.abs().max(1.0),
    }
}
