//! Profile-aware Rietveld phase model with a refinable mixing and width
//! (Tasks 14.1.1, 14.1.2 and 14.1.3).
//!
//! [`RietveldOptions`] and [`fit_phase`](super::fit_phase) are the frozen 1.x
//! phase model; they refine a scale and a constant background with a fixed
//! line shape.
//! [`RietveldProfileOptions`] adds the [`InstrumentProfile`] and refines the
//! pseudo-Voigt mixing and a global width multiplier by a coordinate search.
//! Every trial re-solves the scale and background analytically, so the profile
//! search never disturbs the linear part.
//!
//! The forward model is
//!
//! ```text
//! y(x) = scale sum_k I_k P_k(x; eta, w) + background,
//! ```
//!
//! with `I_k` the PXRD reflection intensity (multiplicity and
//! Lorentz-polarisation included), `P_k` the unit-area profile and a constant
//! background. Anisotropic displacement parameters of the structure enter
//! `I_k` through the PXRD engine.

use serde::{Deserialize, Serialize};
use xtalforge_core::{Error, Structure};
use xtalforge_io::xy::PatternPoint;

use crate::instrument::{InstrumentProfile, LineShape};
use crate::progress::RunControl;
use crate::pxrd::enumerate;

use super::{RietveldOptions, RietveldPhase};

/// Number of coarse `eta` samples before the golden-section refinement.
const ETA_GRID: usize = 5;
/// Golden-section iterations.
const REFINEMENTS: usize = 12;
/// Golden-section shrink factor.
const GOLDEN: f64 = 0.618_033_988_749_894_9;

/// Options for [`fit_phase_profile`].
#[derive(Debug, Clone, PartialEq)]
pub struct RietveldProfileOptions {
    /// The frozen phase options. Its `two_theta_range` and `background` are
    /// used; `radiation`, `step`, `fwhm` and `kind` are ignored because the
    /// [`profile`](Self::profile) supplies the radiation, width and shape and
    /// the observed grid supplies the step.
    pub base: RietveldOptions,
    /// Instrument profile used by every trial.
    pub profile: InstrumentProfile,
    /// Refine the pseudo-Voigt mixing (ignored by a TCH shape, whose mixing is
    /// fixed by its width model).
    pub refine_eta: bool,
    /// Refine a global multiplier of every profile width.
    pub refine_width: bool,
    /// Inclusive search range of `eta`.
    pub eta_range: [f64; 2],
    /// Inclusive search range of the width multiplier.
    pub width_range: [f64; 2],
    /// Cancellation and progress controls.
    pub run: RunControl,
}

impl Default for RietveldProfileOptions {
    fn default() -> Self {
        Self {
            base: RietveldOptions::default(),
            profile: InstrumentProfile::default(),
            refine_eta: true,
            refine_width: false,
            eta_range: [0.0, 1.0],
            width_range: [0.5, 2.0],
            run: RunControl::default(),
        }
    }
}

impl RietveldProfileOptions {
    /// Validates the profile, the search windows and the run controls.
    ///
    /// # Errors
    ///
    /// Returns [`Error::InvalidInput`] for an invalid profile, an `eta` or
    /// width window outside its domain, a reversed window, or a cancelled run.
    pub fn validate(&self) -> Result<(), Error> {
        self.profile.validate()?;
        validate_window(self.eta_range, 0.0, 1.0, "eta")?;
        validate_window(self.width_range, 1e-3, 100.0, "width multiplier")?;
        self.run.check()
    }
}

fn validate_window(window: [f64; 2], low: f64, high: f64, name: &str) -> Result<(), Error> {
    let [minimum, maximum] = window;
    if !minimum.is_finite() || !maximum.is_finite() {
        return Err(Error::InvalidInput(format!(
            "the {name} window is not finite"
        )));
    }
    if minimum < low || maximum > high || maximum < minimum {
        return Err(Error::InvalidInput(format!(
            "the {name} window {minimum}..{maximum} is outside {low}..{high}"
        )));
    }
    Ok(())
}

/// Result of a profile-aware Rietveld phase fit.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct RietveldProfileResult {
    /// The phase with its fitted scale, background, Rwp and pattern.
    pub phase: RietveldPhase,
    /// Mixing used by the best trial.
    pub eta: f64,
    /// Width multiplier used by the best trial.
    pub width_scale: f64,
}

/// Fits one phase scale, a constant background and the profile parameters.
///
/// The profile is area-normalised, so the returned `scale` is an
/// integrated-intensity scale and is not numerically comparable to the
/// unit-height scale of [`fit_phase`](super::fit_phase); the Rwp and the fitted
/// pattern are directly comparable.
///
/// # Errors
///
/// Returns [`Error::InvalidInput`] for fewer than three observed points,
/// non-monotone data, invalid options, a structure without reflections in the
/// observed range, or a cancelled run.
pub fn fit_phase_profile(
    structure: &Structure,
    observed: &[PatternPoint],
    options: &RietveldProfileOptions,
) -> Result<RietveldProfileResult, Error> {
    options.validate()?;
    validate_observed(observed)?;
    let first = observed.first().map_or(0.0, |point| point.two_theta);
    let last = observed.last().map_or(0.0, |point| point.two_theta);
    let range = [
        options.base.two_theta_range[0].min(first).max(0.0),
        options.base.two_theta_range[1].max(last).min(180.0),
    ];
    let reflections = enumerate(
        structure,
        options.profile.radiation.wavelength,
        options.profile.radiation.kind,
        range,
        options.profile.monochromator.polarization,
        true,
    )?;
    if reflections.is_empty() {
        return Err(Error::InvalidInput(
            "the structure has no reflections in the observed range".to_owned(),
        ));
    }
    let mut start = match options.base.background {
        crate::profiles::Background::Constant(value) => value,
        _ => 0.0,
    };
    let mut best: Option<(f64, f64, RietveldPhase)> = None;
    let mut evaluations = 0usize;
    let mut evaluate = |eta: f64, width_scale: f64| -> Result<f64, Error> {
        options.run.check()?;
        let unit = unit_curve(&reflections, observed, &options.profile, eta, width_scale);
        let (scale, background) = solve_scale_background(&unit, observed, start)?;
        start = background;
        let pattern: Vec<PatternPoint> = observed
            .iter()
            .zip(unit.iter())
            .map(|(point, value)| PatternPoint {
                two_theta: point.two_theta,
                intensity: scale * value + background,
                uncertainty: None,
            })
            .collect();
        let rwp = weighted_rwp(observed, &pattern);
        evaluations += 1;
        options
            .run
            .report("rietveld_profile", evaluations, Some(rwp));
        if best
            .as_ref()
            .is_none_or(|(_, _, current)| rwp < current.rwp)
        {
            best = Some((
                eta,
                width_scale,
                RietveldPhase {
                    scale,
                    background,
                    rwp,
                    pattern,
                },
            ));
        }
        Ok(rwp)
    };

    let mut eta = if options.profile.shape == LineShape::Tch {
        options.profile.eta_at(30.0)
    } else {
        options.profile.eta.clamp(0.0, 1.0)
    };
    let mut width_scale = 1.0;
    if options.refine_eta || options.refine_width {
        evaluate(eta, width_scale)?;
        for _ in 0..2 {
            if options.refine_eta && options.profile.shape != LineShape::Tch {
                eta = refine_scalar(options.eta_range, |value| evaluate(value, width_scale))?;
            }
            if options.refine_width {
                width_scale = refine_scalar(options.width_range, |value| evaluate(eta, value))?;
            }
        }
    } else {
        evaluate(eta, width_scale)?;
    }
    let (eta, width_scale, phase) = best.ok_or_else(|| {
        Error::InvalidInput("the Rietveld profile fit evaluated no point".to_owned())
    })?;
    Ok(RietveldProfileResult {
        phase,
        eta,
        width_scale,
    })
}

/// Coarse grid followed by a golden-section refinement inside the best
/// bracket.
fn refine_scalar<F>(window: [f64; 2], mut evaluate: F) -> Result<f64, Error>
where
    F: FnMut(f64) -> Result<f64, Error>,
{
    let [minimum, maximum] = window;
    if maximum <= minimum {
        evaluate(minimum)?;
        return Ok(minimum);
    }
    let mut best = minimum;
    let mut best_value = f64::INFINITY;
    let mut probe = |value: f64| -> Result<f64, Error> {
        let objective = evaluate(value)?;
        if objective < best_value {
            best_value = objective;
            best = value;
        }
        Ok(objective)
    };
    let mut samples = Vec::with_capacity(ETA_GRID);
    for index in 0..ETA_GRID {
        let value = minimum + (maximum - minimum) * index as f64 / (ETA_GRID - 1) as f64;
        samples.push((value, probe(value)?));
    }
    let best_index = samples
        .iter()
        .enumerate()
        .min_by(|(_, first), (_, second)| first.1.total_cmp(&second.1))
        .map(|(index, _)| index)
        .unwrap_or(0);
    let low = samples
        .get(best_index.saturating_sub(1))
        .map_or(minimum, |sample| sample.0);
    let high = samples
        .get(best_index + 1)
        .map_or(maximum, |sample| sample.0);
    let (mut low, mut high) = (low.min(high), low.max(high));
    let mut first = high - GOLDEN * (high - low);
    let mut second = low + GOLDEN * (high - low);
    let mut first_value = probe(first)?;
    let mut second_value = probe(second)?;
    for _ in 0..REFINEMENTS {
        if first_value <= second_value {
            high = second;
            second = first;
            second_value = first_value;
            first = high - GOLDEN * (high - low);
            first_value = probe(first)?;
        } else {
            low = first;
            first = second;
            first_value = second_value;
            second = low + GOLDEN * (high - low);
            second_value = probe(second)?;
        }
    }
    Ok(best)
}

/// Unit-area model curve at the observed points.
fn unit_curve(
    reflections: &[crate::pxrd::Reflection],
    observed: &[PatternPoint],
    profile: &InstrumentProfile,
    eta: f64,
    width_scale: f64,
) -> Vec<f64> {
    let mut curve = vec![0.0; observed.len()];
    for reflection in reflections {
        for (slot, point) in curve.iter_mut().zip(observed.iter()) {
            *slot += reflection.intensity
                * profile.area_normalised_scaled(
                    reflection.two_theta,
                    point.two_theta,
                    width_scale,
                    eta,
                );
        }
    }
    curve
}

/// Weighted linear least squares for the scale and a constant background,
/// with `w = 1 / max(|y|, 1)`.
fn solve_scale_background(
    unit: &[f64],
    observed: &[PatternPoint],
    start: f64,
) -> Result<(f64, f64), Error> {
    let (mut swu, mut swuu, mut sw, mut swy, mut swuy) = (0.0, 0.0, 0.0, 0.0, 0.0);
    for (point, value) in observed.iter().zip(unit.iter()) {
        let weight = 1.0 / point.intensity.abs().max(1.0);
        swu += weight * value;
        swuu += weight * value * value;
        sw += weight;
        swy += weight * point.intensity;
        swuy += weight * value * point.intensity;
    }
    if swuu <= 0.0 {
        return Err(Error::InvalidInput(
            "the calculated pattern has no intensity in the observed range".to_owned(),
        ));
    }
    let determinant = swuu * sw - swu * swu;
    if determinant.abs() > 1e-12 * swuu * sw {
        Ok((
            (swuy * sw - swu * swy) / determinant,
            (swuu * swy - swu * swuy) / determinant,
        ))
    } else {
        Ok(((swuy - start * swu) / swuu, start))
    }
}

fn weighted_rwp(observed: &[PatternPoint], pattern: &[PatternPoint]) -> f64 {
    let weighted: f64 = observed
        .iter()
        .zip(pattern.iter())
        .map(|(point, fit)| {
            let weight = 1.0 / point.intensity.abs().max(1.0);
            weight * (point.intensity - fit.intensity).powi(2)
        })
        .sum();
    let total: f64 = observed
        .iter()
        .map(|point| point.intensity.abs().max(1.0))
        .sum();
    if total > 0.0 {
        (weighted / total).sqrt()
    } else {
        0.0
    }
}

fn validate_observed(observed: &[PatternPoint]) -> Result<(), Error> {
    if observed.len() < 3 {
        return Err(Error::InvalidInput(
            "Rietveld needs at least three observed points".to_owned(),
        ));
    }
    for (index, point) in observed.iter().enumerate() {
        if !point.two_theta.is_finite() || !point.intensity.is_finite() {
            return Err(Error::InvalidInput(format!(
                "observed point {index} must be finite"
            )));
        }
        if let Some(previous) = index.checked_sub(1).and_then(|slot| observed.get(slot))
            && point.two_theta <= previous.two_theta
        {
            return Err(Error::InvalidInput(format!(
                "observed 2 theta must be strictly increasing at point {index}"
            )));
        }
    }
    Ok(())
}
