//! Refinable pseudo-Voigt mixing and profile width in the Pawley fit
//! (Task 14.1.3).
//!
//! [`PawleyOptions`] is the frozen 1.x option struct; it keeps its fixed
//! `eta`. This module adds [`PawleyProfileOptions`], which wraps the old
//! options and refines `eta` (and, on request, the constant FWHM through the
//! existing `refine_size` step) with an outer golden-section search. Every
//! trial is a complete Pawley fit, so the linear intensities and background
//! are re-solved at each `eta` and units stay on the documented scale.
//!
//! With [`PawleyProfileOptions::refine_eta`] disabled this entry point returns
//! exactly the old result, which is how callers migrate one flag at a time.

use serde::{Deserialize, Serialize};
use xtalforge_core::Error;
use xtalforge_io::xy::PatternPoint;

use crate::progress::RunControl;

use super::{PawleyOptions, PawleyResult, pawley_refine};

/// Number of coarse `eta` samples before the golden-section refinement.
const ETA_GRID: usize = 5;
/// Golden-section iterations after the coarse grid.
const ETA_REFINEMENTS: usize = 16;
/// Golden-section shrink factor.
const GOLDEN: f64 = 0.618_033_988_749_894_9;

/// Options for [`pawley_refine_profile`].
#[derive(Debug, Clone, PartialEq)]
pub struct PawleyProfileOptions {
    /// The frozen Pawley options, used as the base of every trial.
    pub base: PawleyOptions,
    /// Refine the pseudo-Voigt mixing `eta`.
    pub refine_eta: bool,
    /// Refine the constant FWHM through the existing `refine_size` step.
    pub refine_width: bool,
    /// Inclusive search range of `eta`.
    pub eta_range: [f64; 2],
    /// Cancellation and progress controls.
    pub run: RunControl,
}

impl Default for PawleyProfileOptions {
    fn default() -> Self {
        Self {
            base: PawleyOptions::default(),
            refine_eta: true,
            refine_width: false,
            eta_range: [0.0, 1.0],
            run: RunControl::default(),
        }
    }
}

impl PawleyProfileOptions {
    /// Validates `eta` range and the run controls.
    ///
    /// # Errors
    ///
    /// Returns [`Error::InvalidInput`] for a non-finite or out-of-range `eta`
    /// window, a reversed window, or a cancelled run.
    pub fn validate(&self) -> Result<(), Error> {
        let [minimum, maximum] = self.eta_range;
        if !minimum.is_finite() || !maximum.is_finite() {
            return Err(Error::InvalidInput(
                "the eta range must be finite".to_owned(),
            ));
        }
        if !(0.0..=1.0).contains(&minimum) || !(0.0..=1.0).contains(&maximum) {
            return Err(Error::InvalidInput(format!(
                "the eta range must lie in [0, 1], got {minimum}..{maximum}"
            )));
        }
        if maximum < minimum {
            return Err(Error::InvalidInput(format!(
                "the eta range is reversed: {minimum}..{maximum}"
            )));
        }
        self.run.check()
    }
}

/// Result of a profile-refined Pawley fit.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct PawleyProfileResult {
    /// The underlying Pawley fit at the best profile.
    pub fit: PawleyResult,
    /// Refined pseudo-Voigt mixing.
    pub eta: f64,
    /// Refined constant FWHM in degrees.
    pub fwhm: f64,
}

/// Fits a Pawley model with a refinable pseudo-Voigt mixing.
///
/// # Errors
///
/// Returns the errors of [`pawley_refine`], invalid options, or a cancellation
/// from the attached [`RunControl`].
pub fn pawley_refine_profile(
    structure: &xtalforge_core::Structure,
    observed: &[PatternPoint],
    options: &PawleyProfileOptions,
) -> Result<PawleyProfileResult, Error> {
    options.validate()?;
    let inner_options = |eta: f64| {
        let mut inner = options.base;
        inner.eta = eta;
        if options.refine_width {
            inner.refine_size = true;
        }
        inner
    };
    let mut best: Option<(f64, PawleyResult)> = None;
    let mut evaluations = 0usize;
    let mut evaluate = |eta: f64| -> Result<f64, Error> {
        options.run.check()?;
        let fit = pawley_refine(structure, observed, &inner_options(eta))?;
        evaluations += 1;
        options.run.report("pawley_eta", evaluations, Some(fit.rwp));
        let value = fit.rwp;
        if best.as_ref().is_none_or(|(best_eta, best_fit)| {
            fit.rwp < best_fit.rwp || (fit.rwp == best_fit.rwp && eta < *best_eta)
        }) {
            best = Some((eta, fit));
        }
        Ok(value)
    };

    if !options.refine_eta {
        let eta = options.base.eta.clamp(0.0, 1.0);
        evaluate(eta)?;
    } else {
        let [minimum, maximum] = options.eta_range;
        let mut samples = Vec::with_capacity(ETA_GRID);
        if maximum > minimum {
            for index in 0..ETA_GRID {
                let fraction = index as f64 / (ETA_GRID - 1) as f64;
                samples.push((
                    minimum + fraction * (maximum - minimum),
                    evaluate(minimum + fraction * (maximum - minimum))?,
                ));
            }
        } else {
            samples.push((minimum, evaluate(minimum)?));
        }
        let best_index = samples
            .iter()
            .enumerate()
            .min_by(|first, second| first.1.1.total_cmp(&second.1.1))
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
        let mut first_value = evaluate(first)?;
        let mut second_value = evaluate(second)?;
        for _ in 0..ETA_REFINEMENTS {
            options.run.check()?;
            if first_value <= second_value {
                high = second;
                second = first;
                second_value = first_value;
                first = high - GOLDEN * (high - low);
                first_value = evaluate(first)?;
            } else {
                low = first;
                first = second;
                first_value = second_value;
                second = low + GOLDEN * (high - low);
                second_value = evaluate(second)?;
            }
        }
    }
    let (eta, fit) = best
        .ok_or_else(|| Error::InvalidInput("the Pawley profile fit evaluated no eta".to_owned()))?;
    let fwhm = if options.refine_width || options.base.refine_size {
        fit.size
    } else {
        options.base.fwhm
    };
    Ok(PawleyProfileResult { fit, eta, fwhm })
}
