//! Refinable pseudo-Voigt mixing in the Le Bail extraction (Task 14.1.3).
//!
//! [`LebailOptions`] keeps its frozen `eta`; [`LebailProfileOptions`] wraps it
//! and refines the mixing with an outer golden-section search. Each trial is a
//! full Le Bail iteration, so the partitioned intensities are recomputed at
//! every `eta`.
//!
//! The search minimises the profile R factor `Rp` of the extraction, not the
//! weighted `Rwp`: the Le Bail partition is an unweighted apportionment of the
//! observed counts, and the weighted factor is dominated by the background
//! partition at the ends of the range. `Rwp` is still reported in the result.

use serde::{Deserialize, Serialize};
use xtalforge_core::{Error, Structure};
use xtalforge_io::xy::PatternPoint;

use crate::progress::RunControl;

use super::{LebailOptions, LebailResult, le_bail};

/// Number of coarse `eta` samples before the golden-section refinement.
const ETA_GRID: usize = 5;
/// Golden-section iterations after the coarse grid.
const ETA_REFINEMENTS: usize = 14;
/// Golden-section shrink factor.
const GOLDEN: f64 = 0.618_033_988_749_894_9;

/// Options for [`le_bail_profile`].
#[derive(Debug, Clone, PartialEq)]
pub struct LebailProfileOptions {
    /// The frozen Le Bail options, used as the base of every trial.
    pub base: LebailOptions,
    /// Refine the pseudo-Voigt mixing `eta`.
    pub refine_eta: bool,
    /// Inclusive search range of `eta`.
    pub eta_range: [f64; 2],
    /// Cancellation and progress controls.
    pub run: RunControl,
}

impl Default for LebailProfileOptions {
    fn default() -> Self {
        Self {
            base: LebailOptions::default(),
            refine_eta: true,
            eta_range: [0.0, 1.0],
            run: RunControl::default(),
        }
    }
}

impl LebailProfileOptions {
    /// Validates the `eta` window and the run controls.
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

/// Result of a profile-refined Le Bail extraction.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct LebailProfileResult {
    /// The underlying Le Bail extraction at the best profile.
    pub fit: LebailResult,
    /// Refined pseudo-Voigt mixing.
    pub eta: f64,
}

/// Extracts reflection intensities with a refinable pseudo-Voigt mixing.
///
/// # Errors
///
/// Returns the errors of [`le_bail`], invalid options, or a cancellation from
/// the attached [`RunControl`].
pub fn le_bail_profile(
    structure: &Structure,
    observed: &[PatternPoint],
    options: &LebailProfileOptions,
) -> Result<LebailProfileResult, Error> {
    options.validate()?;
    let inner_options = |eta: f64| {
        let mut inner = options.base;
        inner.eta = eta;
        inner
    };
    let mut best: Option<(f64, LebailResult)> = None;
    let mut evaluations = 0usize;
    let mut evaluate = |eta: f64| -> Result<f64, Error> {
        options.run.check()?;
        let fit = le_bail(structure, observed, &inner_options(eta))?;
        evaluations += 1;
        options
            .run
            .report("lebail_eta", evaluations, Some(fit.r_factor));
        let value = fit.r_factor;
        if best.as_ref().is_none_or(|(best_eta, best_fit)| {
            fit.r_factor < best_fit.r_factor
                || (fit.r_factor == best_fit.r_factor && eta < *best_eta)
        }) {
            best = Some((eta, fit));
        }
        Ok(value)
    };
    if !options.refine_eta {
        evaluate(options.base.eta.clamp(0.0, 1.0))?;
    } else {
        let [minimum, maximum] = options.eta_range;
        let mut samples = Vec::with_capacity(ETA_GRID);
        if maximum > minimum {
            for index in 0..ETA_GRID {
                let eta = minimum + (maximum - minimum) * index as f64 / (ETA_GRID - 1) as f64;
                samples.push((eta, evaluate(eta)?));
            }
        } else {
            samples.push((minimum, evaluate(minimum)?));
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
    let (eta, fit) = best.ok_or_else(|| {
        Error::InvalidInput("the Le Bail profile fit evaluated no eta".to_owned())
    })?;
    Ok(LebailProfileResult { fit, eta })
}
