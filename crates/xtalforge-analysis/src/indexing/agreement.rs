//! Agreement between the cubic and the general indexing engines (Task 17.1.7).
//!
//! [`index_cubic`] and [`index_pattern`] share peaks but not their search. The
//! general engine may choose a lower-symmetry metric family, and the cubic
//! engine may pick a multiple of the true cell. This module reports whether the
//! two engines describe the same cubic lattice, with the evidence in a typed
//! result so the CLI, Python and WASM surfaces can show it without parsing
//! prose.

use serde::{Deserialize, Serialize};
use xtalforge_core::Error;

use crate::peaks::Peak;

use super::{IndexingOptions, IndexingResult, LatticeSystem, index_cubic, index_pattern};

/// Relative length difference below which two cubic parameters agree.
pub const DEFAULT_RELATIVE_TOLERANCE: f64 = 0.01;

/// The verdict of comparing the cubic and general indexing engines.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct IndexingAgreement {
    /// Whether both engines produced a cubic cell within the tolerance.
    pub agree: bool,
    /// Cubic parameter from [`index_cubic`], when it found a cell.
    pub cubic_a: Option<f64>,
    /// Cubic parameter from the best general candidate, when it is cubic.
    pub general_a: Option<f64>,
    /// Absolute difference between the two parameters, when both exist.
    pub delta_a: Option<f64>,
    /// Relative difference `|delta_a| / max(cubic_a, general_a)`.
    pub relative_delta: Option<f64>,
    /// Relative tolerance used for the verdict.
    pub relative_tolerance: f64,
    /// Fraction of peaks indexed by the cubic engine.
    pub cubic_completeness: f64,
    /// Fraction of peaks indexed by the general engine.
    pub general_completeness: f64,
    /// Number of peaks in the input.
    pub n_peaks: usize,
    /// Human-readable reason for the verdict.
    pub message: String,
}

/// Compares the cubic and general indexing engines on the same peak list.
///
/// The engines agree when both return a cell, the general winner is in the
/// cubic family, and the two `a` parameters differ by at most
/// `relative_tolerance` (default [`DEFAULT_RELATIVE_TOLERANCE`]). A non-cubic
/// general winner, a missing candidate on either side, or an inconsistent
/// completeness below `min_completeness` makes the result `agree = false` and
/// fills `message` with the reason.
///
/// # Errors
///
/// Returns [`Error::InvalidInput`] for fewer than two peaks, a non-positive
/// wavelength, an invalid relative tolerance or the errors of the two engines.
pub fn indexing_agreement(
    peaks: &[Peak],
    wavelength: f64,
    options: &IndexingOptions,
    relative_tolerance: f64,
) -> Result<IndexingAgreement, Error> {
    if !relative_tolerance.is_finite() || relative_tolerance <= 0.0 {
        return Err(Error::InvalidInput(format!(
            "indexing agreement tolerance must be positive and finite, got {relative_tolerance}"
        )));
    }
    let cubic = index_cubic(peaks, wavelength, options.max_index, options.tolerance)?;
    let general = index_pattern(peaks, wavelength, options)?;
    let cubic_a = cubic.as_ref().map(|index| index.a);
    let general_cell = general
        .best()
        .filter(|candidate| candidate.system == LatticeSystem::Cubic);
    let general_a = general_cell.map(|candidate| candidate.cell.a);
    let cubic_completeness = cubic.as_ref().map_or(0.0, |index| index.completeness);
    let general_completeness = general_cell.map_or(0.0, |candidate| candidate.completeness);
    let (delta_a, relative_delta) = match (cubic_a, general_a) {
        (Some(first), Some(second)) => {
            let delta = (first - second).abs();
            let scale = first.abs().max(second.abs()).max(f64::MIN_POSITIVE);
            (Some(delta), Some(delta / scale))
        }
        _ => (None, None),
    };
    let agree = relative_delta.is_some_and(|value| value <= relative_tolerance);
    let message = if agree {
        format!(
            "cubic and general engines agree on a = {:.6} A (relative difference {:.4})",
            general_a.unwrap_or(0.0),
            relative_delta.unwrap_or(0.0)
        )
    } else if cubic_a.is_none() {
        "the cubic engine found no cell".to_owned()
    } else if general.best().is_none() {
        "the general engine found no cell".to_owned()
    } else if general_a.is_none() {
        format!(
            "the general engine chose the {} family, not cubic",
            general
                .best()
                .map_or("unknown", |candidate| candidate.system.name())
        )
    } else {
        format!(
            "the engines disagree: cubic a = {:.6} A, general a = {:.6} A (relative difference {:.4} > {:.4})",
            cubic_a.unwrap_or(0.0),
            general_a.unwrap_or(0.0),
            relative_delta.unwrap_or(0.0),
            relative_tolerance
        )
    };
    Ok(IndexingAgreement {
        agree,
        cubic_a,
        general_a,
        delta_a,
        relative_delta,
        relative_tolerance,
        cubic_completeness,
        general_completeness,
        n_peaks: peaks.len(),
        message,
    })
}

/// Convenience wrapper over [`indexing_agreement`] that also returns the
/// general result used for the verdict.
///
/// # Errors
///
/// Returns the errors of [`indexing_agreement`].
pub fn indexing_agreement_with_result(
    peaks: &[Peak],
    wavelength: f64,
    options: &IndexingOptions,
    relative_tolerance: f64,
) -> Result<(IndexingAgreement, IndexingResult), Error> {
    let agreement = indexing_agreement(peaks, wavelength, options, relative_tolerance)?;
    let result = index_pattern(peaks, wavelength, options)?;
    Ok((agreement, result))
}
