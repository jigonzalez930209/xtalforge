//! Shared line shapes and stick broadening for vibrational spectra.
//!
//! IR and Raman spectra use the same profiles and the same grid policy, so
//! the shapes live here rather than in either spectrum module.  The profile is
//! normalised to unit height at the stick position; the area under a peak is
//! therefore proportional to `intensity * fwhm` and is not renormalised.  That
//! is the convention the reference implementation uses for the simulated
//! spectra, and it keeps a print of the stick intensities meaningful.

use serde::{Deserialize, Serialize};

use crate::Error;

/// Shape used to broaden a stick spectrum.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default, Serialize, Deserialize)]
pub enum Broadening {
    /// Gaussian profile, `exp(-4 ln 2 ((x - x0) / fwhm)^2)`.
    Gaussian,
    /// Lorentzian profile, `1 / (1 + 4 ((x - x0) / fwhm)^2)`.
    #[default]
    Lorentzian,
}

fn line_shape(broadening: Broadening, offset: f64) -> f64 {
    match broadening {
        Broadening::Gaussian => (-4.0 * std::f64::consts::LN_2 * offset * offset).exp(),
        Broadening::Lorentzian => 1.0 / (1.0 + 4.0 * offset * offset),
    }
}

/// Broadens `(frequency, intensity)` sticks onto a linear grid.
///
/// The grid starts five `fwhm` below the lowest stick (clamped at zero) and
/// ends five `fwhm` above the highest one, with the requested spacing.  An
/// empty stick list gives an empty grid.  `label` names the spectrum in error
/// messages and `max_points` caps the number of returned points.
///
/// # Errors
///
/// Returns [`Error::InvalidInput`] for a non-finite stick frequency, a
/// non-finite range, or a grid larger than `max_points`.
pub(crate) fn broaden(
    sticks: &[(f64, f64)],
    fwhm: f64,
    step: f64,
    broadening: Broadening,
    max_points: usize,
    label: &str,
) -> Result<Vec<[f64; 2]>, Error> {
    if sticks.is_empty() {
        return Ok(Vec::new());
    }
    if sticks.iter().any(|(frequency, _)| !frequency.is_finite()) {
        return Err(Error::InvalidInput(format!(
            "{label} stick frequencies must be finite"
        )));
    }
    let first = sticks
        .iter()
        .map(|(frequency, _)| *frequency)
        .fold(f64::INFINITY, f64::min);
    let last = sticks
        .iter()
        .map(|(frequency, _)| *frequency)
        .fold(f64::NEG_INFINITY, f64::max);
    let margin = 5.0 * fwhm;
    let low = (first - margin).max(0.0);
    let high = last + margin;
    if !high.is_finite() || high <= low {
        return Err(Error::InvalidInput(format!(
            "{label} spectrum range is not finite"
        )));
    }
    let intervals = ((high - low) / step).ceil();
    if !intervals.is_finite() || intervals > (max_points.saturating_sub(1)) as f64 {
        return Err(Error::InvalidInput(format!(
            "{label} spectrum would need more than {max_points} points; increase the step"
        )));
    }
    let count = intervals as usize + 1;
    let mut points = Vec::with_capacity(count);
    let mut index = 0;
    while index < count {
        let x = low + step * index as f64;
        let mut intensity = 0.0;
        for (frequency, stick) in sticks {
            let offset = (x - frequency) / fwhm;
            intensity += stick * line_shape(broadening, offset);
        }
        points.push([x, intensity]);
        index += 1;
    }
    Ok(points)
}
