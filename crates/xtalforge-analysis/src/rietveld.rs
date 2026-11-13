//! A minimal Rietveld phase model (Task 7.8).

pub mod profile;

pub use profile::{RietveldProfileOptions, RietveldProfileResult, fit_phase_profile};

use serde::{Deserialize, Serialize};
use xtalforge_core::{Error, Structure};
use xtalforge_io::xy::PatternPoint;

use crate::instrument::RadiationKind;
use crate::profiles::{Background, Profile, ProfileKind};
use crate::pxrd::{Radiation, enumerate};

/// Options for the phase simulation and scale fit.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct RietveldOptions {
    /// Radiation.
    pub radiation: Radiation,
    /// Range in degrees.
    pub two_theta_range: [f64; 2],
    /// Step in degrees.
    pub step: f64,
    /// Profile FWHM in degrees.
    pub fwhm: f64,
    /// Background model.
    pub background: Background,
    /// Profile shape.
    pub kind: ProfileKind,
}

impl Default for RietveldOptions {
    fn default() -> Self {
        Self {
            radiation: Radiation::CuKa,
            two_theta_range: [5.0, 80.0],
            step: 0.02,
            fwhm: 0.1,
            background: Background::Constant(0.0),
            kind: ProfileKind::PseudoVoigt,
        }
    }
}

/// A simulated phase and its refined scale.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct RietveldPhase {
    /// Scale applied to the calculated intensities.
    pub scale: f64,
    /// Background level.
    pub background: f64,
    /// Weighted profile R factor.
    pub rwp: f64,
    /// Calculated points.
    pub pattern: Vec<PatternPoint>,
}

/// Fits one phase scale and a constant background to an observed pattern.
///
/// This is the structural-model hook: the reflection positions and relative
/// intensities come from `Structure` through the PXRD engine, and only the
/// scale and background are refined here. Atomic positions and ADPs are the
/// caller's model and can be changed between runs.
///
/// # Errors
///
/// Returns [`Error::InvalidInput`] for no observed points or a degenerate
/// fit.
pub fn fit_phase(
    structure: &Structure,
    observed: &[PatternPoint],
    options: &RietveldOptions,
) -> Result<RietveldPhase, Error> {
    if observed.len() < 3 {
        return Err(Error::InvalidInput(
            "Rietveld needs at least three observed points".to_owned(),
        ));
    }
    let calculated = enumerate(
        structure,
        options.radiation.wavelength()?,
        RadiationKind::XRay,
        options.two_theta_range,
        0.5,
        true,
    )?;
    let mut unit = vec![0.0; observed.len()];
    for (slot, point) in unit.iter_mut().zip(observed.iter()) {
        let mut value = 0.0;
        for reflection in &calculated {
            value += Profile {
                kind: options.kind,
                center: reflection.two_theta,
                fwhm: options.fwhm,
                eta: 0.5,
                height: reflection.intensity,
            }
            .value(point.two_theta);
        }
        *slot = value;
    }
    // Weighted linear least squares for the scale and a constant background,
    // w = 1 / max(y_obs, 1): minimise sum w (y_obs - s u - b)^2. A
    // `Background::Constant` start value is the fallback when the system is
    // singular (a calculated pattern that is flat across the range).
    let start = match options.background {
        Background::Constant(value) => value,
        _ => 0.0,
    };
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
    let (scale, background) = if determinant.abs() > 1e-12 * swuu * sw {
        (
            (swuy * sw - swu * swy) / determinant,
            (swuu * swy - swu * swuy) / determinant,
        )
    } else {
        ((swuy - start * swu) / swuu, start)
    };
    let pattern: Vec<PatternPoint> = observed
        .iter()
        .zip(unit.iter())
        .map(|(point, value)| PatternPoint {
            two_theta: point.two_theta,
            intensity: scale * value + background,
            uncertainty: None,
        })
        .collect();
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
    let rwp = (weighted / total).sqrt();
    Ok(RietveldPhase {
        scale,
        background,
        rwp,
        pattern,
    })
}
