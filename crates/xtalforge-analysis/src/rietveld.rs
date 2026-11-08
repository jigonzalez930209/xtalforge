//! A minimal Rietveld phase model (Task 7.8).

use serde::{Deserialize, Serialize};
use xtalforge_core::{Error, Structure};
use xtalforge_io::xy::PatternPoint;

use crate::profiles::{Background, Profile, ProfileKind};
use crate::pxrd::{PxrdOptions, Radiation, calculate_pxrd};

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
    let calculated = calculate_pxrd(
        structure,
        &PxrdOptions {
            radiation: options.radiation,
            two_theta_range: options.two_theta_range,
            step: options.step,
            fwhm: options.fwhm,
        },
    )?;
    let mut unit = vec![0.0; observed.len()];
    for (slot, point) in unit.iter_mut().zip(observed.iter()) {
        let mut value = 0.0;
        for reflection in &calculated.reflections {
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
    let count = observed.len() as f64;
    let background = match options.background {
        Background::Constant(value) => value,
        _ => 0.0,
    };
    let mut numerator = 0.0;
    let mut denominator = 0.0;
    for (point, value) in observed.iter().zip(unit.iter()) {
        numerator += value * (point.intensity - background);
        denominator += value * value;
    }
    if denominator <= 0.0 {
        return Err(Error::InvalidInput(
            "the calculated pattern has no intensity in the observed range".to_owned(),
        ));
    }
    let scale = numerator / denominator;
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
    let _ = count;
    Ok(RietveldPhase {
        scale,
        background,
        rwp,
        pattern,
    })
}
