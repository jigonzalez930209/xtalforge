//! Profile-aware calculated powder patterns (Task 14.1.1).
//!
//! [`calculate_pxrd`](crate::calculate_pxrd) keeps the historical Gaussian
//! entry point. This module adds [`PxrdProfileOptions`] and
//! [`calculate_pxrd_profile`], which use an [`InstrumentProfile`]: any
//! radiation kind and wavelength, a constant, Caglioti or TCH width law,
//! Gaussian, Lorentzian, pseudo-Voigt or TCH line shapes, FCJ-style
//! asymmetry and the monochromator polarisation factor.
//!
//! With the profile built by [`PxrdProfileOptions::from_pxrd`] the sampled
//! pattern equals the historical one up to floating-point rounding: same
//! wavelength, Gaussian shape, constant FWHM and unmonochromated polarisation.
//!
//! # Approximation
//!
//! The sampler evaluates each line only within twelve FWHM of its centre, like
//! the historical Gaussian entry point. The far wings of a Lorentzian or TCH
//! shape are therefore not represented in the sampled points, although the
//! listed reflection intensities are unaffected. A fit that must see those
//! wings (Le Bail, Pawley, Rietveld) uses the reflection list and its own
//! profile evaluation, not the sampled points.

use serde::{Deserialize, Serialize};
use xtalforge_core::{Error, Structure};

use crate::instrument::{InstrumentProfile, RadiationSource};
use crate::pxrd::{Pattern, PatternPoint, PxrdOptions};

/// Number of samples above which the requested pattern is refused.
const MAX_PATTERN_POINTS: usize = 2_000_000;
/// How many FWHM from a peak centre the profile is evaluated.
const PROFILE_CUTOFF_FWHM: f64 = 12.0;

/// Options for [`calculate_pxrd_profile`].
#[derive(Debug, Clone, Copy, PartialEq, Serialize, Deserialize)]
pub struct PxrdProfileOptions {
    /// Instrument profile: radiation, width law, shape, asymmetry and
    /// polarisation.
    pub profile: InstrumentProfile,
    /// Inclusive `2 theta` range in degrees.
    pub two_theta_range: [f64; 2],
    /// Sampling step in degrees.
    pub step: f64,
}

impl Default for PxrdProfileOptions {
    fn default() -> Self {
        Self {
            profile: InstrumentProfile::default(),
            two_theta_range: [5.0, 80.0],
            step: 0.02,
        }
    }
}

impl PxrdProfileOptions {
    /// Builds profile options from the historical [`PxrdOptions`].
    ///
    /// The radiation keeps its wavelength (and becomes an X-ray source), the
    /// width is the constant `fwhm`, the shape is a Gaussian and the
    /// polarisation is the unmonochromated `0.5`.
    ///
    /// # Errors
    ///
    /// Returns the wavelength error of [`PxrdOptions::radiation`].
    pub fn from_pxrd(options: &PxrdOptions) -> Result<Self, Error> {
        Ok(Self {
            profile: InstrumentProfile {
                radiation: RadiationSource::from_radiation(options.radiation)?,
                shape: crate::instrument::LineShape::Gaussian,
                fwhm: options.fwhm,
                caglioti: None,
                tch: None,
                eta: 0.5,
                asymmetry: None,
                monochromator: crate::instrument::Monochromator::default(),
            },
            two_theta_range: options.two_theta_range,
            step: options.step,
        })
    }

    /// Copy with a different instrument profile.
    #[must_use]
    pub fn with_profile(mut self, profile: InstrumentProfile) -> Self {
        self.profile = profile;
        self
    }

    /// Copy with a different `2 theta` range.
    #[must_use]
    pub fn with_range(mut self, minimum: f64, maximum: f64) -> Self {
        self.two_theta_range = [minimum, maximum];
        self
    }

    /// Copy with a different sampling step.
    #[must_use]
    pub fn with_step(mut self, step: f64) -> Self {
        self.step = step;
        self
    }

    /// Validates the profile, range and step.
    ///
    /// # Errors
    ///
    /// Returns [`Error::InvalidInput`] for an invalid profile, a range outside
    /// `[0, 180]` or reversed, a non-positive step, or a pattern with more than
    /// two million points.
    pub fn validate(&self) -> Result<(), Error> {
        self.profile.validate()?;
        let [minimum, maximum] = self.two_theta_range;
        if !minimum.is_finite() || !maximum.is_finite() {
            return Err(Error::InvalidInput(
                "2 theta range must be finite".to_owned(),
            ));
        }
        if minimum < 0.0 || maximum <= minimum || maximum > 180.0 {
            return Err(Error::InvalidInput(format!(
                "2 theta range must satisfy 0 <= min < max <= 180, got {minimum}..{maximum}"
            )));
        }
        if !self.step.is_finite() || self.step <= 0.0 {
            return Err(Error::InvalidInput(format!(
                "pattern step must be positive and finite, got {}",
                self.step
            )));
        }
        let points = ((maximum - minimum) / self.step).floor() + 1.0;
        if !points.is_finite() || points > MAX_PATTERN_POINTS as f64 {
            return Err(Error::InvalidInput(format!(
                "pattern would have {points} points; use a larger step or a shorter range"
            )));
        }
        Ok(())
    }
}

/// Calculates the profile-aware powder pattern of a structure.
///
/// # Errors
///
/// Returns [`Error::InvalidInput`] for invalid options, a structure without
/// scattering atoms, a cell that would require too many reflections, or an
/// element outside the reference table (neutron radiation also needs a
/// scattering length for every element).
pub fn calculate_pxrd_profile(
    structure: &Structure,
    options: &PxrdProfileOptions,
) -> Result<Pattern, Error> {
    options.validate()?;
    let profile = &options.profile;
    let reflections = crate::pxrd::enumerate(
        structure,
        profile.radiation.wavelength,
        profile.radiation.kind,
        options.two_theta_range,
        profile.monochromator.polarization,
        true,
    )?;
    let [minimum, maximum] = options.two_theta_range;
    let count = ((maximum - minimum) / options.step).floor() as usize + 1;
    let mut points = Vec::with_capacity(count);
    for index in 0..count {
        let two_theta = minimum + index as f64 * options.step;
        let mut intensity = 0.0;
        for reflection in &reflections {
            let width = profile.fwhm_at(reflection.two_theta);
            if (two_theta - reflection.two_theta).abs() > PROFILE_CUTOFF_FWHM * width {
                continue;
            }
            intensity +=
                reflection.intensity * profile.area_normalised(reflection.two_theta, two_theta);
        }
        points.push(PatternPoint {
            two_theta,
            intensity,
        });
    }
    Ok(Pattern {
        points,
        reflections,
    })
}
