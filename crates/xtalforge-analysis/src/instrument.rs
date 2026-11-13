//! Radiation, geometry and instrument corrections (Task 7.2.2).
//!
//! The model turns a nominal Bragg position into the position a diffractometer
//! reports. It carries the two aberrations that move every line with angle,
//! plus the K alpha doublet and the polarisation fraction used by the profile
//! and intensity models downstream.
//!
//! # Documented formulas
//!
//! * Zero shift, in degrees: `2 theta_reported = 2 theta_nominal + zero_shift`.
//! * Sample displacement, Bragg-Brentano only, in radians:
//!
//!   ```text
//!   delta(2 theta) = -2 s cos(theta) / R
//!   ```
//!
//!   with `s` the displacement in millimetres and `R` the goniometer radius in
//!   millimetres. A positive displacement (the specimen surface above the
//!   focusing plane) therefore moves the reported line to a *lower* angle, the
//!   documented sign of the aberration. The value is converted to degrees
//!   before it is added.
//! * Debye-Scherrer (capillary) geometry has no specimen-displacement term:
//!   the reference bridge states that a capillary has no specimen
//!   displacement, whatever the flag says. Only the zero shift is applied.
//!
//! # Approximations
//!
//! * No axial divergence, specimen transparency, flat-plate illumination or
//!   receiving-slit correction is modelled here; those are profile terms of
//!   [`InstrumentProfile`] (Task 14.1.1).
//! * The K alpha doublet is placed by Bragg's law from the K alpha 1 position
//!   and its intensity ratio, not by a full two-wavelength convolution.
//! * Polarisation is stored and validated here and is propagated to the
//!   profile by [`Instrument::profile`]; the profile applies it to the
//!   Lorentz-polarisation factor.

pub mod presets;
pub mod profile;

mod width;

pub use presets::{
    INSTRUMENT_PRESET_COUNT, INSTRUMENT_PRESET_NAMES, InstrumentPreset, instrument_presets,
    preset_by_name, tch_preset_by_name,
};
pub use profile::{
    Caglioti, FcjAsymmetry, InstrumentProfile, LineShape, Monochromator, RadiationKind,
    RadiationSource, Tch, electron_wavelength,
};

use serde::{Deserialize, Serialize};
use xtalforge_core::Error;

/// Diffractometer geometry.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default, Serialize, Deserialize)]
pub enum Geometry {
    /// Flat specimen, focusing circle; specimen displacement applies.
    #[default]
    BraggBrentano,
    /// Capillary; specimen displacement does not apply.
    DebyeScherrer,
}

/// Instrument description: geometry, zero shift, displacement and radiation.
#[derive(Debug, Clone, Copy, PartialEq, Serialize, Deserialize)]
pub struct Instrument {
    /// Diffractometer geometry.
    pub geometry: Geometry,
    /// Zero-point offset in degrees two-theta, added to the nominal angle.
    pub zero_shift: f64,
    /// Specimen displacement in millimetres, positive above the focusing plane.
    pub displacement: f64,
    /// Goniometer radius in millimetres.
    pub sample_radius: f64,
    /// Intensity ratio `I(K alpha 2) / I(K alpha 1)`.
    pub ka2_ratio: f64,
    /// Polarisation fraction in `[0, 1]` used by the intensity model.
    pub polarization: f64,
}

impl Default for Instrument {
    fn default() -> Self {
        Self {
            geometry: Geometry::BraggBrentano,
            zero_shift: 0.0,
            displacement: 0.0,
            sample_radius: 240.0,
            ka2_ratio: 0.5,
            polarization: 0.5,
        }
    }
}

impl Instrument {
    /// Validates every numeric field.
    ///
    /// # Errors
    ///
    /// Returns [`Error::InvalidInput`] for a non-finite offset, a negative or
    /// non-finite radius, a radius that is not positive, a K alpha 2 ratio
    /// outside `[0, 1]`, or a polarisation fraction outside `[0, 1]`.
    pub fn validate(&self) -> Result<(), Error> {
        if !self.zero_shift.is_finite() {
            return Err(Error::InvalidInput(
                "instrument zero shift must be finite".to_owned(),
            ));
        }
        if !self.displacement.is_finite() {
            return Err(Error::InvalidInput(
                "instrument sample displacement must be finite".to_owned(),
            ));
        }
        if !self.sample_radius.is_finite() || self.sample_radius <= 0.0 {
            return Err(Error::InvalidInput(format!(
                "instrument sample radius must be positive and finite, got {}",
                self.sample_radius
            )));
        }
        if !self.ka2_ratio.is_finite() || !(0.0..=1.0).contains(&self.ka2_ratio) {
            return Err(Error::InvalidInput(format!(
                "instrument K alpha 2 ratio must be in [0, 1], got {}",
                self.ka2_ratio
            )));
        }
        if !self.polarization.is_finite() || !(0.0..=1.0).contains(&self.polarization) {
            return Err(Error::InvalidInput(format!(
                "instrument polarisation must be in [0, 1], got {}",
                self.polarization
            )));
        }
        Ok(())
    }

    /// Maps a nominal Bragg angle to the angle the instrument reports.
    ///
    /// The zero shift is added and the Bragg-Brentano displacement term
    /// `-2 s cos(theta) / R` is applied in degrees. A positive `displacement`
    /// lowers the reported angle. Non-finite input returns `NaN` rather than
    /// panicking.
    #[must_use]
    pub fn correct(&self, two_theta: f64) -> f64 {
        if !two_theta.is_finite() {
            return f64::NAN;
        }
        let theta = 0.5 * two_theta.to_radians();
        let displacement = match self.geometry {
            Geometry::BraggBrentano if self.sample_radius > 0.0 => {
                -(2.0 * self.displacement * theta.cos() / self.sample_radius)
            }
            _ => 0.0,
        };
        two_theta + self.zero_shift + displacement.to_degrees()
    }

    /// K alpha 1 / K alpha 2 pair for a K alpha 1 position.
    ///
    /// # Errors
    ///
    /// Returns [`Error::InvalidInput`] when the instrument is invalid or the
    /// position or wavelengths are not positive and finite.
    pub fn doublet(
        &self,
        two_theta: f64,
        wavelength_ka1: f64,
        wavelength_ka2: f64,
    ) -> Result<Vec<(f64, f64)>, Error> {
        self.validate()?;
        doublet_checked(two_theta, wavelength_ka1, wavelength_ka2, self.ka2_ratio)
    }

    /// Builds the profile half of the instrument model.
    ///
    /// The profile carries the radiation, a pseudo-Voigt line shape with a
    /// constant `0.1` degree FWHM as its seed, and the polarisation fraction of
    /// this instrument. Callers replace the width law, shape, `eta` or
    /// asymmetry with the values of their diffractometer.
    ///
    /// # Errors
    ///
    /// Returns [`Error::InvalidInput`] when the instrument or the radiation is
    /// invalid.
    pub fn profile(&self, radiation: RadiationSource) -> Result<InstrumentProfile, Error> {
        self.validate()?;
        radiation.validate()?;
        Ok(InstrumentProfile {
            radiation,
            monochromator: Monochromator {
                polarization: self.polarization,
            },
            ..InstrumentProfile::default()
        })
    }
}

/// K alpha 1 / K alpha 2 pair for a K alpha 1 position.
///
/// The K alpha 2 position follows from Bragg's law,
/// `sin(theta_2) = (lambda_2 / lambda_1) sin(theta_1)`, and the returned
/// weights are `1 / (1 + ratio)` and `ratio / (1 + ratio)`, so they always sum
/// to one. When the K alpha 2 line would have no Bragg solution (or an input
/// is invalid) both entries keep the input position; the K alpha 2 weight is
/// still the documented ratio. Use [`doublet_checked`] when an error is
/// preferable to a saturated value.
#[must_use]
pub fn doublet(
    two_theta: f64,
    wavelength_ka1: f64,
    wavelength_ka2: f64,
    ratio: f64,
) -> Vec<(f64, f64)> {
    let ratio = if ratio.is_finite() {
        ratio.clamp(0.0, 1.0)
    } else {
        0.0
    };
    let ka1_weight = 1.0 / (1.0 + ratio);
    let ka2_weight = ratio / (1.0 + ratio);
    let valid = two_theta.is_finite()
        && wavelength_ka1.is_finite()
        && wavelength_ka1 > 0.0
        && wavelength_ka2.is_finite()
        && wavelength_ka2 > 0.0;
    if !valid {
        return vec![(two_theta, ka1_weight), (two_theta, ka2_weight)];
    }
    let sine = (0.5 * two_theta.to_radians()).sin() * wavelength_ka2 / wavelength_ka1;
    if sine.abs() >= 1.0 {
        return vec![(two_theta, ka1_weight), (two_theta, ka2_weight)];
    }
    let two_theta_ka2 = 2.0 * sine.asin().to_degrees();
    vec![(two_theta, ka1_weight), (two_theta_ka2, ka2_weight)]
}

/// Fallible variant of [`doublet`].
///
/// # Errors
///
/// Returns [`Error::InvalidInput`] for a non-positive or non-finite position or
/// wavelength, or a ratio outside `[0, 1]`.
pub fn doublet_checked(
    two_theta: f64,
    wavelength_ka1: f64,
    wavelength_ka2: f64,
    ratio: f64,
) -> Result<Vec<(f64, f64)>, Error> {
    if !two_theta.is_finite() || two_theta <= 0.0 {
        return Err(Error::InvalidInput(format!(
            "K alpha 1 position must be positive and finite, got {two_theta}"
        )));
    }
    if !wavelength_ka1.is_finite() || wavelength_ka1 <= 0.0 {
        return Err(Error::InvalidInput(format!(
            "K alpha 1 wavelength must be positive and finite, got {wavelength_ka1}"
        )));
    }
    if !wavelength_ka2.is_finite() || wavelength_ka2 <= 0.0 {
        return Err(Error::InvalidInput(format!(
            "K alpha 2 wavelength must be positive and finite, got {wavelength_ka2}"
        )));
    }
    if !ratio.is_finite() || !(0.0..=1.0).contains(&ratio) {
        return Err(Error::InvalidInput(format!(
            "K alpha 2 ratio must be in [0, 1], got {ratio}"
        )));
    }
    Ok(doublet(two_theta, wavelength_ka1, wavelength_ka2, ratio))
}
