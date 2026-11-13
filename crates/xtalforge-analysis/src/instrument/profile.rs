//! Instrument profile model: radiation, width law, line shape, asymmetry and
//! polarisation (Task 14.1.1).
//!
//! This module is the profile half of a diffractometer model. The
//! [`Instrument`](super::Instrument) keeps the geometric aberrations (zero
//! shift, specimen displacement); `InstrumentProfile` keeps everything that
//! changes the *shape* of a line:
//!
//! * the radiation, with its kind and wavelength;
//! * the width law, either a constant FWHM, the Caglioti law or the
//!   Thompson-Cox-Hastings (TCH) law of the width module;
//! * the line shape: Gaussian, Lorentzian, pseudo-Voigt or TCH;
//! * an optional Finger-Cox-Jephcoat (FCJ) style low-angle asymmetry;
//! * the monochromator polarisation fraction.
//!
//! The width laws, line shapes and asymmetry formulas live in the private
//! width module; this module adds radiation and polarisation.
//!
//! # Polarisation
//!
//! An X-ray beam with polarised fraction `p` uses
//! `LP = 2 (1 - p + p cos^2 2theta) / (sin^2 theta cos theta)`; at `p = 0.5`
//! this is the traditional `(1 + cos^2 2theta)` factor divided by
//! `sin^2 theta cos theta`. Neutron and electron profiles use the Lorentz
//! factor alone.

use serde::{Deserialize, Serialize};
use xtalforge_core::Error;

pub use super::width::{Caglioti, FcjAsymmetry, LineShape, Tch};
use super::width::{MIN_FWHM, shape_value};

/// Kind of probe radiation.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default, Serialize, Deserialize)]
pub enum RadiationKind {
    /// X-rays; atomic form factors and the polarisation factor apply.
    #[default]
    XRay,
    /// Neutrons; coherent scattering lengths are used instead of form factors.
    Neutron,
    /// Electrons; the Mott-Bethe electron scattering factor applies.
    Electron,
}

/// Radiation with its kind and wavelength in Angstrom.
#[derive(Debug, Clone, Copy, PartialEq, Serialize, Deserialize)]
pub struct RadiationSource {
    /// Probe kind.
    pub kind: RadiationKind,
    /// Wavelength in Angstrom.
    pub wavelength: f64,
}

impl Default for RadiationSource {
    fn default() -> Self {
        Self {
            kind: RadiationKind::XRay,
            wavelength: 1.541_84,
        }
    }
}

impl RadiationSource {
    /// X-ray source at a wavelength in Angstrom.
    ///
    /// # Errors
    ///
    /// Returns [`Error::InvalidInput`] for a non-positive or non-finite
    /// wavelength.
    pub fn xray(wavelength: f64) -> Result<Self, Error> {
        Self::checked(RadiationKind::XRay, wavelength)
    }

    /// Neutron source at a wavelength in Angstrom.
    ///
    /// # Errors
    ///
    /// Returns [`Error::InvalidInput`] for a non-positive or non-finite
    /// wavelength.
    pub fn neutron(wavelength: f64) -> Result<Self, Error> {
        Self::checked(RadiationKind::Neutron, wavelength)
    }

    /// Electron source from an accelerating voltage in kilovolt.
    ///
    /// # Errors
    ///
    /// Returns [`Error::InvalidInput`] for a non-positive or non-finite
    /// voltage.
    pub fn electron(voltage_kv: f64) -> Result<Self, Error> {
        Ok(Self {
            kind: RadiationKind::Electron,
            wavelength: electron_wavelength(voltage_kv)?,
        })
    }

    /// Converts the named laboratory radiation of the PXRD engine.
    ///
    /// # Errors
    ///
    /// Returns the error of [`Radiation::wavelength`](crate::pxrd::Radiation).
    pub fn from_radiation(radiation: crate::pxrd::Radiation) -> Result<Self, Error> {
        Self::xray(radiation.wavelength()?)
    }

    fn checked(kind: RadiationKind, wavelength: f64) -> Result<Self, Error> {
        let source = Self { kind, wavelength };
        source.validate()?;
        Ok(source)
    }

    /// Validates the wavelength.
    ///
    /// # Errors
    ///
    /// Returns [`Error::InvalidInput`] for a non-positive or non-finite
    /// wavelength.
    pub fn validate(&self) -> Result<(), Error> {
        if !self.wavelength.is_finite() || self.wavelength <= 0.0 {
            return Err(Error::InvalidInput(format!(
                "radiation wavelength must be positive and finite, got {}",
                self.wavelength
            )));
        }
        Ok(())
    }
}

/// Relativistic electron wavelength in Angstrom for an accelerating voltage in
/// kilovolt.
///
/// The formula is `lambda = 1.22639 / sqrt(V (1 + 0.9785e-6 V))` nanometre
/// with `V` in volt; the result is converted to Angstrom. It reproduces the
/// tabulated 200 kV value of `0.02508 Angstrom` within `2e-4`.
///
/// # Errors
///
/// Returns [`Error::InvalidInput`] for a non-positive or non-finite voltage.
pub fn electron_wavelength(voltage_kv: f64) -> Result<f64, Error> {
    if !voltage_kv.is_finite() || voltage_kv <= 0.0 {
        return Err(Error::InvalidInput(format!(
            "electron accelerating voltage must be positive and finite, got {voltage_kv}"
        )));
    }
    let volts = voltage_kv * 1000.0;
    let wavelength_nm = 1.226_39 / (volts * (1.0 + 9.785e-7 * volts)).sqrt();
    let wavelength = wavelength_nm * 10.0;
    if !wavelength.is_finite() || wavelength <= 0.0 {
        return Err(Error::InvalidInput(
            "the electron wavelength is not finite".to_owned(),
        ));
    }
    Ok(wavelength)
}

/// Monochromator polarisation fraction.
#[derive(Debug, Clone, Copy, PartialEq, Serialize, Deserialize)]
pub struct Monochromator {
    /// Fraction of the beam polarised in the scattering plane, in `[0, 1]`;
    /// `0.5` is an unmonochromated beam.
    pub polarization: f64,
}

impl Default for Monochromator {
    fn default() -> Self {
        Self { polarization: 0.5 }
    }
}

impl Monochromator {
    /// Validates the fraction.
    ///
    /// # Errors
    ///
    /// Returns [`Error::InvalidInput`] for a non-finite or out-of-range
    /// fraction.
    pub fn validate(&self) -> Result<(), Error> {
        if !self.polarization.is_finite() || !(0.0..=1.0).contains(&self.polarization) {
            return Err(Error::InvalidInput(format!(
                "monochromator polarisation must be in [0, 1], got {}",
                self.polarization
            )));
        }
        Ok(())
    }
}

/// Instrument profile: radiation, width law, line shape, asymmetry and
/// polarisation.
#[derive(Debug, Clone, Copy, PartialEq, Serialize, Deserialize)]
pub struct InstrumentProfile {
    /// Radiation source.
    pub radiation: RadiationSource,
    /// Line shape.
    pub shape: LineShape,
    /// Constant FWHM in degrees; used when neither `caglioti` nor `tch` is set.
    pub fwhm: f64,
    /// Caglioti width law; overrides `fwhm` when present.
    pub caglioti: Option<Caglioti>,
    /// TCH width law; required when `shape` is [`LineShape::Tch`].
    pub tch: Option<Tch>,
    /// Pseudo-Voigt mixing in `[0, 1]` (0 = Gaussian, 1 = Lorentzian).
    pub eta: f64,
    /// Low-angle asymmetry, when present.
    pub asymmetry: Option<FcjAsymmetry>,
    /// Monochromator polarisation.
    pub monochromator: Monochromator,
}

impl Default for InstrumentProfile {
    fn default() -> Self {
        Self {
            radiation: RadiationSource::default(),
            shape: LineShape::PseudoVoigt,
            fwhm: 0.1,
            caglioti: None,
            tch: None,
            eta: 0.5,
            asymmetry: None,
            monochromator: Monochromator::default(),
        }
    }
}

impl InstrumentProfile {
    /// Validates every field.
    ///
    /// # Errors
    ///
    /// Returns [`Error::InvalidInput`] for an invalid radiation, a non-positive
    /// FWHM, an out-of-range `eta`, an invalid width law, asymmetry,
    /// polarisation fraction, or a TCH shape without a TCH model.
    pub fn validate(&self) -> Result<(), Error> {
        self.radiation.validate()?;
        if !self.fwhm.is_finite() || self.fwhm <= 0.0 {
            return Err(Error::InvalidInput(format!(
                "profile FWHM must be positive and finite, got {}",
                self.fwhm
            )));
        }
        if !self.eta.is_finite() || !(0.0..=1.0).contains(&self.eta) {
            return Err(Error::InvalidInput(format!(
                "profile mixing must be in [0, 1], got {}",
                self.eta
            )));
        }
        if let Some(caglioti) = self.caglioti {
            caglioti.validate()?;
        }
        if let Some(tch) = self.tch {
            tch.validate()?;
        }
        if self.shape == LineShape::Tch && self.tch.is_none() {
            return Err(Error::InvalidInput(
                "the TCH line shape needs a TCH width model".to_owned(),
            ));
        }
        if let Some(asymmetry) = self.asymmetry {
            asymmetry.validate()?;
        }
        self.monochromator.validate()
    }

    /// FWHM in degrees at a Bragg angle.
    ///
    /// A TCH shape uses the TCH model, any other shape uses the Caglioti law
    /// when it is set and the constant `fwhm` otherwise.
    #[must_use]
    pub fn fwhm_at(&self, two_theta: f64) -> f64 {
        if self.shape == LineShape::Tch
            && let Some(tch) = self.tch
        {
            return tch.width(two_theta).0.max(MIN_FWHM);
        }
        if let Some(caglioti) = self.caglioti {
            return caglioti.fwhm(two_theta).max(MIN_FWHM);
        }
        if let Some(tch) = self.tch {
            return tch.width(two_theta).0.max(MIN_FWHM);
        }
        if self.fwhm.is_finite() && self.fwhm > 0.0 {
            self.fwhm.max(MIN_FWHM)
        } else {
            MIN_FWHM
        }
    }

    /// Pseudo-Voigt mixing at a Bragg angle: TCH mixing for a TCH shape, the
    /// configured `eta` for a pseudo-Voigt, and the pure-shape limits
    /// otherwise.
    #[must_use]
    pub fn eta_at(&self, two_theta: f64) -> f64 {
        match self.shape {
            LineShape::Tch => self.tch.map_or(0.0, |tch| tch.width(two_theta).1),
            LineShape::PseudoVoigt => self.eta.clamp(0.0, 1.0),
            LineShape::Gaussian => 0.0,
            LineShape::Lorentzian => 1.0,
        }
    }

    /// Unit-height profile value at `x` for a line centred at `center`.
    #[must_use]
    pub fn value(&self, center: f64, x: f64) -> f64 {
        self.value_scaled(center, x, 1.0, self.eta_at(center))
    }

    /// Unit-height profile value with an explicit width multiplier and mixing.
    ///
    /// `width_scale` multiplies the width law of this profile. For a TCH shape
    /// the mixing is always the TCH one; `eta` is used by the Gaussian,
    /// Lorentzian and pseudo-Voigt shapes.
    #[must_use]
    pub(crate) fn value_scaled(&self, center: f64, x: f64, width_scale: f64, eta: f64) -> f64 {
        let fwhm = (self.fwhm_at(center) * width_scale).max(MIN_FWHM);
        let eta = if self.shape == LineShape::Tch {
            self.eta_at(center)
        } else {
            eta
        };
        shape_value(self.shape, eta, center, fwhm, x)
            * self
                .asymmetry
                .map_or(1.0, |asymmetry| asymmetry.factor(x - center, center, fwhm))
    }

    /// Profile value scaled so the Gaussian part has unit area.
    ///
    /// The scale is `2 / (FWHM sqrt(pi / ln 2))`, the same convention as the
    /// Pawley and Le Bail engines, so a fitted intensity is on the PXRD
    /// reflection scale.
    #[must_use]
    pub fn area_normalised(&self, center: f64, x: f64) -> f64 {
        let fwhm = self.fwhm_at(center);
        self.value(center, x) * 2.0 / fwhm * (std::f64::consts::LN_2 / std::f64::consts::PI).sqrt()
    }

    /// Area-normalised value with an explicit width multiplier and mixing.
    #[must_use]
    pub(crate) fn area_normalised_scaled(
        &self,
        center: f64,
        x: f64,
        width_scale: f64,
        eta: f64,
    ) -> f64 {
        let fwhm = (self.fwhm_at(center) * width_scale).max(MIN_FWHM);
        self.value_scaled(center, x, width_scale, eta) * 2.0 / fwhm
            * (std::f64::consts::LN_2 / std::f64::consts::PI).sqrt()
    }

    /// Lorentz-polarisation factor at a Bragg angle in degrees.
    ///
    /// At the default `0.5` polarisation and an X-ray source this equals the
    /// traditional `(1 + cos^2 2theta) / (sin^2 theta cos theta)`.
    #[must_use]
    pub fn lorentz_polarization(&self, two_theta: f64) -> f64 {
        lorentz_polarization(
            self.radiation.kind,
            self.monochromator.polarization,
            two_theta,
        )
    }
}

/// Lorentz-polarisation factor of a probe at a Bragg angle in degrees.
///
/// X-rays use `2 (1 - p + p cos^2 2theta) / (sin^2 theta cos theta)`; neutrons
/// and electrons use the Lorentz factor alone.
pub(crate) fn lorentz_polarization(kind: RadiationKind, polarization: f64, two_theta: f64) -> f64 {
    let theta = 0.5 * two_theta.to_radians();
    let sine = theta.sin();
    let denominator = sine * sine * theta.cos();
    if denominator.abs() < 1e-12 {
        return 0.0;
    }
    match kind {
        RadiationKind::XRay => {
            let p = if polarization.is_finite() {
                polarization.clamp(0.0, 1.0)
            } else {
                0.5
            };
            let cos_two = two_theta.to_radians().cos();
            2.0 * (1.0 - p + p * cos_two * cos_two) / denominator
        }
        RadiationKind::Neutron | RadiationKind::Electron => 1.0 / denominator,
    }
}
