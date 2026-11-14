//! Geometry, axis and options of the neutron powder engine (Task 19.3.1).

use serde::{Deserialize, Serialize};
use xtalforge_core::Error;

use crate::instrument::{InstrumentProfile, LineShape, RadiationKind, RadiationSource};

use super::sample::{MAX_PATTERN_POINTS, point_count};

/// Geometry of a neutron powder measurement.
#[derive(Debug, Clone, Copy, PartialEq, Serialize, Deserialize)]
pub enum NeutronGeometry {
    /// Constant wavelength in Angstrom; the pattern is sampled on `2 theta`.
    ConstantWavelength {
        /// Incident wavelength in Angstrom.
        wavelength: f64,
    },
    /// Time of flight of one detector bank; the pattern is sampled in
    /// microseconds.
    TimeOfFlight {
        /// Total flight path in metre.
        total_path_m: f64,
        /// Detector bank angle `2 theta` in degrees.
        two_theta_deg: f64,
        /// Relative `Delta d / d` resolution.
        relative_d_fwhm: f64,
    },
}

impl NeutronGeometry {
    pub(crate) fn validate(&self) -> Result<(), Error> {
        match *self {
            Self::ConstantWavelength { wavelength } => {
                if !wavelength.is_finite() || wavelength <= 0.0 {
                    return Err(Error::InvalidInput(format!(
                        "neutron wavelength must be positive and finite, got {wavelength}"
                    )));
                }
            }
            Self::TimeOfFlight {
                total_path_m,
                two_theta_deg,
                relative_d_fwhm,
            } => {
                if !total_path_m.is_finite() || total_path_m <= 0.0 {
                    return Err(Error::InvalidInput(format!(
                        "time-of-flight path length must be positive and finite, got {total_path_m}"
                    )));
                }
                if !two_theta_deg.is_finite() || two_theta_deg <= 0.0 || two_theta_deg >= 180.0 {
                    return Err(Error::InvalidInput(format!(
                        "time-of-flight bank angle must be in (0, 180) degrees, got {two_theta_deg}"
                    )));
                }
                if !relative_d_fwhm.is_finite() || relative_d_fwhm <= 0.0 {
                    return Err(Error::InvalidInput(format!(
                        "relative d resolution must be positive and finite, got {relative_d_fwhm}"
                    )));
                }
            }
        }
        Ok(())
    }
}

/// Sampling axis of a [`NeutronPattern`](super::NeutronPattern).
#[derive(Debug, Clone, Copy, PartialEq, Serialize, Deserialize)]
pub enum NeutronAxis {
    /// `2 theta` in degrees at a constant wavelength.
    TwoTheta {
        /// Wavelength in Angstrom.
        wavelength: f64,
    },
    /// Time of flight in microseconds at one detector bank.
    TimeOfFlight {
        /// Total flight path in metre.
        total_path_m: f64,
        /// Bank angle `2 theta` in degrees.
        two_theta_deg: f64,
    },
}

impl NeutronAxis {
    /// Unit of the `x` coordinate.
    #[must_use]
    pub fn unit(&self) -> &'static str {
        match self {
            Self::TwoTheta { .. } => "degrees 2theta",
            Self::TimeOfFlight { .. } => "microseconds",
        }
    }
}

/// Options for [`neutron_powder`](super::neutron_powder).
#[derive(Debug, Clone, Copy, PartialEq, Serialize, Deserialize)]
pub struct NeutronOptions {
    /// Measurement geometry.
    pub geometry: NeutronGeometry,
    /// Phase 7 instrument profile. Its radiation must be
    /// [`RadiationKind::Neutron`]; the wavelength of the geometry sets the
    /// line positions.
    pub profile: InstrumentProfile,
    /// Inclusive sampling range in the axis unit.
    pub range: [f64; 2],
    /// Sampling step in the axis unit.
    pub step: f64,
}

impl Default for NeutronOptions {
    fn default() -> Self {
        Self {
            geometry: NeutronGeometry::ConstantWavelength { wavelength: 1.54 },
            profile: InstrumentProfile {
                radiation: RadiationSource {
                    kind: RadiationKind::Neutron,
                    wavelength: 1.54,
                },
                shape: LineShape::PseudoVoigt,
                fwhm: 0.2,
                ..InstrumentProfile::default()
            },
            range: [5.0, 120.0],
            step: 0.05,
        }
    }
}

impl NeutronOptions {
    /// Validates the geometry, profile, range and step.
    ///
    /// # Errors
    ///
    /// Returns [`Error::InvalidInput`] for a non-neutron profile, an invalid
    /// geometry, a reversed or out-of-range window, a non-positive step, or a
    /// pattern with more than two million points.
    pub fn validate(&self) -> Result<(), Error> {
        self.profile.validate()?;
        if self.profile.radiation.kind != RadiationKind::Neutron {
            return Err(Error::InvalidInput(
                "neutron powder diffraction needs a neutron instrument profile".to_owned(),
            ));
        }
        self.geometry.validate()?;
        let [minimum, maximum] = self.range;
        if !minimum.is_finite() || !maximum.is_finite() || maximum <= minimum {
            return Err(Error::InvalidInput(format!(
                "pattern range must be finite and increasing, got {minimum}..{maximum}"
            )));
        }
        if let NeutronGeometry::ConstantWavelength { .. } = self.geometry
            && (minimum <= 0.0 || maximum > 180.0)
        {
            return Err(Error::InvalidInput(format!(
                "2 theta range must satisfy 0 < min < max <= 180, got {minimum}..{maximum}"
            )));
        }
        if let NeutronGeometry::TimeOfFlight { .. } = self.geometry
            && minimum <= 0.0
        {
            return Err(Error::InvalidInput(format!(
                "time-of-flight range must be positive, got {minimum}..{maximum}"
            )));
        }
        if !self.step.is_finite() || self.step <= 0.0 {
            return Err(Error::InvalidInput(format!(
                "pattern step must be positive and finite, got {}",
                self.step
            )));
        }
        let count = point_count(self.range, self.step)?;
        if count > MAX_PATTERN_POINTS {
            return Err(Error::InvalidInput(format!(
                "pattern would have {count} points; use a larger step or a shorter range"
            )));
        }
        Ok(())
    }
}
