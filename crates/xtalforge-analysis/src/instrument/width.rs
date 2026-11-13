//! Width laws, line shapes and low-angle asymmetry of the instrument profile
//! (Task 14.1.1).
//!
//! # Documented formulas
//!
//! * Gaussian, unit height: `exp(-4 ln 2 d^2)` with `d = (x - x0) / FWHM`.
//! * Lorentzian, unit height: `1 / (1 + 4 d^2)`.
//! * Pseudo-Voigt: `eta L + (1 - eta) G`.
//! * TCH: with `r = Gamma_L / Gamma_G`,
//!
//!   ```text
//!   Gamma^5 = G^5 + 2.69269 G^4 L + 2.42843 G^3 L^2
//!             + 4.47163 G^2 L^3 + 0.07842 G L^4 + L^5
//!   eta = 1.36603 r - 0.47719 r^2 + 0.11116 r^3
//!   ```
//!
//!   `Gamma_G` comes from the Caglioti law and
//!   `Gamma_L = X tan(theta) + Y / cos(theta)`.
//! * FCJ-style asymmetry: for `x` on the low-angle side of the centre the
//!   shape is multiplied by `1 + a t` with `t = (x0 - x) / FWHM` clamped to
//!   four widths and `a = clamp(fcj1 + fcj2 tan(theta), 0, 1)`.
//!
//! # Approximations
//!
//! * The FCJ correction is a multiplicative low-angle tail, not the published
//!   convolution with the axial-divergence distribution; it has the correct
//!   sign and angle growth but does not renormalise the integrated intensity.
//!   The reference implementation documents the same simplification.
//! * The Caglioti law is evaluated at the line centre; a full profile with a
//!   width that varies inside the line is not modelled.

use serde::{Deserialize, Serialize};
use xtalforge_core::Error;

use crate::profiles::ProfileKind;

/// Smallest FWHM in degrees accepted by the profile evaluator.
pub(crate) const MIN_FWHM: f64 = 1e-4;
/// Largest low-angle tail, in FWHM units, before the asymmetry clamps.
const MAX_TAIL: f64 = 4.0;

/// Line shape of one reflection.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default, Serialize, Deserialize)]
pub enum LineShape {
    /// Pure Gaussian.
    #[default]
    Gaussian,
    /// Pure Lorentzian.
    Lorentzian,
    /// Pseudo-Voigt with a mixing parameter.
    PseudoVoigt,
    /// Thompson-Cox-Hastings; requires a [`Tch`] width model.
    Tch,
}

impl From<ProfileKind> for LineShape {
    fn from(kind: ProfileKind) -> Self {
        match kind {
            ProfileKind::Gaussian => Self::Gaussian,
            ProfileKind::Lorentzian => Self::Lorentzian,
            ProfileKind::PseudoVoigt => Self::PseudoVoigt,
        }
    }
}

/// Caglioti full-width law, in degrees squared:
/// `FWHM^2 = U tan^2(theta) + V tan(theta) + W`.
#[derive(Debug, Clone, Copy, PartialEq, Default, Serialize, Deserialize)]
pub struct Caglioti {
    /// `U` coefficient in degrees squared.
    pub u: f64,
    /// `V` coefficient in degrees squared.
    pub v: f64,
    /// `W` coefficient in degrees squared.
    pub w: f64,
}

impl Caglioti {
    /// FWHM in degrees at a Bragg angle in degrees two-theta.
    ///
    /// A non-finite angle or a negative discriminant returns the floor of
    /// `1e-4` degrees rather than `NaN`.
    #[must_use]
    pub fn fwhm(&self, two_theta: f64) -> f64 {
        if !two_theta.is_finite() {
            return MIN_FWHM;
        }
        let theta = 0.5 * two_theta.to_radians();
        let tangent = theta.tan();
        let squared = self.u * tangent * tangent + self.v * tangent + self.w;
        squared.max(MIN_FWHM * MIN_FWHM).sqrt().max(MIN_FWHM)
    }

    /// Validates the coefficients.
    ///
    /// # Errors
    ///
    /// Returns [`Error::InvalidInput`] for a non-finite coefficient.
    pub fn validate(&self) -> Result<(), Error> {
        if !self.u.is_finite() || !self.v.is_finite() || !self.w.is_finite() {
            return Err(Error::InvalidInput(
                "Caglioti coefficients must be finite".to_owned(),
            ));
        }
        Ok(())
    }
}

/// Thompson-Cox-Hastings width model.
///
/// The Gaussian part is a [`Caglioti`] law; the Lorentzian part is
/// `X tan(theta) + Y / cos(theta)`. Both parts are in degrees.
#[derive(Debug, Clone, Copy, PartialEq, Default, Serialize, Deserialize)]
pub struct Tch {
    /// Gaussian width law.
    pub gaussian: Caglioti,
    /// `X` coefficient in degrees.
    pub x: f64,
    /// `Y` coefficient in degrees.
    pub y: f64,
}

impl Tch {
    /// Thompson-Cox-Hastings `A` constant.
    const A: f64 = 2.692_69;
    /// `B` constant.
    const B: f64 = 2.428_43;
    /// `C` constant.
    const C: f64 = 4.471_63;
    /// `D` constant.
    const D: f64 = 0.078_42;

    /// Effective FWHM and pseudo-Voigt mixing at a Bragg angle.
    ///
    /// Returns `(MIN_FWHM, 0)` for a non-finite angle.
    #[must_use]
    pub fn width(&self, two_theta: f64) -> (f64, f64) {
        if !two_theta.is_finite() {
            return (MIN_FWHM, 0.0);
        }
        let theta = 0.5 * two_theta.to_radians();
        let gaussian = self.gaussian.fwhm(two_theta).max(MIN_FWHM);
        let lorentzian = (self.x * theta.tan() + self.y / theta.cos())
            .abs()
            .max(MIN_FWHM);
        let ratio = lorentzian / gaussian;
        let polynomial = 1.0
            + Self::A * ratio
            + Self::B * ratio * ratio
            + Self::C * ratio.powi(3)
            + Self::D * ratio.powi(4)
            + ratio.powi(5);
        let fwhm = gaussian * polynomial.powf(0.2);
        let eta = 1.366_03 * ratio - 0.477_19 * ratio * ratio + 0.111_16 * ratio.powi(3);
        (fwhm.max(MIN_FWHM), eta.clamp(0.0, 1.0))
    }

    /// Validates the model.
    ///
    /// # Errors
    ///
    /// Returns [`Error::InvalidInput`] for a non-finite coefficient.
    pub fn validate(&self) -> Result<(), Error> {
        self.gaussian.validate()?;
        if !self.x.is_finite() || !self.y.is_finite() {
            return Err(Error::InvalidInput(
                "TCH Lorentzian coefficients must be finite".to_owned(),
            ));
        }
        Ok(())
    }
}

/// Finger-Cox-Jephcoat style low-angle asymmetry.
///
/// Both parameters are dimensionless: `fcj1` is the constant part of the tail
/// and `fcj2` multiplies `tan(theta)`, so the asymmetry grows with angle as
/// the published correction does.
#[derive(Debug, Clone, Copy, PartialEq, Default, Serialize, Deserialize)]
pub struct FcjAsymmetry {
    /// Constant tail coefficient.
    pub fcj1: f64,
    /// `tan(theta)` tail coefficient.
    pub fcj2: f64,
}

impl FcjAsymmetry {
    /// Multiplier of the shape at an offset `delta = x - centre`.
    ///
    /// The multiplier is one on the high-angle side and at the centre, and
    /// `1 + a t` on the low-angle side, with `t` the offset in FWHM units
    /// clamped to four widths.
    #[must_use]
    pub fn factor(&self, delta: f64, two_theta: f64, fwhm: f64) -> f64 {
        if delta >= 0.0 || !delta.is_finite() || !two_theta.is_finite() || fwhm <= 0.0 {
            return 1.0;
        }
        let tail = (-delta / fwhm).min(MAX_TAIL);
        let theta = 0.5 * two_theta.to_radians();
        let slope = (self.fcj1 + self.fcj2 * theta.tan()).clamp(0.0, 1.0);
        (1.0 + slope * tail).max(0.0)
    }

    /// Validates the coefficients.
    ///
    /// # Errors
    ///
    /// Returns [`Error::InvalidInput`] for a non-finite or negative
    /// coefficient.
    pub fn validate(&self) -> Result<(), Error> {
        if !self.fcj1.is_finite() || !self.fcj2.is_finite() {
            return Err(Error::InvalidInput(
                "FCJ asymmetry coefficients must be finite".to_owned(),
            ));
        }
        if self.fcj1 < 0.0 || self.fcj2 < 0.0 {
            return Err(Error::InvalidInput(
                "FCJ asymmetry coefficients must be non-negative".to_owned(),
            ));
        }
        Ok(())
    }
}

/// Evaluates one line shape with a unit height.
pub(crate) fn shape_value(shape: LineShape, eta: f64, center: f64, fwhm: f64, x: f64) -> f64 {
    let width = fwhm.max(MIN_FWHM);
    let delta = (x - center) / width;
    let gaussian = (-4.0 * std::f64::consts::LN_2 * delta * delta).exp();
    let lorentzian = 1.0 / (1.0 + 4.0 * delta * delta);
    match shape {
        LineShape::Gaussian => gaussian,
        LineShape::Lorentzian => lorentzian,
        LineShape::PseudoVoigt | LineShape::Tch => {
            let eta = eta.clamp(0.0, 1.0);
            eta * lorentzian + (1.0 - eta) * gaussian
        }
    }
}
