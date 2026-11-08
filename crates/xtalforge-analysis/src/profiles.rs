//! Peak profiles and background models (Task 7.4.2/7.4.3).

use serde::{Deserialize, Serialize};
use xtalforge_core::Error;

/// Peak shape.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default, Serialize, Deserialize)]
pub enum ProfileKind {
    /// Gaussian.
    #[default]
    Gaussian,
    /// Lorentzian.
    Lorentzian,
    /// Pseudo-Voigt with a mixing parameter.
    PseudoVoigt,
}

/// A peak profile at a position with a FWHM.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Profile {
    /// Shape.
    pub kind: ProfileKind,
    /// Centre in degrees two-theta.
    pub center: f64,
    /// Full width at half maximum in degrees.
    pub fwhm: f64,
    /// Pseudo-Voigt mixing (0 = Gaussian, 1 = Lorentzian).
    pub eta: f64,
    /// Height at the centre.
    pub height: f64,
}

impl Profile {
    /// Value at `x`.
    #[must_use]
    pub fn value(&self, x: f64) -> f64 {
        let width = self.fwhm.max(1e-9);
        let delta = (x - self.center) / width;
        let gaussian = (-4.0 * std::f64::consts::LN_2 * delta * delta).exp();
        let lorentzian = 1.0 / (1.0 + 4.0 * delta * delta);
        let shape = match self.kind {
            ProfileKind::Gaussian => gaussian,
            ProfileKind::Lorentzian => lorentzian,
            ProfileKind::PseudoVoigt => {
                let eta = self.eta.clamp(0.0, 1.0);
                eta * lorentzian + (1.0 - eta) * gaussian
            }
        };
        self.height * shape
    }

    /// Samples the profile on a range.
    ///
    /// # Errors
    ///
    /// Returns [`Error::InvalidInput`] for a non-positive step or a reversed
    /// range.
    pub fn sample(&self, range: [f64; 2], step: f64) -> Result<Vec<(f64, f64)>, Error> {
        if !step.is_finite() || step <= 0.0 || range[0] > range[1] {
            return Err(Error::InvalidInput(
                "profile range and step must be valid".to_owned(),
            ));
        }
        let count = ((range[1] - range[0]) / step).floor() as usize + 1;
        Ok((0..count)
            .map(|index| {
                let x = range[0] + index as f64 * step;
                (x, self.value(x))
            })
            .collect())
    }
}

/// Background model.
#[derive(Debug, Clone, Copy, PartialEq, Serialize, Deserialize)]
pub enum Background {
    /// Constant.
    Constant(f64),
    /// Linear `c0 + c1 * x`.
    Linear(f64, f64),
    /// Quadratic `c0 + c1 * x + c2 * x^2`.
    Quadratic(f64, f64, f64),
}

impl Background {
    /// Value at `x`.
    #[must_use]
    pub fn value(&self, x: f64) -> f64 {
        match *self {
            Self::Constant(c0) => c0,
            Self::Linear(c0, c1) => c0 + c1 * x,
            Self::Quadratic(c0, c1, c2) => c0 + c1 * x + c2 * x * x,
        }
    }
}
