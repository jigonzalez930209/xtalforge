//! Powder-averaged Raman spectra from engine-provided tensors (Task 19.4.2).
//!
//! **There is no native polarizability model in this crate.**  A Raman
//! intensity needs the derivative of the polarizability with respect to a
//! normal coordinate, which a force field does not provide.  This module reads
//! per-mode Raman tensors produced by an electronic-structure code and turns
//! them into the powder-averaged spectrum:
//!
//! ```text
//! alpha = (R_xx + R_yy + R_zz) / 3
//! beta^2 = ((R_xx - R_yy)^2 + (R_yy - R_zz)^2 + (R_zz - R_xx)^2
//!           + 6 (R_xy^2 + R_xz^2 + R_yz^2)) / 2
//! activity = (45 alpha^2 + 7 beta^2) / 45
//! ```
//!
//! The tensor is assumed symmetric; an engine tensor that is not is used as
//! written and the caller should symmetrise it first.  `activity` is the
//! orientation average of the scattering intensity in arbitrary units, not a
//! calibrated cross-section: no laser frequency, temperature or instrument
//! response enters, and the fourth-power scattering factor
//! `(nu_0 - nu)^4` is deliberately omitted.  Broadening is the same
//! unit-height Lorentzian or Gaussian used by the IR spectrum.
//!
//! Tensors are read with [`RamanTensorSet::from_text`] (the text/JSON formats
//! are documented on [`parse_tensor_file`](super::parse_tensor_file))
//! or supplied directly with [`RamanTensorSet::from_tensors`].  The set may
//! contain one tensor per mode (`3 N`, indices matching
//! [`VibrationalModes::frequencies`]) or one per real mode in ascending
//! frequency order; both layouts are unambiguous because the second count can
//! only match when no zero or imaginary mode exists.  Imaginary and zero modes
//! are never drawn.  The tensors are the mode polarizability derivatives as
//! reported by the engine, so no mass weighting is applied here.

use serde::{Deserialize, Serialize};

use crate::Error;

use super::modes::VibrationalModes;
use super::spectrum::{Broadening, broaden};
use super::tensors::parse_tensor_file;

/// Default full width at half maximum, in cm^-1.
pub const DEFAULT_RAMAN_FWHM: f64 = 10.0;

/// Default spacing of the broadened spectrum, in cm^-1.
pub const DEFAULT_RAMAN_STEP: f64 = 5.0;

/// Largest number of spectrum points a request may produce.
pub const MAX_RAMAN_POINTS: usize = 1_000_000;

/// Options for one Raman spectrum.
#[derive(Debug, Clone, Copy, PartialEq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct RamanOptions {
    /// Full width at half maximum in cm^-1.
    pub fwhm: f64,
    /// Spacing of the returned points in cm^-1.
    pub step: f64,
    /// Line shape.
    pub broadening: Broadening,
}

impl Default for RamanOptions {
    fn default() -> Self {
        Self {
            fwhm: DEFAULT_RAMAN_FWHM,
            step: DEFAULT_RAMAN_STEP,
            broadening: Broadening::default(),
        }
    }
}

impl RamanOptions {
    /// Checks that the options are usable.
    ///
    /// # Errors
    ///
    /// Returns [`Error::InvalidInput`] for a non-positive or non-finite width
    /// or spacing.
    pub fn validate(&self) -> Result<(), Error> {
        if !self.fwhm.is_finite() || self.fwhm <= 0.0 {
            return Err(Error::InvalidInput(format!(
                "Raman fwhm must be positive and finite, got {}",
                self.fwhm
            )));
        }
        if !self.step.is_finite() || self.step <= 0.0 {
            return Err(Error::InvalidInput(format!(
                "Raman step must be positive and finite, got {}",
                self.step
            )));
        }
        Ok(())
    }
}

/// Per-mode Raman tensors.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct RamanTensorSet {
    /// One `3 x 3` tensor per mode, row-major and symmetric in practice.
    pub tensors: Vec<[[f64; 3]; 3]>,
}

impl RamanTensorSet {
    /// Wraps explicit per-mode tensors.
    #[must_use]
    pub fn from_tensors(tensors: Vec<[[f64; 3]; 3]>) -> Self {
        Self { tensors }
    }

    /// Reads tensors from a simple text or JSON tensor file.
    ///
    /// The formats are documented on
    /// [`parse_tensor_file`](super::parse_tensor_file).  A full VASP
    /// or Quantum ESPRESSO output is not parsed; export its tensor block with
    /// one `3 x 3` tensor per mode in mode order.
    ///
    /// # Errors
    ///
    /// Returns the parser errors of
    /// [`parse_tensor_file`](super::parse_tensor_file).
    pub fn from_text(text: &str) -> Result<Self, Error> {
        Ok(Self {
            tensors: parse_tensor_file(text)?.tensors,
        })
    }

    /// Number of stored tensors.
    #[must_use]
    pub fn len(&self) -> usize {
        self.tensors.len()
    }

    /// Whether no tensor is stored.
    #[must_use]
    pub fn is_empty(&self) -> bool {
        self.tensors.is_empty()
    }
}

/// One Raman stick.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct RamanPeak {
    /// Index of the mode in the source [`VibrationalModes`].
    pub mode: usize,
    /// Frequency in cm^-1.
    pub frequency: f64,
    /// Powder-averaged Raman activity in arbitrary units.
    pub activity: f64,
}

/// Raman stick spectrum and its broadened representation.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct RamanSpectrum {
    /// Broadened spectrum as `[wavenumber (cm^-1), activity]` points.
    pub points: Vec<[f64; 2]>,
    /// One stick per real mode, in ascending frequency order.
    pub peaks: Vec<RamanPeak>,
}

impl RamanSpectrum {
    /// Largest broadened activity, or zero for an empty spectrum.
    #[must_use]
    pub fn max_intensity(&self) -> f64 {
        self.points
            .iter()
            .filter_map(|point| point.get(1).copied())
            .fold(0.0_f64, f64::max)
    }

    /// Largest stick activity, or zero when there are no peaks.
    #[must_use]
    pub fn max_activity(&self) -> f64 {
        self.peaks
            .iter()
            .map(|peak| peak.activity)
            .fold(0.0_f64, f64::max)
    }
}

/// Powder-orientation average of one Raman tensor.
///
/// Uses the standard perpendicular/parallel invariant formula documented in
/// the module.  The result is the scattering activity in arbitrary units, so
/// only ratios are meaningful.
#[must_use]
pub fn powder_average(tensor: &[[f64; 3]; 3]) -> f64 {
    let xx = tensor[0][0];
    let yy = tensor[1][1];
    let zz = tensor[2][2];
    let alpha = (xx + yy + zz) / 3.0;
    let beta_squared = 0.5
        * ((xx - yy) * (xx - yy)
            + (yy - zz) * (yy - zz)
            + (zz - xx) * (zz - xx)
            + 6.0
                * (tensor[0][1] * tensor[0][1]
                    + tensor[0][2] * tensor[0][2]
                    + tensor[1][2] * tensor[1][2]));
    (45.0 * alpha * alpha + 7.0 * beta_squared) / 45.0
}

/// Builds a powder-averaged Raman spectrum from mode frequencies and tensors.
///
/// `tensors` holds one tensor per mode or one per real mode; the mapping rule
/// is documented in the module comment.
///
/// # Errors
///
/// Returns [`Error::InvalidInput`] for invalid options, a tensor count that
/// matches neither layout, a non-finite tensor component, or a grid larger
/// than [`MAX_RAMAN_POINTS`].
pub fn raman_spectrum(
    modes: &VibrationalModes,
    tensors: &RamanTensorSet,
    options: &RamanOptions,
) -> Result<RamanSpectrum, Error> {
    options.validate()?;
    let real = modes.real_modes().collect::<Vec<_>>();
    let total = modes.n_modes();
    if tensors.is_empty() {
        if real.is_empty() {
            return Ok(RamanSpectrum {
                points: Vec::new(),
                peaks: Vec::new(),
            });
        }
        return Err(Error::InvalidInput(
            "Raman tensors are required for the real modes".to_owned(),
        ));
    }
    let per_mode = tensors.len() == total;
    if !per_mode && tensors.len() != real.len() {
        return Err(Error::InvalidInput(format!(
            "got {} Raman tensors for {total} modes ({} real)",
            tensors.len(),
            real.len()
        )));
    }
    let mut peaks = Vec::with_capacity(real.len());
    for (position, (mode, frequency)) in real.iter().enumerate() {
        let index = if per_mode { *mode } else { position };
        let Some(tensor) = tensors.tensors.get(index) else {
            return Err(Error::InvalidInput(format!(
                "Raman tensor {index} is missing"
            )));
        };
        if tensor.iter().flatten().any(|value| !value.is_finite()) {
            return Err(Error::InvalidInput(format!(
                "Raman tensor {index} has a non-finite component"
            )));
        }
        peaks.push(RamanPeak {
            mode: *mode,
            frequency: *frequency,
            activity: powder_average(tensor),
        });
    }
    peaks.sort_by(|first, second| {
        first
            .frequency
            .partial_cmp(&second.frequency)
            .unwrap_or(std::cmp::Ordering::Equal)
            .then(first.mode.cmp(&second.mode))
    });
    let sticks = peaks
        .iter()
        .map(|peak| (peak.frequency, peak.activity))
        .collect::<Vec<_>>();
    let points = broaden(
        &sticks,
        options.fwhm,
        options.step,
        options.broadening,
        MAX_RAMAN_POINTS,
        "Raman",
    )?;
    Ok(RamanSpectrum { points, peaks })
}
