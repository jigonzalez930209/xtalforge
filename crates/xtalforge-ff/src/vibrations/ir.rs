//! IR intensities and a broadened spectrum from Gamma-point modes (Tasks 4.8.2
//! and 19.4.1).
//!
//! The intensity of a normal mode is proportional to the square of the dipole
//! derivative along that mode.  [`ir_spectrum`] uses the documented fixed-charge
//! approximation
//!
//! ```text
//! mu = sum_i q_i r_i
//! ```
//!
//! with QEq charges `q_i` from [`crate::qeq_charges`] and the derivative taken
//! at the reference geometry.  This is the reference implementation's `fixed`
//! charge model: it captures the charge separation of a polar bond but not the
//! charge that flows when that bond stretches.  Because the charges are held
//! fixed the dipole is linear in the coordinates, so this derivative is
//! analytic; a central finite difference of the same dipole is numerically
//! identical.  A finite-difference `flux` model that re-solves EQeq at every
//! displaced geometry remains future work.
//!
//! [`ir_spectrum_with_born`](super::ir_spectrum_with_born) replaces the scalar
//! charges by per-atom Born effective-charge tensors
//! [`BornCharges`](super::BornCharges), the modern form of the same derivative:
//!
//! ```text
//! p_a = sum_i sum_b Z*_{i,a b} u_{i,b}
//! ```
//!
//! so the intensity is `|p|^2 / sum_i m_i |u_i|^2` in the same units.  The
//! tensors can come from a periodic engine (VASP, Quantum ESPRESSO) through
//! [`BornCharges::from_text`](super::BornCharges::from_text), from a JSON
//! tensor file, or from diagonal QEq or EQeq charges through
//! [`BornCharges::qeq`](super::BornCharges::qeq) and
//! [`BornCharges::eqeq`](super::BornCharges::eqeq).  The diagonal shortcut is
//! the fixed-charge model written as a tensor; a real engine tensor also
//! carries the charge flow that the fixed-charge path misses.
//!
//! The fixed-charge path shifts the charges by their mean before use and the
//! Born path offers the tensor analogue as
//! [`BornCharges::neutralised`](super::BornCharges::neutralised).  A uniform
//! translation then moves no dipole and the three acoustic modes carry no
//! intensity.  Neither shift is applied silently: each changes a non-neutral
//! cell and is therefore an explicit approximation of the neutral-cell limit.
//!
//! Units: frequencies in cm^-1, stick intensities in km/mol, broadened points
//! as `[wavenumber, intensity]` with the profile area not renormalised.  Only
//! real modes above [`super::ZERO_MODE_TOLERANCE_CM`] are included; imaginary
//! and zero modes are not drawn, because a saddle or a rigid motion is not an
//! absorption line.

use serde::{Deserialize, Serialize};
use xtalforge_core::{Structure, element_info};

use crate::Error;
use crate::charges::{ChargeOptions, qeq_charges};

use super::constants::E2_PER_AMU_TO_KM_PER_MOL;
use super::modes::VibrationalModes;
use super::spectrum::{Broadening, broaden};

/// Default full width at half maximum, in cm^-1.
pub const DEFAULT_IR_FWHM: f64 = 10.0;

/// Default spacing of the broadened spectrum, in cm^-1.
pub const DEFAULT_IR_STEP: f64 = 5.0;

/// Largest number of spectrum points a request may produce.
pub const MAX_IR_POINTS: usize = 1_000_000;

/// Options for one IR spectrum.
#[derive(Debug, Clone, Copy, PartialEq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct IrOptions {
    /// Full width at half maximum in cm^-1.
    pub fwhm: f64,
    /// Spacing of the returned points in cm^-1.
    pub step: f64,
    /// Line shape.
    pub broadening: Broadening,
}

impl Default for IrOptions {
    fn default() -> Self {
        Self {
            fwhm: DEFAULT_IR_FWHM,
            step: DEFAULT_IR_STEP,
            broadening: Broadening::default(),
        }
    }
}

impl IrOptions {
    /// Checks that the options are usable.
    ///
    /// # Errors
    ///
    /// Returns [`Error::InvalidInput`] for a non-positive or non-finite width
    /// or spacing.
    pub fn validate(&self) -> Result<(), Error> {
        if !self.fwhm.is_finite() || self.fwhm <= 0.0 {
            return Err(Error::InvalidInput(format!(
                "IR fwhm must be positive and finite, got {}",
                self.fwhm
            )));
        }
        if !self.step.is_finite() || self.step <= 0.0 {
            return Err(Error::InvalidInput(format!(
                "IR step must be positive and finite, got {}",
                self.step
            )));
        }
        Ok(())
    }
}

/// One stick in the IR spectrum.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct IrPeak {
    /// Index of the mode in the source [`VibrationalModes`].
    pub mode: usize,
    /// Frequency in cm^-1.
    pub frequency: f64,
    /// Integrated intensity in km/mol.
    pub intensity: f64,
}

/// IR stick spectrum and its broadened representation.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct IrSpectrum {
    /// Broadened spectrum as `[wavenumber (cm^-1), intensity]` points.
    pub points: Vec<[f64; 2]>,
    /// One stick per real mode, in ascending frequency order.
    pub peaks: Vec<IrPeak>,
}

impl IrSpectrum {
    /// Largest broadened intensity, or zero for an empty spectrum.
    #[must_use]
    pub fn max_intensity(&self) -> f64 {
        self.points
            .iter()
            .filter_map(|point| point.get(1).copied())
            .fold(0.0_f64, f64::max)
    }

    /// Largest stick intensity, or zero when there are no peaks.
    #[must_use]
    pub fn max_stick_intensity(&self) -> f64 {
        self.peaks
            .iter()
            .map(|peak| peak.intensity)
            .fold(0.0_f64, f64::max)
    }
}

/// Builds an IR spectrum from a molecule's Gamma-point modes.
///
/// The charges are the default QEq ones; a future charge-model option can
/// supply EQeq or a force field's own charges without changing this result
/// type.  The structure and the modes must describe the same P1 atoms.
///
/// # Errors
///
/// Returns [`Error::InvalidInput`] for invalid options, a non-P1 structure, a
/// mode/atom mismatch, a non-positive mass or an excessively large grid, and
/// the charge solver's errors when QEq cannot be evaluated.
pub fn ir_spectrum(
    structure: &Structure,
    modes: &VibrationalModes,
    options: &IrOptions,
) -> Result<IrSpectrum, Error> {
    options.validate()?;
    if !structure.is_p1() {
        return Err(Error::InvalidInput(
            "IR intensities need a P1 structure; expand the asymmetric unit first".to_owned(),
        ));
    }
    if modes.n_atoms() != structure.sites().len() {
        return Err(Error::InvalidInput(format!(
            "modes describe {} atoms but the structure has {}",
            modes.n_atoms(),
            structure.sites().len()
        )));
    }

    let charges = qeq_charges(structure, &ChargeOptions::default())?;
    if charges.per_atom.len() != structure.sites().len() {
        return Err(Error::InvalidInput(format!(
            "QEq returned {} charges for {} modes",
            charges.per_atom.len(),
            structure.sites().len()
        )));
    }
    let count = charges.per_atom.len();
    let mean = if count == 0 {
        0.0
    } else {
        charges.per_atom.iter().sum::<f64>() / count as f64
    };
    let charges: Vec<f64> = charges
        .per_atom
        .iter()
        .map(|charge| charge - mean)
        .collect();

    let mut masses = Vec::with_capacity(count);
    for site in structure.sites() {
        let element = element_info(site.element())?;
        if element.mass <= 0.0 {
            return Err(Error::InvalidInput(format!(
                "element {} has no positive mass",
                site.element()
            )));
        }
        masses.push(element.mass);
    }

    let mut peaks = Vec::new();
    for (mode, frequency) in modes.real_modes() {
        let Some(displacements) = modes.modes.get(mode) else {
            continue;
        };
        if displacements.len() != count {
            return Err(Error::InvalidInput(format!(
                "mode {mode} has {} atoms, expected {count}",
                displacements.len()
            )));
        }
        let mut dipole = [0.0; 3];
        let mut mass_norm = 0.0;
        for ((point, charge), mass) in displacements.iter().zip(&charges).zip(&masses) {
            let length_squared = point[0] * point[0] + point[1] * point[1] + point[2] * point[2];
            mass_norm += mass * length_squared;
            dipole[0] += charge * point[0];
            dipole[1] += charge * point[1];
            dipole[2] += charge * point[2];
        }
        if mass_norm <= 0.0 {
            return Err(Error::InvalidInput(format!(
                "mode {mode} has zero mass-weighted norm"
            )));
        }
        let intensity = E2_PER_AMU_TO_KM_PER_MOL
            * (dipole[0] * dipole[0] + dipole[1] * dipole[1] + dipole[2] * dipole[2])
            / mass_norm;
        peaks.push(IrPeak {
            mode,
            frequency,
            intensity,
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
        .map(|peak| (peak.frequency, peak.intensity))
        .collect::<Vec<_>>();
    let points = broaden(
        &sticks,
        options.fwhm,
        options.step,
        options.broadening,
        MAX_IR_POINTS,
        "IR",
    )?;
    Ok(IrSpectrum { points, peaks })
}
