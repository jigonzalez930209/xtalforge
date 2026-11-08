//! Gamma-point phonons (Task 4.8.3, core subset).
//!
//! The full task calls for a symmetry-reduced supercell, force constants folded
//! back onto the primitive cell, a Monkhorst-Pack q-mesh and a phonon DOS.  That
//! is not implemented here: this module exposes the Gamma-point dynamical
//! matrix, which is exactly the mass-weighted Hessian of the primitive cell,
//! and reports it through [`PhononResult`].  The q-grid is therefore restricted
//! to `[1, 1, 1]`; asking for anything else returns a typed error instead of a
//! silently wrong dispersion.
//!
//! **Why this is still useful.**  Gamma-point frequencies are what a
//! thermochemistry run at a single q-point consumes, and they are the acoustic
//! sum-rule check a supercell calculation is validated against.  The three
//! acoustic modes are projected to zero by the mode analysis.
//!
//! **Pending, documented, not approximated.**  A q-mesh needs force constants
//! between atoms in neighbouring cells; the primitive-cell Hessian alone
//! cannot produce them.  The next implementation step is a supercell Hessian
//! folded onto the primitive cell, as Crystal-Builder delegates to phonopy.
//! Until then this module reports only Gamma and says so in every result note.

use serde::{Deserialize, Serialize};
use xtalforge_core::Structure;

use crate::Error;
use crate::calculator::Calculator;

use super::hessian::{DisplacementMode, HessianOptions, hessian};
use super::modes::{VibrationalModes, vibrational_modes};

/// Options for a phonon calculation.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct PhononOptions {
    /// Monkhorst-Pack grid repeats.  Only `[1, 1, 1]` is implemented.
    pub q_grid: [u32; 3],
    /// Hessian settings, including the finite-difference step.
    pub hessian: HessianOptions,
}

impl Default for PhononOptions {
    fn default() -> Self {
        Self {
            q_grid: [1, 1, 1],
            hessian: HessianOptions::crystal(),
        }
    }
}

impl PhononOptions {
    /// Checks that the options are usable.
    ///
    /// # Errors
    ///
    /// Returns [`Error::InvalidInput`] when a repeat is zero, when the grid is
    /// not the implemented `[1, 1, 1]` Gamma point, or when the Hessian
    /// settings are invalid.
    pub fn validate(&self) -> Result<(), Error> {
        self.hessian.validate()?;
        if self.q_grid.contains(&0) {
            return Err(Error::InvalidInput(
                "phonon q-grid repeats must be positive".to_owned(),
            ));
        }
        if self.q_grid != [1, 1, 1] {
            return Err(Error::InvalidInput(format!(
                "phonon q-grids beyond Gamma are not implemented yet; got {:?}, use [1, 1, 1]",
                self.q_grid
            )));
        }
        if self.hessian.displacement != DisplacementMode::Crystal {
            return Err(Error::InvalidInput(
                "phonons are a periodic-crystal calculation; use DisplacementMode::Crystal"
                    .to_owned(),
            ));
        }
        Ok(())
    }
}

/// Gamma-point phonon result.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct PhononResult {
    /// Sampled q-points in fractional reciprocal coordinates.
    ///
    /// Only `[[0, 0, 0]]` is returned by the implemented Gamma-point path.
    pub q_points: Vec<[f64; 3]>,
    /// Frequencies in cm^-1 for each sampled q-point, one vector per q-point.
    pub frequencies: Vec<Vec<f64>>,
    /// Gamma-point modes, including the zero-mode and imaginary indices.
    pub modes: VibrationalModes,
    /// Human-readable note about what was calculated.
    pub note: String,
}

impl PhononResult {
    /// Frequencies of the first (Gamma) q-point, or an empty slice.
    #[must_use]
    pub fn gamma_frequencies(&self) -> &[f64] {
        self.frequencies.first().map(Vec::as_slice).unwrap_or(&[])
    }
}

/// Calculates Gamma-point phonons of a periodic structure.
///
/// # Errors
///
/// Returns [`Error::InvalidInput`] for a non-P1 structure, a q-grid other than
/// `[1, 1, 1]`, invalid Hessian options, or any error from the Hessian or mode
/// analysis.
pub fn phonons(
    calculator: &dyn Calculator,
    structure: &Structure,
    options: &PhononOptions,
) -> Result<PhononResult, Error> {
    options.validate()?;
    if !structure.is_p1() {
        return Err(Error::InvalidInput(
            "Gamma-point phonons need a P1 structure; expand the asymmetric unit first".to_owned(),
        ));
    }
    let hessian = hessian(calculator, structure, &options.hessian)?;
    let modes = vibrational_modes(&hessian, structure, &options.hessian)?;
    Ok(PhononResult {
        q_points: vec![[0.0, 0.0, 0.0]],
        frequencies: vec![modes.frequencies.clone()],
        modes,
        note: "Gamma point only: the primitive-cell dynamical matrix. A symmetry-reduced supercell, q-mesh and phonon DOS are pending (Task 4.8.3).".to_owned(),
    })
}
