//! Phonons on a q-mesh (Task 4.8.3).
//!
//! The Gamma-point path uses the primitive-cell Hessian directly. A q-mesh
//! builds a supercell, takes the finite-difference Hessian of that supercell,
//! folds the force constants back onto the primitive cell and diagonalizes the
//! dynamical matrix at every Monkhorst-Pack point. The phonon DOS is a
//! Gaussian-smeared histogram over those frequencies.
//!
//! The Gamma result remains the reference for the acoustic sum rule; a
//! supercell calculation is checked against it by setting `supercell` to
//! `[1, 1, 1]` and `q_grid` to `[1, 1, 1]`.

use serde::{Deserialize, Serialize};
use xtalforge_core::Structure;

use crate::Error;
use crate::calculator::Calculator;

use super::hessian::{DisplacementMode, HessianOptions, hessian};
use super::modes::{VibrationalModes, vibrational_modes};
use super::qmesh::ForceConstants;

pub use super::qmesh::{DosOptions, PhononDos};

/// Options for a phonon calculation.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct PhononOptions {
    /// Monkhorst-Pack grid repeats.
    pub q_grid: [u32; 3],
    /// Supercell repeats used for the force constants.
    pub supercell: [u32; 3],
    /// Explicit q-points in fractional reciprocal coordinates.
    ///
    /// When non-empty these replace the Monkhorst-Pack grid, so a caller can
    /// pass the q-points of a band path. A supercell larger than `[1, 1, 1]`
    /// is still required.
    #[serde(default)]
    pub q_path: Vec<[f64; 3]>,
    /// Phonon DOS histogram, when requested.
    pub dos: Option<DosOptions>,
    /// Hessian settings, including the finite-difference step.
    pub hessian: HessianOptions,
}

impl Default for PhononOptions {
    fn default() -> Self {
        Self {
            q_grid: [1, 1, 1],
            supercell: [1, 1, 1],
            q_path: Vec::new(),
            dos: None,
            hessian: HessianOptions::crystal(),
        }
    }
}

impl PhononOptions {
    /// Checks that the options are usable.
    ///
    /// # Errors
    ///
    /// Returns [`Error::InvalidInput`] when a repeat is zero, when a q-mesh
    /// has no supercell to interpolate from, when the Hessian settings are
    /// invalid, or when the DOS options are invalid.
    pub fn validate(&self) -> Result<(), Error> {
        self.hessian.validate()?;
        if self.q_grid.contains(&0) {
            return Err(Error::InvalidInput(
                "phonon q-grid repeats must be positive".to_owned(),
            ));
        }
        if self.supercell.contains(&0) {
            return Err(Error::InvalidInput(
                "phonon supercell repeats must be positive".to_owned(),
            ));
        }
        if (!self.q_path.is_empty() || self.q_grid != [1, 1, 1]) && self.supercell == [1, 1, 1] {
            return Err(Error::InvalidInput(
                "a phonon q-mesh or q-path needs a supercell larger than [1, 1, 1]".to_owned(),
            ));
        }
        if !self.q_path.iter().flatten().all(|value| value.is_finite()) {
            return Err(Error::InvalidInput(
                "phonon q-path values must be finite".to_owned(),
            ));
        }
        if self.supercell == [1, 1, 1] && !self.q_path.is_empty() {
            return Err(Error::InvalidInput(
                "a phonon q-path needs a supercell larger than [1, 1, 1]".to_owned(),
            ));
        }
        if let Some(dos) = &self.dos {
            dos.validate()?;
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

/// Phonon result over one or more q-points.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct PhononResult {
    /// Sampled q-points in fractional reciprocal coordinates.
    pub q_points: Vec<[f64; 3]>,
    /// Frequencies in cm^-1 for each sampled q-point, one vector per q-point.
    pub frequencies: Vec<Vec<f64>>,
    /// Gamma-point modes, including the zero-mode and imaginary indices.
    pub modes: VibrationalModes,
    /// Phonon density of states, when requested.
    pub dos: Option<PhononDos>,
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

/// Calculates phonons of a periodic structure.
///
/// A `[1, 1, 1]` grid with a `[1, 1, 1]` supercell is the primitive-cell
/// Gamma point. Any other grid builds the supercell force constants and
/// samples the Monkhorst-Pack mesh.
///
/// # Errors
///
/// Returns [`Error::InvalidInput`] for a non-P1 structure, invalid options, or
/// any error from the Hessian, force constants, mode analysis or DOS.
pub fn phonons(
    calculator: &dyn Calculator,
    structure: &Structure,
    options: &PhononOptions,
) -> Result<PhononResult, Error> {
    options.validate()?;
    if !structure.is_p1() {
        return Err(Error::InvalidInput(
            "phonons need a P1 structure; expand the asymmetric unit first".to_owned(),
        ));
    }
    let hessian = hessian(calculator, structure, &options.hessian)?;
    let modes = vibrational_modes(&hessian, structure, &options.hessian)?;
    if options.q_path.is_empty() && options.q_grid == [1, 1, 1] && options.supercell == [1, 1, 1] {
        let dos = options
            .dos
            .as_ref()
            .map(|dos| super::qmesh::dos(std::slice::from_ref(&modes.frequencies), dos))
            .transpose()?;
        return Ok(PhononResult {
            q_points: vec![[0.0, 0.0, 0.0]],
            frequencies: vec![modes.frequencies.clone()],
            modes,
            dos,
            note: "Gamma point: the primitive-cell dynamical matrix.".to_owned(),
        });
    }

    let constants: ForceConstants =
        super::qmesh::force_constants(calculator, structure, options.supercell, &options.hessian)?;
    let q_points = if options.q_path.is_empty() {
        super::qmesh::q_points(options.q_grid)
    } else {
        options.q_path.clone()
    };
    let mut frequencies = Vec::with_capacity(q_points.len());
    for q in &q_points {
        frequencies.push(super::qmesh::frequencies(&constants, *q)?);
    }
    let dos = options
        .dos
        .as_ref()
        .map(|dos| super::qmesh::dos(&frequencies, dos))
        .transpose()?;
    Ok(PhononResult {
        q_points,
        frequencies,
        modes,
        dos,
        note: format!(
            "{} q-points from a {:?} grid on a {:?} supercell.",
            options
                .q_grid
                .iter()
                .map(|value| *value as usize)
                .product::<usize>(),
            options.q_grid,
            options.supercell
        ),
    })
}
