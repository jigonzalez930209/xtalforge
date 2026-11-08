//! Mass weighting, rigid-motion projection and normal modes (Task 4.8.1).
//!
//! The physical Hessian in eV/Angstrom^2 is divided by the square roots of the
//! atomic masses (amu) and diagonalized.  The eigenvalues are in
//! eV/(Angstrom^2 amu); a positive one is a real harmonic frequency, a negative
//! one is reported as a negative wavenumber and is an imaginary mode.
//!
//! **Rigid motions are projected, not left to finite differences.**  For a
//! crystal the three uniform translations are exact acoustic zero modes; for a
//! molecule the three translations and (for a non-linear molecule) three
//! rotations are removed.  The projection uses the masses and the Cartesian
//! geometry, not a symmetry library, and is the same idea as the reference
//! implementation in Crystal-Builder `xtal/ff/vibrations.py`.
//!
//! **Documented approximation.**  Projecting rotations is physically exact
//! for an isolated molecule but only approximate for a molecule in a periodic
//! padded box, because the periodic images break exact rotational invariance.
//! The reference makes the same choice and this crate documents it rather than
//! hiding it.  A crystal must use [`DisplacementMode::Crystal`]; no rotations
//! are then removed.
//!
//! Modes are returned as Cartesian displacement vectors normalised so that the
//! largest atomic displacement in each mode is one Angstrom.  The mass-weighted
//! eigenvectors are not exposed; the thermochemistry and IR code reconstruct
//! the normal-coordinate scaling from the mass-weighted norm when it needs it.

use serde::{Deserialize, Serialize};
use xtalforge_core::{Structure, element_info};

use crate::Error;
use crate::math::{norm, symmetric_eigen};

use super::constants::CM_PER_SQRT_EV_PER_A2_AMU;
use super::hessian::{Hessian, HessianOptions};
use super::rigid::{project_out, rigid_basis};

/// A frequency this close to zero is a zero mode, in cm^-1.
///
/// The value is the reference implementation's `ZERO_BELOW`.
pub const ZERO_MODE_TOLERANCE_CM: f64 = 5.0;

/// A mode below this frequency is imaginary and counted as one, in cm^-1.
///
/// The value is the reference implementation's `IMAGINARY_BELOW`.
pub const IMAGINARY_BELOW_CM: f64 = -5.0;

/// Numerical round-off of a projected rigid motion is clamped to zero.
const RIGID_FREQUENCY_CM: f64 = 1.0e-6;

/// Frequencies and Cartesian normal modes of one P1 cell.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct VibrationalModes {
    /// Harmonic frequencies in cm^-1, ascending; an imaginary mode is negative.
    pub frequencies: Vec<f64>,
    /// One Cartesian displacement field per mode, `modes[mode][atom][axis]`.
    ///
    /// Each mode is normalised so its largest atomic displacement is one
    /// Angstrom.  A zero or imaginary mode still has a well-defined direction;
    /// the caller decides whether to draw it.
    pub modes: Vec<Vec<[f64; 3]>>,
    /// Indices of modes within [`ZERO_MODE_TOLERANCE_CM`] of zero.
    pub zero_modes: Vec<usize>,
    /// Indices of modes below [`IMAGINARY_BELOW_CM`].
    pub imaginary: Vec<usize>,
}

impl VibrationalModes {
    /// Mass-weights and diagonalizes a Hessian.
    ///
    /// A convenience wrapper for [`vibrational_modes`].
    ///
    /// # Errors
    ///
    /// Returns the errors of [`vibrational_modes`].
    pub fn from_hessian(
        hessian: &Hessian,
        structure: &Structure,
        options: &HessianOptions,
    ) -> Result<Self, Error> {
        vibrational_modes(hessian, structure, options)
    }

    /// Number of modes (`3 N`).
    #[must_use]
    pub fn n_modes(&self) -> usize {
        self.frequencies.len()
    }

    /// Number of P1 atoms, or zero when there are no modes.
    #[must_use]
    pub fn n_atoms(&self) -> usize {
        self.modes
            .first()
            .map(Vec::len)
            .unwrap_or_else(|| self.frequencies.len() / 3)
    }

    /// Whether any mode has an imaginary frequency.
    #[must_use]
    pub fn has_imaginary(&self) -> bool {
        !self.imaginary.is_empty()
    }

    /// The real, non-zero modes as `(index, frequency)` pairs.
    pub fn real_modes(&self) -> impl Iterator<Item = (usize, f64)> + '_ {
        self.frequencies
            .iter()
            .enumerate()
            .filter(|(_, frequency)| **frequency > ZERO_MODE_TOLERANCE_CM)
            .map(|(index, frequency)| (index, *frequency))
    }
}

/// Mass-weights and diagonalizes a Gamma-point Hessian.
///
/// `structure` must be P1 and have the same atoms as the Hessian.  The
/// `options` argument supplies which rigid motions are projected.
///
/// # Errors
///
/// Returns [`Error::InvalidInput`] for a non-P1 structure, a Hessian whose
/// dimensions do not match, a dummy atom, a non-finite matrix element, or a
/// Jacobi diagonalization that does not converge.
pub fn vibrational_modes(
    hessian: &Hessian,
    structure: &Structure,
    options: &HessianOptions,
) -> Result<VibrationalModes, Error> {
    if !structure.is_p1() {
        return Err(Error::InvalidInput(
            "normal modes need a P1 structure; expand the asymmetric unit first".to_owned(),
        ));
    }
    let n_atoms = structure.sites().len();
    if hessian.n_atoms != n_atoms {
        return Err(Error::InvalidInput(format!(
            "Hessian has {} atoms but the structure has {n_atoms}",
            hessian.n_atoms
        )));
    }
    if structure.p1()?.n_atoms() != n_atoms {
        return Err(Error::InvalidInput(
            "the structure expands to a different P1 atom count; build modes from an expanded P1 cell"
                .to_owned(),
        ));
    }
    let size = 3 * n_atoms;
    if hessian.matrix.len() != size || hessian.matrix.iter().any(|row| row.len() != size) {
        return Err(Error::InvalidInput(format!(
            "Hessian must be a square {size} x {size} matrix"
        )));
    }
    if n_atoms == 0 {
        return Err(Error::InvalidInput(
            "cannot build modes for a structure with no atoms".to_owned(),
        ));
    }

    let mut masses = Vec::with_capacity(n_atoms);
    for (index, site) in structure.sites().iter().enumerate() {
        let element = element_info(site.element())?;
        if element.is_dummy() {
            return Err(Error::InvalidInput(format!(
                "atom {index} is the dummy species {}; a dummy atom has no mass and no vibration",
                site.element()
            )));
        }
        masses.push(element.mass);
    }

    let mut weighted = vec![vec![0.0; size]; size];
    for row in 0..size {
        for column in 0..size {
            let value = hessian
                .matrix
                .get(row)
                .and_then(|values| values.get(column))
                .copied()
                .unwrap_or(0.0);
            if !value.is_finite() {
                return Err(Error::InvalidInput(format!(
                    "Hessian element ({row}, {column}) is not finite"
                )));
            }
            let first_mass = masses
                .get(row / 3)
                .copied()
                .ok_or_else(|| Error::InvalidInput("mass index out of range".to_owned()))?;
            let second_mass = masses
                .get(column / 3)
                .copied()
                .ok_or_else(|| Error::InvalidInput("mass index out of range".to_owned()))?;
            let denominator = (first_mass * second_mass).sqrt();
            if denominator <= 0.0 {
                return Err(Error::InvalidInput(
                    "a Hessian mass must be positive".to_owned(),
                ));
            }
            let entry = weighted
                .get_mut(row)
                .and_then(|values| values.get_mut(column))
                .ok_or_else(|| Error::InvalidInput("Hessian index out of range".to_owned()))?;
            *entry = value / denominator;
        }
    }

    let basis = rigid_basis(structure, &masses, options.displacement)?;
    project_out(&mut weighted, &basis);
    symmetrize(&mut weighted);

    let (eigenvalues, eigenvectors) = symmetric_eigen(&weighted).ok_or_else(|| {
        Error::InvalidInput(
            "the mass-weighted Hessian could not be diagonalized; check that it is finite and symmetric"
                .to_owned(),
        )
    })?;

    let frequencies = eigenvalues
        .iter()
        .map(|value| {
            let frequency = value.signum()
                * (value.abs() * CM_PER_SQRT_EV_PER_A2_AMU * CM_PER_SQRT_EV_PER_A2_AMU).sqrt();
            if frequency.abs() < RIGID_FREQUENCY_CM {
                0.0
            } else {
                frequency
            }
        })
        .collect::<Vec<_>>();

    let mut modes = Vec::with_capacity(frequencies.len());
    for vector in &eigenvectors {
        let mut cartesian = vec![[0.0; 3]; n_atoms];
        for index in 0..size {
            let mass = masses
                .get(index / 3)
                .copied()
                .ok_or_else(|| Error::InvalidInput("mass index out of range".to_owned()))?;
            let component =
                vector.get(index).copied().ok_or_else(|| {
                    Error::InvalidInput("eigenvector index out of range".to_owned())
                })? / mass.sqrt();
            let atom = index / 3;
            let axis = index % 3;
            if let Some(point) = cartesian.get_mut(atom)
                && let Some(slot) = point.get_mut(axis)
            {
                *slot = component;
            }
        }
        normalise_mode(&mut cartesian);
        modes.push(cartesian);
    }

    let zero_modes = frequencies
        .iter()
        .enumerate()
        .filter(|(_, frequency)| frequency.abs() <= ZERO_MODE_TOLERANCE_CM)
        .map(|(index, _)| index)
        .collect();
    let imaginary = frequencies
        .iter()
        .enumerate()
        .filter(|(_, frequency)| **frequency < IMAGINARY_BELOW_CM)
        .map(|(index, _)| index)
        .collect();

    Ok(VibrationalModes {
        frequencies,
        modes,
        zero_modes,
        imaginary,
    })
}

fn symmetrize(matrix: &mut [Vec<f64>]) {
    let n = matrix.len();
    for row in 0..n {
        for column in (row + 1)..n {
            let first = matrix
                .get(row)
                .and_then(|values| values.get(column))
                .copied()
                .unwrap_or(0.0);
            let second = matrix
                .get(column)
                .and_then(|values| values.get(row))
                .copied()
                .unwrap_or(0.0);
            let mean = 0.5 * (first + second);
            if let Some(values) = matrix.get_mut(row)
                && let Some(slot) = values.get_mut(column)
            {
                *slot = mean;
            }
            if let Some(values) = matrix.get_mut(column)
                && let Some(slot) = values.get_mut(row)
            {
                *slot = mean;
            }
        }
    }
}

fn normalise_mode(mode: &mut [[f64; 3]]) {
    let mut largest: f64 = 0.0;
    for point in mode.iter() {
        largest = largest.max(norm(*point));
    }
    if largest <= 1.0e-14 {
        return;
    }
    for point in mode.iter_mut() {
        for value in point.iter_mut() {
            *value /= largest;
        }
    }
}
