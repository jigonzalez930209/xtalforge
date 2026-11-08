//! Gamma-point Hessian by central differences of analytic forces (Task 4.8.1).
//!
//! Every P1 atom is displaced by `+step` and `-step` along each Cartesian
//! axis and the force is evaluated at both points; one column of the Hessian is
//! the central difference of the forces.  The implementation needs
//! `6 N + 1` calculator evaluations and no analytic second derivatives, so it
//! works with UFF today and with any external engine that returns forces.
//!
//! Units are the crate boundary: displacements in Angstrom, forces in
//! eV/Angstrom, and the Hessian in eV/Angstrom^2.  The result is the physical
//! (mass-independent) Hessian; mass weighting and the rigid-mode projection
//! belong to [`crate::vibrations::vibrational_modes`].
//!
//! **Documented approximations.**  The calculation is at the Gamma point of a
//! P1 cell.  A structure with space-group operations must be expanded by the
//! caller first.  Dummy `X` atoms are refused because they carry no mass and
//! no physical vibration.  A finite difference is only as good as the chosen
//! step: the default is 0.01 Angstrom, the value phonopy and ASE use, which is
//! small enough for a force field's quadratic region and large enough to stay
//! above double-precision noise.

use serde::{Deserialize, Serialize};
use xtalforge_core::{Structure, element_info};

use crate::Error;
use crate::calculator::Calculator;
use crate::structure_util::with_frac;

/// Default central-difference displacement, in Angstrom.
pub const DEFAULT_HESSIAN_STEP: f64 = 0.01;

/// Which rigid motions are removed from the mass-weighted Hessian.
///
/// A periodic crystal has three uniform translations as exact zero modes (the
/// acoustic sum rule); an isolated molecule has three translations and, for a
/// non-linear molecule, three rotations.  Projecting them out avoids reporting
/// finite-difference round-off as a soft vibration.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default, Serialize, Deserialize)]
pub enum DisplacementMode {
    /// A periodic crystal: project the three translations.
    #[default]
    Crystal,
    /// An isolated molecule: project translations and rotations.
    Molecule,
}

/// Options for one finite-difference Hessian.
#[derive(Debug, Clone, Copy, PartialEq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct HessianOptions {
    /// Central-difference displacement in Angstrom.
    pub step: f64,
    /// Rigid motions the mode analysis removes.
    pub displacement: DisplacementMode,
}

impl Default for HessianOptions {
    fn default() -> Self {
        Self::crystal()
    }
}

impl HessianOptions {
    /// Options for a periodic crystal.
    #[must_use]
    pub fn crystal() -> Self {
        Self {
            step: DEFAULT_HESSIAN_STEP,
            displacement: DisplacementMode::Crystal,
        }
    }

    /// Options for an isolated molecule.
    #[must_use]
    pub fn molecule() -> Self {
        Self {
            step: DEFAULT_HESSIAN_STEP,
            displacement: DisplacementMode::Molecule,
        }
    }

    /// Options with an explicit step, keeping the crystal projection.
    ///
    /// # Errors
    ///
    /// Returns [`Error::InvalidInput`] when the step is not positive and
    /// finite.
    pub fn with_step(step: f64) -> Result<Self, Error> {
        let options = Self {
            step,
            ..Self::crystal()
        };
        options.validate()?;
        Ok(options)
    }

    /// Checks that the options are usable.
    ///
    /// # Errors
    ///
    /// Returns [`Error::InvalidInput`] for a non-positive or non-finite step.
    pub fn validate(&self) -> Result<(), Error> {
        if !self.step.is_finite() || self.step <= 0.0 {
            return Err(Error::InvalidInput(format!(
                "Hessian step must be positive and finite, got {}",
                self.step
            )));
        }
        Ok(())
    }
}

/// A physical Hessian at the Gamma point.
///
/// `matrix[row][column]` is `d^2 E / d x_row d x_column` in
/// eV/Angstrom^2, with row and column `3 * atom + axis` over the P1 atoms.
/// The matrix is symmetric to rounding.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct Hessian {
    /// Row-major `3 N x 3 N` matrix in eV/Angstrom^2.
    pub matrix: Vec<Vec<f64>>,
    /// Number of P1 atoms.
    pub n_atoms: usize,
}

impl Hessian {
    /// Number of rows and columns (`3 * n_atoms`).
    #[must_use]
    pub fn size(&self) -> usize {
        3 * self.n_atoms
    }

    /// One matrix element, or `None` when either index is out of range.
    #[must_use]
    pub fn value(&self, row: usize, column: usize) -> Option<f64> {
        self.matrix
            .get(row)
            .and_then(|values| values.get(column))
            .copied()
    }
}

/// Builds the Gamma-point Hessian of `structure`.
///
/// The structure must be P1: the central difference moves one atom at a time
/// and a chain of symmetry operations would otherwise duplicate the
/// displacement.  Expand the asymmetric unit first with
/// [`xtalforge_core::Structure::p1`] and rebuild it, or use a P1 input.
///
/// # Errors
///
/// Returns [`Error::InvalidInput`] for invalid options, a non-P1 or empty
/// structure, a dummy atom, a force vector of the wrong length, or a
/// non-finite energy or force.  Returns the calculator's own error when an
/// evaluation fails.
pub fn hessian(
    calculator: &dyn Calculator,
    structure: &Structure,
    options: &HessianOptions,
) -> Result<Hessian, Error> {
    options.validate()?;
    if !calculator.provides_forces() {
        return Err(Error::InvalidInput(format!(
            "calculator {} does not provide forces, so a finite-difference Hessian is impossible",
            calculator.name()
        )));
    }
    if !structure.is_p1() {
        return Err(Error::InvalidInput(
            "a Gamma-point Hessian needs a P1 structure; expand the asymmetric unit first"
                .to_owned(),
        ));
    }
    let n_atoms = structure.sites().len();
    if n_atoms == 0 {
        return Err(Error::InvalidInput(
            "cannot take the Hessian of a structure with no atoms".to_owned(),
        ));
    }
    if structure.p1()?.n_atoms() != n_atoms {
        return Err(Error::InvalidInput(
            "the structure expands to a different P1 atom count; give the Hessian an already expanded P1 cell"
                .to_owned(),
        ));
    }
    for (index, site) in structure.sites().iter().enumerate() {
        if element_info(site.element())?.is_dummy() {
            return Err(Error::InvalidInput(format!(
                "atom {index} is the dummy species {}; a dummy atom has no mass and no vibration",
                site.element()
            )));
        }
    }

    let lattice = *structure.lattice();
    let cart: Vec<[f64; 3]> = structure
        .sites()
        .iter()
        .map(|site| lattice.to_cart(site.frac()))
        .collect();
    let (_, reference) = calculator.energy_forces(structure)?;
    validate_forces(&reference, n_atoms, "reference geometry")?;

    let size = 3 * n_atoms;
    let mut matrix = vec![vec![0.0; size]; size];
    for atom in 0..n_atoms {
        for axis in 0..3 {
            let plus = displaced(calculator, structure, &cart, atom, axis, options.step)?;
            let minus = displaced(calculator, structure, &cart, atom, axis, -options.step)?;
            for (other, (plus_atom, minus_atom)) in plus.iter().zip(&minus).enumerate() {
                for component in 0..3 {
                    let (Some(plus_force), Some(minus_force)) =
                        (plus_atom.get(component), minus_atom.get(component))
                    else {
                        return Err(Error::InvalidInput(
                            "calculator returned a force vector shorter than three components"
                                .to_owned(),
                        ));
                    };
                    let value = -(plus_force - minus_force) / (2.0 * options.step);
                    if let Some(row) = matrix.get_mut(3 * atom + axis)
                        && let Some(slot) = row.get_mut(3 * other + component)
                    {
                        *slot = value;
                    }
                }
            }
        }
    }

    for row in 0..size {
        for column in (row + 1)..size {
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
            if !mean.is_finite() {
                return Err(Error::InvalidInput(format!(
                    "Hessian element ({row}, {column}) is not finite"
                )));
            }
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
    for (row, values) in matrix.iter().enumerate() {
        for (column, value) in values.iter().enumerate() {
            if !value.is_finite() {
                return Err(Error::InvalidInput(format!(
                    "Hessian element ({row}, {column}) is not finite"
                )));
            }
        }
    }

    Ok(Hessian { matrix, n_atoms })
}

fn displaced(
    calculator: &dyn Calculator,
    structure: &Structure,
    cart: &[[f64; 3]],
    atom: usize,
    axis: usize,
    delta: f64,
) -> Result<Vec<[f64; 3]>, Error> {
    let lattice = *structure.lattice();
    let mut moved = cart.to_vec();
    let Some(point) = moved.get_mut(atom) else {
        return Err(Error::InvalidInput(format!(
            "atom index {atom} is outside the P1 cell"
        )));
    };
    let Some(component) = point.get_mut(axis) else {
        return Err(Error::InvalidInput(format!(
            "Cartesian component {axis} is outside 0..3"
        )));
    };
    *component += delta;
    let frac: Vec<[f64; 3]> = moved
        .iter()
        .map(|position| lattice.to_frac(*position))
        .collect();
    let moved = with_frac(structure, lattice, &frac)?;
    let (_, forces) = calculator.energy_forces(&moved)?;
    validate_forces(&forces, cart.len(), "displaced geometry")?;
    Ok(forces)
}

fn validate_forces(forces: &[[f64; 3]], n_atoms: usize, place: &str) -> Result<(), Error> {
    if forces.len() != n_atoms {
        return Err(Error::InvalidInput(format!(
            "calculator returned {} force vectors for {n_atoms} atoms at the {place}",
            forces.len()
        )));
    }
    for (index, force) in forces.iter().enumerate() {
        if force.iter().any(|value| !value.is_finite()) {
            return Err(Error::InvalidInput(format!(
                "calculator returned a non-finite force at atom {index} of the {place}"
            )));
        }
    }
    Ok(())
}
