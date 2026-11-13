//! Coordinate constraints by a quadratic penalty (Task 4.5.2).
//!
//! A relaxed scan holds one coordinate while everything else moves.  This
//! crate's core implementation adds a quadratic penalty,
//!
//! ```text
//! E_c = 1/2 k (value - target)^2,
//! ```
//!
//! and its gradient to the calculator's energy and force.  That is the
//! "penalty" half of Task 4.5.2; the exact holonomic projection of
//! Crystal-Builder `xtal/ff/constraints.py` (project the search direction and
//! restore the coordinate by Newton steps) is a larger follow-up.  The penalty
//! is documented here so no caller mistakes it for an exact constraint: the
//! reported `achieved` value is the number to inspect.

use serde::{Deserialize, Serialize};
use xtalforge_core::Lattice;

use crate::Error;
use crate::math::scale;
use crate::scan::coordinates::ScanCoordinate;

/// One coordinate held near a target value.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct Constraint {
    /// Coordinate being held.
    pub coordinate: ScanCoordinate,
    /// Target value in the coordinate's units.
    pub target: f64,
}

impl Constraint {
    /// A constraint for any coordinate.
    ///
    /// # Errors
    ///
    /// Returns [`Error::InvalidInput`] when the target is not finite, or when
    /// the coordinate is a cell parameter: a Cartesian penalty cannot hold a
    /// lattice parameter, so cell parameters are imposed by a scan axis
    /// instead.
    pub fn new(coordinate: ScanCoordinate, target: f64) -> Result<Self, Error> {
        if !target.is_finite() {
            return Err(Error::InvalidInput(
                "constraint target must be finite".to_owned(),
            ));
        }
        if coordinate.is_cell_parameter() {
            return Err(Error::InvalidInput(format!(
                "{} is a cell parameter and cannot be held by a Cartesian penalty; \
                 scan it with a fixed cell instead",
                coordinate.label()
            )));
        }
        Ok(Self { coordinate, target })
    }

    /// A distance constraint in Angstrom.
    ///
    /// # Errors
    ///
    /// Returns the errors of [`ScanCoordinate::distance`] or
    /// [`Constraint::new`].
    pub fn distance(i: usize, j: usize, target: f64) -> Result<Self, Error> {
        Self::new(ScanCoordinate::distance(i, j)?, target)
    }

    /// An angle constraint in degrees.
    ///
    /// # Errors
    ///
    /// Returns the errors of [`ScanCoordinate::angle`] or
    /// [`Constraint::new`].
    pub fn angle(i: usize, j: usize, k: usize, target: f64) -> Result<Self, Error> {
        Self::new(ScanCoordinate::angle(i, j, k)?, target)
    }

    /// A torsion constraint in degrees.
    ///
    /// # Errors
    ///
    /// Returns the errors of [`ScanCoordinate::torsion`] or
    /// [`Constraint::new`].
    pub fn torsion(i: usize, j: usize, k: usize, l: usize, target: f64) -> Result<Self, Error> {
        Self::new(ScanCoordinate::torsion(i, j, k, l)?, target)
    }

    /// Current value of the coordinate.
    ///
    /// # Errors
    ///
    /// Returns the errors of [`ScanCoordinate::value`].
    pub fn value(&self, positions: &[[f64; 3]], lattice: &Lattice) -> Result<f64, Error> {
        self.coordinate.value(positions, lattice)
    }

    /// Signed residual wrapped onto the nearest turn for angular coordinates.
    ///
    /// # Errors
    ///
    /// Returns the errors of [`ScanCoordinate::value`].
    pub fn residual(&self, positions: &[[f64; 3]], lattice: &Lattice) -> Result<f64, Error> {
        let mut residual = self.value(positions, lattice)? - self.target;
        if self.coordinate.units() == "deg" {
            residual -= 360.0 * (residual / 360.0).round();
        }
        Ok(residual)
    }

    /// Penalty energy and its gradient with respect to Cartesian positions.
    ///
    /// `constant` is `k` in eV per squared coordinate unit.  The returned
    /// energy is `1/2 k residual^2` and the gradient is `k residual d(value)`.
    ///
    /// # Errors
    ///
    /// Returns [`Error::InvalidInput`] for a negative or non-finite constant,
    /// for a cell-parameter coordinate, and the errors of
    /// [`ScanCoordinate::gradient`].
    pub fn penalty(
        &self,
        positions: &[[f64; 3]],
        lattice: &Lattice,
        constant: f64,
    ) -> Result<(f64, Vec<[f64; 3]>), Error> {
        if !constant.is_finite() || constant < 0.0 {
            return Err(Error::InvalidInput(format!(
                "penalty constant must be non-negative and finite, got {constant}"
            )));
        }
        if self.coordinate.is_cell_parameter() {
            return Err(Error::InvalidInput(format!(
                "{} is a cell parameter and has a zero Cartesian gradient",
                self.coordinate.label()
            )));
        }
        let residual = self.residual(positions, lattice)?;
        let gradient = self.coordinate.gradient(positions, lattice)?;
        let energy = 0.5 * constant * residual * residual;
        let scaled: Vec<[f64; 3]> = gradient
            .iter()
            .map(|value| scale(*value, constant * residual))
            .collect();
        Ok((energy, scaled))
    }
}
