//! Cell parameters as scan coordinates (Task 14.2.3).
//!
//! A cell-parameter scan varies one of the six parameters and lets the sites
//! relax in the resulting cell.  This is the workflow behind an equation-of-
//! state or pressure sweep:
//!
//! ```text
//! for a in 5.0..=6.0 { set a; relax the sites at fixed a; record E(a) }
//! ```
//!
//! The parameter is imposed on the seeded lattice *before* the relaxation and,
//! when the optimizer also relaxes the cell, held there through
//! [`CellFreedom`](crate::CellFreedom) by the scan driver.  Every fractional
//! coordinate is kept when the parameter changes, so the atoms follow the cell
//! (the affine deformation an equation of state needs).  Imposing a parameter
//! rebuilds the cell in its standard orientation; a rotated input cell is
//! therefore reoriented, which changes no energy or distance.

use serde::{Deserialize, Serialize};
use xtalforge_core::Lattice;

use crate::Error;
use crate::optimize::CellFreedom;

/// One of the six cell parameters.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash, Serialize, Deserialize)]
#[serde(rename_all = "snake_case")]
pub enum CellParameter {
    /// Edge `a` in Angstrom.
    A,
    /// Edge `b` in Angstrom.
    B,
    /// Edge `c` in Angstrom.
    C,
    /// Angle `alpha` in degrees.
    Alpha,
    /// Angle `beta` in degrees.
    Beta,
    /// Angle `gamma` in degrees.
    Gamma,
}

impl CellParameter {
    /// Every parameter, in the order of [`Lattice::parameters`].
    pub const ALL: [Self; 6] = [
        Self::A,
        Self::B,
        Self::C,
        Self::Alpha,
        Self::Beta,
        Self::Gamma,
    ];

    /// Index in [`Lattice::parameters`], `0..6`.
    #[must_use]
    pub const fn index(self) -> usize {
        match self {
            Self::A => 0,
            Self::B => 1,
            Self::C => 2,
            Self::Alpha => 3,
            Self::Beta => 4,
            Self::Gamma => 5,
        }
    }

    /// Short name: `"a"`, `"b"`, `"c"`, `"alpha"`, `"beta"` or `"gamma"`.
    #[must_use]
    pub const fn name(self) -> &'static str {
        match self {
            Self::A => "a",
            Self::B => "b",
            Self::C => "c",
            Self::Alpha => "alpha",
            Self::Beta => "beta",
            Self::Gamma => "gamma",
        }
    }

    /// Units of the parameter: `"A"` for a length, `"deg"` for an angle.
    #[must_use]
    pub const fn units(self) -> &'static str {
        match self {
            Self::A | Self::B | Self::C => "A",
            Self::Alpha | Self::Beta | Self::Gamma => "deg",
        }
    }

    /// Current value of the parameter in `lattice`.
    #[must_use]
    pub fn value(self, lattice: &Lattice) -> f64 {
        lattice
            .parameters()
            .get(self.index())
            .copied()
            .unwrap_or(0.0)
    }

    /// Cell rebuilt with this parameter replaced, in standard orientation.
    ///
    /// # Errors
    ///
    /// Returns [`Error::InvalidInput`] for a non-finite value, a non-positive
    /// length, an angle outside `(0, 180)` degrees, or a value that makes the
    /// cell geometrically impossible.
    pub fn with_value(self, lattice: &Lattice, value: f64) -> Result<Lattice, Error> {
        if !value.is_finite() {
            return Err(Error::InvalidInput(format!(
                "cell parameter {} must be finite, got {value}",
                self.name()
            )));
        }
        let current = lattice.parameters();
        let mut next = [0.0_f64; 6];
        for (index, slot) in next.iter_mut().enumerate() {
            *slot = current.get(index).copied().unwrap_or(0.0);
        }
        if let Some(slot) = next.get_mut(self.index()) {
            *slot = value;
        }
        let [a, b, c, alpha, beta, gamma] = next;
        Ok(Lattice::from_parameters(a, b, c, alpha, beta, gamma)?)
    }

    /// Adds this parameter to the held quantities of `freedom`.
    ///
    /// The scan driver uses it so a variable-cell relaxation cannot move the
    /// parameter that is being scanned.
    #[must_use]
    pub fn hold(self, freedom: CellFreedom) -> CellFreedom {
        match self {
            Self::A => freedom.hold_length(0),
            Self::B => freedom.hold_length(1),
            Self::C => freedom.hold_length(2),
            Self::Alpha => freedom.hold_angle(0),
            Self::Beta => freedom.hold_angle(1),
            Self::Gamma => freedom.hold_angle(2),
        }
    }
}
