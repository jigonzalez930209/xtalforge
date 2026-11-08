//! Undoable unit-cell edits (Task 3.10.3).

use serde::{Deserialize, Serialize};
use xtalforge_core::{Error, Lattice, Structure};

use crate::command::Command;
use crate::edit::rebuild_with;

/// What stays fixed when the unit cell is replaced.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
#[serde(rename_all = "snake_case")]
pub enum PositionPolicy {
    /// Keep fractional coordinates; atoms move with the cell.
    KeepFractional,
    /// Keep Cartesian positions; fractional coordinates rescale.
    KeepCartesian,
}

/// Replaces the unit cell under one of the position policies.
#[derive(Debug, Clone)]
pub struct SetLattice {
    lattice: Lattice,
    keep: PositionPolicy,
}

impl SetLattice {
    /// Creates the command.
    #[must_use]
    pub fn new(lattice: Lattice, keep: PositionPolicy) -> Self {
        Self { lattice, keep }
    }

    /// Cell the command installs.
    #[must_use]
    pub fn lattice(&self) -> &Lattice {
        &self.lattice
    }

    /// Position policy.
    #[must_use]
    pub fn keep(&self) -> PositionPolicy {
        self.keep
    }
}

impl Command for SetLattice {
    fn apply(&self, structure: &Structure) -> Result<Structure, Error> {
        let sites = match self.keep {
            PositionPolicy::KeepFractional => structure.sites().to_vec(),
            PositionPolicy::KeepCartesian => {
                let mut sites = Vec::with_capacity(structure.sites().len());
                for site in structure.sites() {
                    let cart = structure.lattice().to_cart(site.frac());
                    let frac = self.lattice.to_frac(cart);
                    sites.push(site.clone().with_frac(frac)?);
                }
                sites
            }
        };
        rebuild_with(
            structure,
            sites,
            self.lattice,
            structure.hall(),
            structure.operations(),
        )
    }

    fn description(&self) -> String {
        match self.keep {
            PositionPolicy::KeepFractional => "Set lattice (keep fractional)".to_owned(),
            PositionPolicy::KeepCartesian => "Set lattice (keep Cartesian)".to_owned(),
        }
    }
}

/// Builds an `na x nb x nc` supercell in P1.
#[derive(Debug, Clone, Copy)]
pub struct Supercell {
    na: u32,
    nb: u32,
    nc: u32,
}

impl Supercell {
    /// Creates the command.
    #[must_use]
    pub fn new(na: u32, nb: u32, nc: u32) -> Self {
        Self { na, nb, nc }
    }

    /// Repeat counts along a, b and c.
    #[must_use]
    pub fn counts(&self) -> (u32, u32, u32) {
        (self.na, self.nb, self.nc)
    }
}

impl Command for Supercell {
    fn apply(&self, structure: &Structure) -> Result<Structure, Error> {
        structure.supercell(self.na, self.nb, self.nc)
    }

    fn description(&self) -> String {
        format!("Supercell {}x{}x{}", self.na, self.nb, self.nc)
    }
}

/// Shifts the origin of a P1 structure and wraps the sites into the cell.
#[derive(Debug, Clone, Copy)]
pub struct ShiftOrigin {
    delta: [f64; 3],
}

impl ShiftOrigin {
    /// Creates the command from a fractional shift.
    #[must_use]
    pub fn new(delta: [f64; 3]) -> Self {
        Self { delta }
    }

    /// Fractional shift.
    #[must_use]
    pub fn delta(&self) -> [f64; 3] {
        self.delta
    }
}

impl Command for ShiftOrigin {
    fn apply(&self, structure: &Structure) -> Result<Structure, Error> {
        if !self.delta.iter().all(|value| value.is_finite()) {
            return Err(Error::InvalidInput(
                "origin shift must be finite".to_owned(),
            ));
        }
        if !structure.is_p1() {
            return Err(Error::InvalidInput(
                "moving the origin requires a P1 structure; reduce to P1 first".to_owned(),
            ));
        }
        let mut sites = Vec::with_capacity(structure.sites().len());
        for site in structure.sites() {
            let frac = [
                site.frac()[0] + self.delta[0],
                site.frac()[1] + self.delta[1],
                site.frac()[2] + self.delta[2],
            ];
            sites.push(site.clone().with_frac(Lattice::wrap(frac))?);
        }
        rebuild_with(
            structure,
            sites,
            *structure.lattice(),
            structure.hall(),
            structure.operations(),
        )
    }

    fn description(&self) -> String {
        format!(
            "Shift origin by ({}, {}, {})",
            self.delta[0], self.delta[1], self.delta[2]
        )
    }
}

/// Folds every site back into `[0, 1)`.
#[derive(Debug, Clone, Copy)]
pub struct WrapIntoCell;

impl Command for WrapIntoCell {
    fn apply(&self, structure: &Structure) -> Result<Structure, Error> {
        let mut sites = Vec::with_capacity(structure.sites().len());
        for site in structure.sites() {
            sites.push(site.clone().with_frac(Lattice::wrap(site.frac()))?);
        }
        rebuild_with(
            structure,
            sites,
            *structure.lattice(),
            structure.hall(),
            structure.operations(),
        )
    }

    fn description(&self) -> String {
        "Wrap into the cell".to_owned()
    }
}
