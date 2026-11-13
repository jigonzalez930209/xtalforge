//! Undoable unit-cell edits (Task 3.10.3).

use serde::{Deserialize, Serialize};
use xtalforge_core::{Error, Lattice, Slab, Structure, make_slab};

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

/// Re-expresses the crystal on a general integer basis (Task 14.3.3).
///
/// The rows of the matrix are the new lattice vectors in units of the old
/// ones, so a determinant of one is a change of basis and a larger determinant
/// is a supercell no `na x nb x nc` repeat can express. The result is P1 and
/// its bonds are perceived afresh, exactly like [`Supercell`].
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct TransformCell {
    matrix: [[i32; 3]; 3],
}

impl TransformCell {
    /// Creates the command from an integer transformation matrix.
    #[must_use]
    pub fn new(matrix: [[i32; 3]; 3]) -> Self {
        Self { matrix }
    }

    /// Transformation matrix.
    #[must_use]
    pub fn matrix(&self) -> [[i32; 3]; 3] {
        self.matrix
    }

    /// Determinant of the transformation matrix.
    #[must_use]
    pub fn determinant(&self) -> i32 {
        let [[a, b, c], [d, e, f], [g, h, i]] = self.matrix;
        a * (e * i - f * h) - b * (d * i - f * g) + c * (d * h - e * g)
    }
}

impl Command for TransformCell {
    fn apply(&self, structure: &Structure) -> Result<Structure, Error> {
        let p = self.matrix.map(|row| row.map(f64::from));
        structure.supercell_matrix(p)
    }

    fn description(&self) -> String {
        format!("Transform cell (det P = {})", self.determinant())
    }
}

/// Builds a slab cut along `(hkl)` (Task 14.3.3).
///
/// The parent bonds are carried to the copies they join and the bonds the two
/// surfaces sever are counted in [`Slab::cut`]; see
/// [`xtalforge_core::make_slab`].
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct MakeSlab {
    hkl: [i32; 3],
    layers: u32,
    vacuum: f64,
    shift: f64,
}

impl MakeSlab {
    /// Creates the command.
    #[must_use]
    pub fn new(hkl: [i32; 3], layers: u32, vacuum: f64, shift: f64) -> Self {
        Self {
            hkl,
            layers,
            vacuum,
            shift,
        }
    }

    /// Creates the command with one layer, 15 Angstrom of vacuum and no shift.
    #[must_use]
    pub fn with_defaults(hkl: [i32; 3]) -> Self {
        Self::new(hkl, 1, 15.0, 0.0)
    }

    /// Miller index of the plane.
    #[must_use]
    pub fn hkl(&self) -> [i32; 3] {
        self.hkl
    }

    /// Number of layers.
    #[must_use]
    pub fn layers(&self) -> u32 {
        self.layers
    }

    /// Vacuum above the slab in Angstrom.
    #[must_use]
    pub fn vacuum(&self) -> f64 {
        self.vacuum
    }

    /// Termination shift as a fraction of one layer.
    #[must_use]
    pub fn shift(&self) -> f64 {
        self.shift
    }

    /// Builds the slab, including the severed-bond count.
    ///
    /// # Errors
    ///
    /// Returns the errors of [`xtalforge_core::make_slab`].
    pub fn make(&self, structure: &Structure) -> Result<Slab, Error> {
        make_slab(structure, self.hkl, self.layers, self.vacuum, self.shift)
    }
}

impl Command for MakeSlab {
    fn apply(&self, structure: &Structure) -> Result<Structure, Error> {
        Ok(self.make(structure)?.structure)
    }

    fn description(&self) -> String {
        format!("Slab ({} {} {})", self.hkl[0], self.hkl[1], self.hkl[2])
    }
}
