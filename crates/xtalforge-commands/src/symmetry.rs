//! Undoable symmetry operations (Task 3.10.4).

use xtalforge_core::{Error, Structure};
use xtalforge_symmetry::{
    SetMode, SpaceGroup, SymmetryOptions, asymmetrize, primitive_cell, set_space_group,
    standardize, to_p1,
};

use crate::command::Command;

/// Expands every orbit and discards the space group.
#[derive(Debug, Clone, Copy)]
pub struct ReduceToP1;

impl Command for ReduceToP1 {
    fn apply(&self, structure: &Structure) -> Result<Structure, Error> {
        to_p1(structure)
    }

    fn description(&self) -> String {
        "Reduce to P1".to_owned()
    }
}

/// Detects the space group and keeps the resulting asymmetric unit.
#[derive(Debug, Clone, Copy, Default)]
pub struct FindSymmetry {
    options: SymmetryOptions,
}

impl FindSymmetry {
    /// Creates the command with the default detection tolerance.
    #[must_use]
    pub fn new() -> Self {
        Self::default()
    }

    /// Creates the command with explicit options.
    #[must_use]
    pub fn with_options(options: SymmetryOptions) -> Self {
        Self { options }
    }

    /// Creates the command with a position tolerance in Angstrom.
    #[must_use]
    pub fn with_tolerance(symprec: f64) -> Self {
        Self {
            options: SymmetryOptions {
                symprec,
                ..SymmetryOptions::default()
            },
        }
    }

    /// Detection options.
    #[must_use]
    pub fn options(&self) -> SymmetryOptions {
        self.options
    }
}

impl Command for FindSymmetry {
    fn apply(&self, structure: &Structure) -> Result<Structure, Error> {
        asymmetrize(structure, self.options)
    }

    fn description(&self) -> String {
        "Find symmetry".to_owned()
    }
}

/// Reinterprets or imposes a space group.
#[derive(Debug, Clone)]
pub struct SetSpaceGroup {
    group: SpaceGroup,
    mode: SetMode,
}

impl SetSpaceGroup {
    /// Creates the command from a resolved space group.
    #[must_use]
    pub fn new(group: SpaceGroup, mode: SetMode) -> Self {
        Self { group, mode }
    }

    /// Creates the command for the default setting of an ITA number.
    ///
    /// # Errors
    ///
    /// Returns [`Error::Symmetry`] when the number is outside 1 to 230.
    pub fn from_number(number: u32, mode: SetMode) -> Result<Self, Error> {
        Ok(Self::new(SpaceGroup::from_number(number)?, mode))
    }

    /// Space group the command installs.
    #[must_use]
    pub fn group(&self) -> &SpaceGroup {
        &self.group
    }

    /// Application mode.
    #[must_use]
    pub fn mode(&self) -> SetMode {
        self.mode
    }
}

impl Command for SetSpaceGroup {
    fn apply(&self, structure: &Structure) -> Result<Structure, Error> {
        Ok(set_space_group(structure, &self.group, self.mode)?.structure)
    }

    fn description(&self) -> String {
        format!("Set space group {}", self.group.hm_short)
    }
}

/// Standardizes the cell in the conventional or primitive setting.
#[derive(Debug, Clone, Copy, Default)]
pub struct Standardize {
    options: SymmetryOptions,
    primitive: bool,
}

impl Standardize {
    /// Creates the command with the default tolerance and conventional cell.
    #[must_use]
    pub fn new() -> Self {
        Self::default()
    }

    /// Creates the command with explicit detection options.
    #[must_use]
    pub fn with_options(options: SymmetryOptions) -> Self {
        Self {
            options,
            primitive: false,
        }
    }

    /// Selects the primitive instead of the conventional cell.
    #[must_use]
    pub fn with_primitive(mut self, primitive: bool) -> Self {
        self.primitive = primitive;
        self
    }

    /// Detection options.
    #[must_use]
    pub fn options(&self) -> SymmetryOptions {
        self.options
    }

    /// Whether the primitive cell is requested.
    #[must_use]
    pub fn primitive(&self) -> bool {
        self.primitive
    }
}

impl Command for Standardize {
    fn apply(&self, structure: &Structure) -> Result<Structure, Error> {
        let result = if self.primitive {
            primitive_cell(structure, self.options)?
        } else {
            standardize(structure, self.options)?
        };
        Ok(result.structure)
    }

    fn description(&self) -> String {
        if self.primitive {
            "Standardize to primitive cell".to_owned()
        } else {
            "Standardize cell".to_owned()
        }
    }
}
