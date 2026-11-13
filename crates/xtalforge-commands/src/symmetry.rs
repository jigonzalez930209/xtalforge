//! Undoable symmetry operations (Task 3.10.4).

use std::collections::BTreeMap;

use serde_json::json;
use xtalforge_core::{Error, Site, Structure};
use xtalforge_symmetry::{
    SetMode, SpaceGroup, SymmetryOptions, WyckoffAssignment, asymmetrize, lower_to_subgroup,
    primitive_cell, set_space_group, standardize, to_p1, wyckoff_assignments,
};

use crate::command::Command;
use crate::edit::rebuild_with;

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

/// Keys written by [`AssignWyckoff`].
pub const WYCKOFF_KEYS: [&str; 4] = [
    "wyckoff",
    "wyckoff_letter",
    "wyckoff_multiplicity",
    "site_symmetry",
];

/// Fills in the Wyckoff letter of every site (Task 14.3.3).
///
/// The assignment is informational: no coordinate moves. The letter and
/// multiplicity are stored as site props (`wyckoff` as `4a`,
/// `wyckoff_letter`, `wyckoff_multiplicity` and `site_symmetry`) so they travel
/// with the site through a save and can be cleared by undoing the command.
#[derive(Debug, Clone, Copy, Default)]
pub struct AssignWyckoff {
    options: SymmetryOptions,
}

impl AssignWyckoff {
    /// Creates the command with the default detection tolerance.
    #[must_use]
    pub fn new() -> Self {
        Self::default()
    }

    /// Creates the command with explicit detection options.
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

impl Command for AssignWyckoff {
    fn apply(&self, structure: &Structure) -> Result<Structure, Error> {
        let assignments = wyckoff_assignments(structure, self.options)?;
        let p1 = structure.p1()?;
        let mut by_site: BTreeMap<usize, &WyckoffAssignment> = BTreeMap::new();
        for assignment in &assignments {
            let Some(atom) = p1.atoms.get(assignment.atom) else {
                continue;
            };
            by_site.entry(atom.site_idx).or_insert(assignment);
        }
        let mut sites: Vec<Site> = Vec::with_capacity(structure.sites().len());
        for (index, site) in structure.sites().iter().enumerate() {
            let mut updated = site.clone();
            for key in WYCKOFF_KEYS {
                updated = updated.without_prop(key);
            }
            if let Some(assignment) = by_site.get(&index) {
                updated = updated
                    .with_prop(
                        "wyckoff",
                        json!(format!("{}{}", assignment.multiplicity, assignment.letter)),
                    )
                    .with_prop("wyckoff_letter", json!(assignment.letter.to_string()))
                    .with_prop("wyckoff_multiplicity", json!(assignment.multiplicity))
                    .with_prop("site_symmetry", json!(assignment.site_symmetry));
            }
            sites.push(updated);
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
        "Assign Wyckoff letters".to_owned()
    }
}

/// The same crystal in the other hand (Task 14.3.3).
///
/// A centrosymmetric structure is already its own enantiomer and comes back
/// unchanged. A chiral setting that cannot describe the inverted coordinates is
/// expanded to P1; see [`xtalforge_core::invert()`].
#[derive(Debug, Clone, Copy, Default)]
pub struct Invert;

impl Command for Invert {
    fn apply(&self, structure: &Structure) -> Result<Structure, Error> {
        xtalforge_core::invert(structure)
    }

    fn description(&self) -> String {
        "Invert the structure".to_owned()
    }
}

/// Drops to a subgroup so that an orbit splits (Task 14.3.3).
///
/// The subgroup operations are taken in the default setting of
/// `subgroup_number`, the cell is left alone, and the result is P1.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct DescendToSubgroup {
    subgroup_number: u32,
    tolerance: f64,
}

impl DescendToSubgroup {
    /// Creates the command with a 0.05 Angstrom matching tolerance.
    #[must_use]
    pub fn new(subgroup_number: u32) -> Self {
        Self {
            subgroup_number,
            tolerance: 0.05,
        }
    }

    /// Sets the matching tolerance in Angstrom.
    #[must_use]
    pub fn with_tolerance(mut self, tolerance: f64) -> Self {
        self.tolerance = tolerance;
        self
    }

    /// ITA number of the target subgroup.
    #[must_use]
    pub fn subgroup_number(&self) -> u32 {
        self.subgroup_number
    }

    /// Matching tolerance in Angstrom.
    #[must_use]
    pub fn tolerance(&self) -> f64 {
        self.tolerance
    }
}

impl Command for DescendToSubgroup {
    fn apply(&self, structure: &Structure) -> Result<Structure, Error> {
        lower_to_subgroup(structure, self.subgroup_number, self.tolerance)
    }

    fn description(&self) -> String {
        format!("Descend to subgroup {}", self.subgroup_number)
    }
}

/// Merges same-element sites that sit on the same position (Task 14.3.3).
///
/// The result is P1, because a duplicate is judged over the whole cell; run
/// [`FindSymmetry`] afterwards to re-impose a group.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct MergeDuplicates {
    tolerance: f64,
}

impl MergeDuplicates {
    /// Creates the command with the reference 0.05 Angstrom tolerance.
    #[must_use]
    pub fn new() -> Self {
        Self { tolerance: 0.05 }
    }

    /// Sets the merge tolerance in Angstrom.
    #[must_use]
    pub fn with_tolerance(mut self, tolerance: f64) -> Self {
        self.tolerance = tolerance;
        self
    }

    /// Merge tolerance in Angstrom.
    #[must_use]
    pub fn tolerance(&self) -> f64 {
        self.tolerance
    }
}

impl Default for MergeDuplicates {
    fn default() -> Self {
        Self::new()
    }
}

impl Command for MergeDuplicates {
    fn apply(&self, structure: &Structure) -> Result<Structure, Error> {
        xtalforge_core::merge_duplicates(structure, self.tolerance)
    }

    fn description(&self) -> String {
        "Merge duplicate sites".to_owned()
    }
}
