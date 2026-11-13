//! Undoable edits that come out of the force field (Task 14.3.3).
//!
//! An optimisation moves every site at once and lands as a single undo step: a
//! relaxation is one thing the user asked for, and Ctrl+Z must give back the
//! structure they started from. An atom-type override is an edit like any
//! other because it changes the answer the next calculation gives. Adding
//! hydrogens is here for the same reason: what it does to the structure is an
//! ordinary edit, and the force field decides *what* to add.

use serde::{Deserialize, Serialize};
use serde_json::json;
use xtalforge_core::{Error, Lattice, Site, Structure};

use crate::command::Command;
use crate::edit::{check_index, rebuild, rebuild_explicit, rebuild_with};

/// Site prop used by the force field's atom-type override.
pub const ATOM_TYPE_KEY: &str = "uff_type";

/// Moves every site to where the optimizer left it (Task 14.3.3).
///
/// A variable-cell relaxation stores the new lattice in the same command, so
/// one undo puts both the atoms and the cell back together. Deliberately not a
/// [`crate::atoms::MoveSites`]: two optimisations in a row are two things that
/// happened, and collapsing them would make the first unreachable.
#[derive(Debug, Clone, PartialEq)]
pub struct ApplyOptimizedGeometry {
    frac: Vec<[f64; 3]>,
    lattice: Option<Lattice>,
    description: String,
}

impl ApplyOptimizedGeometry {
    /// Creates the command from one fractional position per site.
    #[must_use]
    pub fn new(frac: Vec<[f64; 3]>) -> Self {
        Self {
            frac,
            lattice: None,
            description: "Apply optimized geometry".to_owned(),
        }
    }

    /// Carries the relaxed lattice with the coordinates.
    #[must_use]
    pub fn with_lattice(mut self, lattice: Lattice) -> Self {
        self.lattice = Some(lattice);
        self
    }

    /// Overrides the description.
    #[must_use]
    pub fn with_description(mut self, description: impl Into<String>) -> Self {
        self.description = description.into();
        self
    }

    /// New fractional positions.
    #[must_use]
    pub fn fracs(&self) -> &[[f64; 3]] {
        &self.frac
    }

    /// Relaxed lattice, when the cell moved.
    #[must_use]
    pub fn lattice(&self) -> Option<&Lattice> {
        self.lattice.as_ref()
    }

    /// Largest Cartesian distance any atom would move, in Angstrom.
    ///
    /// The new positions are measured in the cell they land in.
    ///
    /// # Errors
    ///
    /// Returns [`Error::InvalidInput`] when the optimizer returned a different
    /// number of sites than the structure holds.
    pub fn max_displacement(&self, structure: &Structure) -> Result<f64, Error> {
        if self.frac.len() != structure.sites().len() {
            return Err(Error::InvalidInput(format!(
                "the optimizer returned {} sites but the structure has {}",
                self.frac.len(),
                structure.sites().len()
            )));
        }
        let lattice = self.lattice.unwrap_or(*structure.lattice());
        let mut largest: f64 = 0.0;
        for (site, frac) in structure.sites().iter().zip(self.frac.iter()) {
            let before = structure.lattice().to_cart(site.frac());
            let after = lattice.to_cart(*frac);
            let distance = ((after[0] - before[0]).powi(2)
                + (after[1] - before[1]).powi(2)
                + (after[2] - before[2]).powi(2))
            .sqrt();
            largest = largest.max(distance);
        }
        Ok(largest)
    }
}

impl Command for ApplyOptimizedGeometry {
    fn apply(&self, structure: &Structure) -> Result<Structure, Error> {
        if self.frac.len() != structure.sites().len() {
            return Err(Error::InvalidInput(format!(
                "the optimizer returned {} sites but the structure has {}; it was edited while the calculation ran",
                self.frac.len(),
                structure.sites().len()
            )));
        }
        let mut sites = Vec::with_capacity(structure.sites().len());
        for (site, frac) in structure.sites().iter().zip(self.frac.iter()) {
            sites.push(site.clone().with_frac(*frac)?);
        }
        let lattice = self.lattice.unwrap_or(*structure.lattice());
        rebuild_with(
            structure,
            sites,
            lattice,
            structure.hall(),
            structure.operations(),
        )
    }

    fn description(&self) -> String {
        self.description.clone()
    }
}

/// Appends hydrogen sites and their bonds to the parent atoms (Task 14.3.3).
///
/// The placement itself belongs to the force field; this command records the
/// result as one undoable edit. Only hydrogen and deuterium sites are accepted.
/// A structure with symmetry must be reduced to P1 before bonds can be named,
/// because an explicit bond indexes P1 atoms.
#[derive(Debug, Clone, PartialEq)]
pub struct AddHydrogens {
    sites: Vec<Site>,
    bonds: Vec<(usize, usize)>,
}

impl AddHydrogens {
    /// Creates the command from the new sites.
    #[must_use]
    pub fn new(sites: Vec<Site>) -> Self {
        Self {
            sites,
            bonds: Vec::new(),
        }
    }

    /// Adds explicit bonds as `(existing site index, new site position)`.
    #[must_use]
    pub fn with_bonds(mut self, bonds: Vec<(usize, usize)>) -> Self {
        self.bonds = bonds;
        self
    }

    /// Adds one explicit bond.
    #[must_use]
    pub fn with_bond(mut self, parent: usize, hydrogen: usize) -> Self {
        self.bonds.push((parent, hydrogen));
        self
    }

    /// New sites the command appends.
    #[must_use]
    pub fn sites(&self) -> &[Site] {
        &self.sites
    }

    /// Bonds to the parent atoms as `(existing site index, new site position)`.
    #[must_use]
    pub fn bonds(&self) -> &[(usize, usize)] {
        &self.bonds
    }
}

impl Command for AddHydrogens {
    fn apply(&self, structure: &Structure) -> Result<Structure, Error> {
        if self.sites.is_empty() {
            return Err(Error::InvalidInput(
                "at least one hydrogen is required".to_owned(),
            ));
        }
        for site in &self.sites {
            if site.element() != "H" && site.element() != "D" {
                return Err(Error::InvalidInput(format!(
                    "add hydrogens only accepts H or D sites, got {}",
                    site.element()
                )));
            }
        }
        for &(parent, hydrogen) in &self.bonds {
            check_index(parent, structure.sites().len(), "site")?;
            check_index(hydrogen, self.sites.len(), "hydrogen")?;
        }
        if !self.bonds.is_empty() && !structure.is_p1() {
            return Err(Error::InvalidInput(
                "adding bonded hydrogens requires a P1 structure; reduce to P1 first".to_owned(),
            ));
        }
        let base = structure.sites().len();
        let mut sites = structure.sites().to_vec();
        sites.extend(self.sites.iter().cloned());
        let with_sites = rebuild(structure, sites, *structure.lattice())?;
        if self.bonds.is_empty() {
            return Ok(with_sites);
        }
        let mut bonds = with_sites.explicit_bonds().to_vec();
        for &(parent, hydrogen) in &self.bonds {
            bonds.push(xtalforge_core::ExplicitBond {
                i: parent,
                j: base + hydrogen,
                image: [0, 0, 0],
                order: 1.0,
                remove: false,
            });
        }
        rebuild_explicit(&with_sites, bonds)
    }

    fn description(&self) -> String {
        format!("Add {} hydrogens", self.sites.len())
    }
}

/// One recorded preparation step.
#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
pub struct PrepareStep {
    /// Stable step name.
    pub name: String,
    /// Whether the step ran.
    pub done: bool,
    /// Human-readable detail.
    pub detail: String,
}

impl PrepareStep {
    /// Creates one recorded step.
    #[must_use]
    pub fn new(name: impl Into<String>, done: bool, detail: impl Into<String>) -> Self {
        Self {
            name: name.into(),
            done,
            detail: detail.into(),
        }
    }
}

/// Installs the result of the preparation pipeline as one undo step
/// (Task 14.3.3).
///
/// The pipeline needs a calculator and lives in `xtalforge_ff::prepare`; a
/// [`Command`] stores no trait objects, so the caller runs it and hands the
/// prepared structure and its step log to this command. One Ctrl+Z gives back
/// the structure as deposited: undoing half the steps would leave a cell nobody
/// chose.
#[derive(Debug, Clone, PartialEq)]
pub struct Prepare {
    prepared: Structure,
    steps: Vec<PrepareStep>,
}

impl Prepare {
    /// Creates the command from a prepared structure and its recorded steps.
    ///
    /// # Errors
    ///
    /// Returns [`Error::InvalidInput`] when no step is recorded or every step
    /// is skipped.
    pub fn new(prepared: Structure, steps: Vec<PrepareStep>) -> Result<Self, Error> {
        if steps.is_empty() {
            return Err(Error::InvalidInput(
                "no preparation step is recorded".to_owned(),
            ));
        }
        if !steps.iter().any(|step| step.done) {
            return Err(Error::InvalidInput(
                "every preparation step is skipped".to_owned(),
            ));
        }
        Ok(Self { prepared, steps })
    }

    /// Structure installed by the command.
    #[must_use]
    pub fn prepared(&self) -> &Structure {
        &self.prepared
    }

    /// Recorded steps in pipeline order.
    #[must_use]
    pub fn steps(&self) -> &[PrepareStep] {
        &self.steps
    }
}

impl Command for Prepare {
    fn apply(&self, _structure: &Structure) -> Result<Structure, Error> {
        Ok(self.prepared.clone())
    }

    fn description(&self) -> String {
        "Prepare for simulation".to_owned()
    }
}

/// Overrides the force field's choice of atom type (Task 14.3.3).
///
/// The type is stored on the site as [`ATOM_TYPE_KEY`], so it applies to the
/// whole symmetry orbit and survives a save. Passing `None` clears the override
/// and hands the site back to the typer.
#[derive(Debug, Clone, PartialEq)]
pub struct SetAtomTypes {
    indices: Vec<usize>,
    type_name: Option<String>,
}

impl SetAtomTypes {
    /// Creates the command from site indexes and an optional type name.
    #[must_use]
    pub fn new(indices: Vec<usize>, type_name: Option<&str>) -> Self {
        Self {
            indices,
            type_name: type_name.map(str::to_owned),
        }
    }

    /// Indexes the command edits.
    #[must_use]
    pub fn indices(&self) -> &[usize] {
        &self.indices
    }

    /// Type name, or `None` when the override is cleared.
    #[must_use]
    pub fn type_name(&self) -> Option<&str> {
        self.type_name.as_deref()
    }
}

impl Command for SetAtomTypes {
    fn apply(&self, structure: &Structure) -> Result<Structure, Error> {
        let mut sites = structure.sites().to_vec();
        for index in &self.indices {
            check_index(*index, sites.len(), "site")?;
            let Some(site) = sites.get(*index) else {
                continue;
            };
            let updated = match &self.type_name {
                Some(name) => site.clone().with_prop(ATOM_TYPE_KEY, json!(name)),
                None => site.clone().without_prop(ATOM_TYPE_KEY),
            };
            if let Some(slot) = sites.get_mut(*index) {
                *slot = updated;
            }
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
        match &self.type_name {
            Some(name) => format!("Set atom type to {name}"),
            None => "Clear atom type override".to_owned(),
        }
    }
}
