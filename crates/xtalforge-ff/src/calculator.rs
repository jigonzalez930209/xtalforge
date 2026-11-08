//! The energy-engine interface, options, registry and run records (Task 4.1).
//!
//! The shape follows Crystal-Builder `xtal/ff/api.py` and `registry.py`: an
//! engine is asked for one energy, or an energy and its forces, for a fixed
//! structure. Topology and atom types belong to the engine implementation and
//! are built once per [`Calculator::energy_forces`] call by the current UFF
//! model; [`crate::UffModel`] exposes the prepared form directly for callers
//! that move the geometry many times.
//!
//! Units at this boundary are eV, eV/Angstrom and GPa. The UFF implementation
//! converts its internal kcal/mol at the last step so no caller has to know
//! which force field produced a number.

use serde::{Deserialize, Serialize};
use xtalforge_core::{BondRules, Structure};

use crate::Error;
use crate::uff::params::ParameterSet;

/// Everything about a calculation that does not come from the structure.
#[derive(Debug, Clone, PartialEq)]
pub struct CalculatorOptions {
    /// UFF parameter table to use.
    pub parameter_set: ParameterSet,
    /// Whether the Lennard-Jones term is evaluated.
    pub vdw: bool,
    /// van der Waals cutoff in Angstrom.
    pub vdw_cutoff: f64,
    /// Rules for distance-based bond perception.
    pub bond_rules: BondRules,
}

impl Default for CalculatorOptions {
    fn default() -> Self {
        Self {
            parameter_set: ParameterSet::default(),
            vdw: true,
            vdw_cutoff: crate::uff::DEFAULT_VDW_CUTOFF,
            bond_rules: BondRules::default(),
        }
    }
}

impl CalculatorOptions {
    /// Checks that the options are usable.
    ///
    /// # Errors
    ///
    /// Returns [`Error::InvalidInput`] for a non-finite or non-positive
    /// cutoff, or for bond rules the core bond perception rejects.
    pub fn validate(&self) -> Result<(), Error> {
        if !self.vdw_cutoff.is_finite() || self.vdw_cutoff <= 0.0 {
            return Err(Error::InvalidInput(format!(
                "van der Waals cutoff must be positive and finite, got {}",
                self.vdw_cutoff
            )));
        }
        if !self.bond_rules.scale.is_finite() || self.bond_rules.scale <= 0.0 {
            return Err(Error::InvalidInput(
                "bond scale must be positive and finite".to_owned(),
            ));
        }
        if !self.bond_rules.minimum_distance.is_finite() || self.bond_rules.minimum_distance < 0.0 {
            return Err(Error::InvalidInput(
                "bond minimum distance must be non-negative and finite".to_owned(),
            ));
        }
        Ok(())
    }
}

/// An engine that answers energy (and optionally forces and stress) queries.
///
/// Implementations are stateless with respect to the structure: every call
/// receives the structure it should evaluate. The prepared UFF form,
/// [`crate::UffModel`], is the exception offered to optimisation loops that
/// keep one topology while positions move.
pub trait Calculator: Send + Sync {
    /// Stable registry name, for example `"uff"`.
    fn name(&self) -> &'static str;

    /// Energy in eV.
    ///
    /// # Errors
    ///
    /// Returns a typing, topology or data error when the structure cannot be
    /// given an energy.
    fn energy(&self, structure: &Structure) -> Result<f64, Error>;

    /// Energy in eV and forces in eV/Angstrom, one vector per P1 atom.
    ///
    /// Dummy `X` atoms are excluded from the energy but keep a zero force in
    /// the returned vector, so the result always indexes like
    /// [`xtalforge_core::P1Cell::atoms`].
    ///
    /// # Errors
    ///
    /// Returns a typing, topology or data error when the structure cannot be
    /// evaluated.
    fn energy_forces(&self, structure: &Structure) -> Result<(f64, Vec<[f64; 3]>), Error>;

    /// Whether [`Calculator::energy_forces`] returns forces.
    fn provides_forces(&self) -> bool {
        true
    }

    /// Whether the engine can return an analytic stress.
    fn provides_stress(&self) -> bool {
        false
    }

    /// Stress tensor in GPa (row `a`, column `b` is the force along `b` on
    /// the face normal to `a`), or `None` when unsupported.
    ///
    /// # Errors
    ///
    /// Returns the errors of [`Calculator::energy`] when the engine evaluates
    /// the stress from the same topology.
    fn stress(&self, _structure: &Structure) -> Result<Option<[[f64; 3]; 3]>, Error> {
        Ok(None)
    }
}

/// Whether an engine can run in this build and environment.
#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
pub enum Availability {
    /// The engine is ready.
    Available,
    /// The engine is registered but cannot run; the string says why.
    Unavailable(String),
}

impl Availability {
    /// Whether the engine is available.
    #[must_use]
    pub fn is_available(&self) -> bool {
        matches!(self, Self::Available)
    }

    /// Reason for unavailability, when there is one.
    #[must_use]
    pub fn reason(&self) -> Option<&str> {
        match self {
            Self::Available => None,
            Self::Unavailable(reason) => Some(reason),
        }
    }
}

/// What an engine can do, advertised before it is constructed.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
pub struct CalculatorCapabilities {
    /// Understands periodic boundary conditions.
    pub periodic: bool,
    /// Returns forces.
    pub forces: bool,
    /// Returns an analytic stress.
    pub stress: bool,
    /// Uses per-atom charges.
    pub charges: bool,
    /// Runs in-process rather than launching a program.
    pub in_process: bool,
    /// Can be called from several threads with one instance.
    pub thread_safe: bool,
}

/// Registry description of one engine.
#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
pub struct CalculatorInfo {
    /// Stable registry name.
    pub name: String,
    /// Human-readable label.
    pub label: String,
    /// One-sentence description.
    pub description: String,
    /// Availability in this build.
    pub availability: Availability,
    /// Advertised capabilities.
    pub capabilities: CalculatorCapabilities,
}

impl CalculatorInfo {
    /// Whether the engine is available.
    #[must_use]
    pub fn available(&self) -> bool {
        self.availability.is_available()
    }
}

/// Builds a calculator from validated options.
pub type CalculatorFactory = fn(CalculatorOptions) -> Result<Box<dyn Calculator>, Error>;

#[derive(Clone, Copy)]
struct RegistryEntry {
    info: &'static CalculatorInfo,
    build: CalculatorFactory,
}

/// Name-to-factory table read by the CLI, bindings and agent tools.
#[derive(Clone)]
pub struct CalculatorRegistry {
    entries: Vec<RegistryEntry>,
}

impl Default for CalculatorRegistry {
    fn default() -> Self {
        crate::registry::builtin_registry()
    }
}

impl CalculatorRegistry {
    /// An empty registry.
    #[must_use]
    pub fn new() -> Self {
        Self {
            entries: Vec::new(),
        }
    }

    /// The registry with every built-in engine registered.
    #[must_use]
    pub fn builtin() -> Self {
        crate::registry::builtin_registry()
    }

    /// Adds one engine.
    ///
    /// # Errors
    ///
    /// Returns [`Error::InvalidInput`] when the name is empty or already
    /// present. Re-registering an engine is a programming mistake, not a
    /// configuration choice.
    pub fn register(
        &mut self,
        info: &'static CalculatorInfo,
        build: CalculatorFactory,
    ) -> Result<(), Error> {
        if info.name.trim().is_empty() {
            return Err(Error::InvalidInput(
                "calculator name must not be empty".to_owned(),
            ));
        }
        if self
            .entries
            .iter()
            .any(|entry| entry.info.name == info.name)
        {
            return Err(Error::InvalidInput(format!(
                "calculator {} is already registered",
                info.name
            )));
        }
        self.entries.push(RegistryEntry { info, build });
        Ok(())
    }

    /// Description of one engine.
    #[must_use]
    pub fn info(&self, name: &str) -> Option<&'static CalculatorInfo> {
        self.entries
            .iter()
            .find(|entry| entry.info.name == name)
            .map(|entry| entry.info)
    }

    /// Adds a built-in entry whose name and description are already known to
    /// be unique.
    pub(crate) fn add_builtin(&mut self, info: &'static CalculatorInfo, build: CalculatorFactory) {
        self.entries.push(RegistryEntry { info, build });
    }

    /// Every registered engine, in registration order.
    pub fn infos(&self) -> impl Iterator<Item = &'static CalculatorInfo> + '_ {
        self.entries.iter().map(|entry| entry.info)
    }

    /// Registered names, in registration order.
    #[must_use]
    pub fn names(&self) -> Vec<&'static str> {
        self.entries
            .iter()
            .map(|entry| entry.info.name.as_str())
            .collect()
    }

    /// Builds one engine.
    ///
    /// # Errors
    ///
    /// Returns [`Error::UnknownCalculator`] for an unregistered name,
    /// [`Error::InvalidInput`] for invalid options, or an error from the
    /// factory when the engine cannot start.
    pub fn build(
        &self,
        name: &str,
        options: &CalculatorOptions,
    ) -> Result<Box<dyn Calculator>, Error> {
        options.validate()?;
        let entry = self
            .entries
            .iter()
            .find(|entry| entry.info.name == name)
            .ok_or_else(|| Error::UnknownCalculator(name.to_owned()))?;
        if let Availability::Unavailable(reason) = &entry.info.availability {
            return Err(Error::InvalidInput(format!(
                "calculator {name} is unavailable: {reason}"
            )));
        }
        (entry.build)(options.clone())
    }
}

impl std::fmt::Debug for CalculatorRegistry {
    fn fmt(&self, formatter: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        formatter
            .debug_struct("CalculatorRegistry")
            .field("names", &self.names())
            .finish()
    }
}

/// What a run writes down for the workspace and report system (Task 4.1.4).
///
/// Ported from Crystal-Builder `xtal/ff/record.py`, reduced to the fields the
/// Phase 4 surfaces consume first: the engine, the final energy, the cell
/// size, convergence, the number of steps and human-readable notes. Energies
/// are eV, as returned by [`Calculator`].
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct RunRecord {
    /// Registry name of the engine that produced the result.
    pub calculator: String,
    /// Final energy in eV.
    pub energy: f64,
    /// Number of atoms the engine evaluated (dummy atoms excluded).
    pub n_atoms: usize,
    /// Whether the run reported convergence.
    pub converged: bool,
    /// Number of optimisation or evaluation steps.
    pub steps: usize,
    /// Warnings and decisions worth keeping with the run.
    pub notes: Vec<String>,
}

impl RunRecord {
    /// A record for one completed energy evaluation.
    #[must_use]
    pub fn new(calculator: impl Into<String>, energy: f64, n_atoms: usize) -> Self {
        Self {
            calculator: calculator.into(),
            energy,
            n_atoms,
            converged: true,
            steps: 1,
            notes: Vec::new(),
        }
    }

    /// Adds a note.
    #[must_use]
    pub fn with_note(mut self, note: impl Into<String>) -> Self {
        self.notes.push(note.into());
        self
    }

    /// Marks the run as unconverged after `steps` steps.
    #[must_use]
    pub fn with_steps(mut self, steps: usize, converged: bool) -> Self {
        self.steps = steps;
        self.converged = converged;
        self
    }

    /// Serializes the record to JSON.
    ///
    /// # Errors
    ///
    /// Returns [`Error::Json`] when serialization fails.
    pub fn to_json(&self) -> Result<String, Error> {
        serde_json::to_string(self).map_err(|error| Error::Json(error.to_string()))
    }

    /// Parses a record from JSON.
    ///
    /// # Errors
    ///
    /// Returns [`Error::Json`] for invalid JSON or unknown fields.
    pub fn from_json(text: &str) -> Result<Self, Error> {
        serde_json::from_str(text).map_err(|error| Error::Json(error.to_string()))
    }
}
