//! Force fields, charges, optimizers and thermochemistry for xtalforge.
//!
//! Phase 4 starts with the energy-engine interface ([`Calculator`]) and the
//! native Universal Force Field. The UFF implementation is ported from
//! Crystal-Builder `xtal/ff/uff/` (MIT) and keeps the reference units
//! (Angstrom, degrees, kcal/mol) internally; the public API reports energies
//! in eV, forces in eV/Angstrom and stress in GPa.

pub mod calculator;
pub mod charges;
pub mod optimize;
pub mod registry;
pub mod scan;
pub mod uff;

mod math;
mod structure_util;

pub use calculator::{
    Availability, Calculator, CalculatorCapabilities, CalculatorInfo, CalculatorOptions,
    CalculatorRegistry, RunRecord,
};
pub use charges::{
    ChargeOptions, Charges, EwaldOptions, EwaldResult, apply_charges, eqeq_charges, ewald_energy,
    ewald_energy_forces, qeq_charges,
};
pub use optimize::{
    EosFit, EosForm, Minimizer, OptimizationResult, OptimizeOptions, RelaxMode,
    apply_optimized_geometry, birch_murnaghan, fire, fit_eos, lbfgs, minimize, murnaghan,
};
pub use registry::builtin_registry;
pub use scan::{
    Constraint, ScanAxis, ScanCoordinate, ScanOptions, ScanPoint, ScanResult, relaxed_scan,
};
pub use uff::calculator::UffCalculator;
pub use uff::model::UffModel;
pub use uff::params::{ParameterSet, UffParams};
pub use uff::typer::{AtomType, Confidence, Typing};

/// Errors returned by force-field operations.
#[derive(Debug, Clone, PartialEq, thiserror::Error)]
pub enum Error {
    /// An argument does not satisfy the documented preconditions.
    #[error("invalid input: {0}")]
    InvalidInput(String),
    /// An atom has no UFF parameters or its environment cannot be typed.
    #[error("atom typing: {0}")]
    Typing(String),
    /// An element has no entry in the UFF parameter table.
    #[error("no UFF parameters for element {0}")]
    UnknownElement(String),
    /// A type name is not present in the UFF parameter table.
    #[error("unknown UFF atom type {0}")]
    UnknownType(String),
    /// The topology could not be built or is inconsistent.
    #[error("topology: {0}")]
    Topology(String),
    /// The embedded parameter data is missing or corrupt.
    #[error("parameter data: {0}")]
    Data(String),
    /// A calculator name is not registered.
    #[error("unknown calculator {0}")]
    UnknownCalculator(String),
    /// Serialization or deserialization of a record failed.
    #[error("JSON: {0}")]
    Json(String),
    /// An error from the core crate.
    #[error("core: {0}")]
    Core(#[from] xtalforge_core::Error),
}

impl Error {
    /// Stable `snake_case` code for this error, suitable for bindings and JSON.
    ///
    /// The string is part of the public API (see `AGENTS.md`, section 5.2).
    #[must_use]
    pub fn code(&self) -> &'static str {
        match self {
            Self::InvalidInput(_) => "invalid_input",
            Self::Typing(_) => "typing",
            Self::UnknownElement(_) => "unknown_element",
            Self::UnknownType(_) => "unknown_type",
            Self::Topology(_) => "topology",
            Self::Data(_) => "data",
            Self::UnknownCalculator(_) => "unknown_calculator",
            Self::Json(_) => "json",
            Self::Core(error) => error.code(),
        }
    }

    /// Human-readable message.
    #[must_use]
    pub fn message(&self) -> String {
        self.to_string()
    }
}
