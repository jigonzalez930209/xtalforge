//! Force fields, charges, optimizers, dynamics and thermochemistry for
//! xtalforge.
//!
//! Phase 4 starts with the energy-engine interface ([`Calculator`]) and the
//! native Universal Force Field. The UFF implementation is ported from
//! Crystal-Builder `xtal/ff/uff/` (MIT) and keeps the reference units
//! (Angstrom, degrees, kcal/mol) internally; the public API reports energies
//! in eV, forces in eV/Angstrom and stress in GPa.
//!
//! [`vibrations`] builds on that boundary: a Gamma-point Hessian from central
//! differences of analytic forces, mass-weighted modes through a Jacobi
//! eigensolver, IR intensities from a documented fixed-charge dipole model or
//! from Born effective-charge tensors, powder-averaged Raman activities from
//! engine tensors, Gamma-point phonons with a DOS and inelastic-neutron
//! weights, harmonic thermochemistry and a sampled-volume quasi-harmonic free
//! energy.  Every approximation is named in the module that makes it.
//!
//! Task 14.2 adds the scientific-completeness layer on top of Phase 4:
//! [`CellFreedom`] holds selected cell lengths, angles or the volume and
//! `OptimizeOptions::pressure_gpa` adds the `P V` term; [`RunControl`] cancels
//! a run and reports progress, and [`Trajectory`] records frames from every
//! optimizer and scan; [`ScanCoordinate`] covers cell parameters and collective
//! variables in addition to internal coordinates; [`ScanCheckpoint`] makes a
//! relaxed scan resumable; and [`md`] provides NVE, NVT Langevin, NVT
//! Nose-Hoover and NPT molecular dynamics with explicit deterministic seeds.
//! [`try_run_engine`] is the opportunistic DFTB+/xTB path: it probes first and
//! runs only when the binary is installed, so no build or test needs one.

pub mod calculator;
pub mod charges;
pub mod elastic;
pub mod external;
pub mod md;
pub mod optimize;
pub mod prepare;
pub mod quantum;
pub mod registry;
pub mod run;
pub mod scan;
pub mod trajectory;
pub mod uff;
pub mod vibrations;

pub mod models;

mod math;
mod structure_util;

pub use calculator::{
    Availability, Calculator, CalculatorCapabilities, CalculatorInfo, CalculatorOptions,
    CalculatorRegistry, RunRecord,
};
pub use charges::{
    ChargeOptions, Charges, EwaldOptions, EwaldResult, apply_charges, eeq_charges, eqeq_charges,
    ewald_energy, ewald_energy_forces, qeq_charges,
};
pub use elastic::{ElasticOptions, ElasticTensor, deformation_set, elastic_tensor};
pub use external::{
    BOHR_TO_ANGSTROM, DftbCalculator, DftbOptions, EngineKind, EngineResult, EngineRunOptions,
    EngineRunStatus, EngineStatus, HARTREE_PER_BOHR_TO_EV_PER_ANGSTROM, HARTREE_PER_BOHR3_TO_GPA,
    HARTREE_TO_EV, OrcaOptions, PROBE_ENGINES, RunOutcome, RunSpec, RunStatus, XtbCalculator,
    XtbOptions, available_engines, orca_catalogue, parse_dftb_output, parse_xtb_output, probe,
    probe_binary, run_in_scratch, try_run_engine, write_dftb_hsd, write_orca_input,
    write_xtb_input,
};
pub use md::{
    AMU_TO_EV_FS2_PER_A2, Ensemble, KB_EV_PER_K, MdOptions, MdResult, MdSample, molecular_dynamics,
    molecular_dynamics_with,
};
pub use models::{CATALOGUE as MODEL_CATALOGUE, ModelRecord, ModelSpec, ModelStatus};
pub use optimize::{
    CellFreedom, EosFit, EosForm, Minimizer, OptimizationResult, OptimizeOptions, RelaxMode,
    apply_optimized_geometry, birch_murnaghan, fire, fire_with, fit_eos, lbfgs, lbfgs_with,
    minimize, minimize_with, murnaghan,
};
pub use prepare::{PrepareOptions, PrepareReport, PrepareStep, prepare};
pub use quantum::{QuantumAdapter, QuantumCapabilities, QuantumMethod, adapter, capabilities};
pub use registry::builtin_registry;
pub use run::{CancelToken, ProgressEvent, RunControl};
pub use scan::{
    CellParameter, CollectiveVariable, Constraint, ScanAxis, ScanCheckpoint, ScanCoordinate,
    ScanGeometry, ScanOptions, ScanPoint, ScanResult, relaxed_scan, relaxed_scan_with,
};
pub use trajectory::{Trajectory, TrajectoryFrame, TrajectoryOptions};
pub use uff::calculator::UffCalculator;
pub use uff::model::UffModel;
pub use uff::params::{ParameterSet, UffParams};
pub use uff::typer::{AtomType, Confidence, Typing};
pub use vibrations::{
    BornCharges, Broadening, CM_PER_SQRT_EV_PER_A2_AMU, CM_TO_EV, DEFAULT_HESSIAN_STEP,
    DEFAULT_IR_FWHM, DEFAULT_IR_STEP, DEFAULT_RAMAN_FWHM, DEFAULT_RAMAN_STEP, DisplacementMode,
    DosOptions, Hessian, HessianOptions, IMAGINARY_BELOW_CM, InelasticNeutronDos, IrOptions,
    IrPeak, IrSpectrum, MAX_IR_POINTS, MAX_RAMAN_POINTS, NeutronWeights, PhononDos, PhononOptions,
    PhononResult, QhaFreeEnergy, QhaPoint, QhaResult, RamanOptions, RamanPeak, RamanSpectrum,
    RamanTensorSet, TensorFile, ThermoOptions, ThermoResult, VibrationalModes,
    ZERO_MODE_TOLERANCE_CM, hessian, ideal_gas_thermochemistry, inelastic_neutron_dos, ir_spectrum,
    ir_spectrum_with_born, mode_neutron_weights, parse_tensor_file, phonon_dos, phonons,
    powder_average, quasi_harmonic_free_energy, raman_spectrum, thermochemistry, vibrational_modes,
};

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
    /// An in-process quantum backend failed.
    #[error("quantum: {0}")]
    Quantum(String),
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
            Self::Quantum(_) => "quantum",
            Self::Json(_) => "json",
            Self::Core(error) => error.code(),
        }
    }

    /// Human-readable message.
    #[must_use]
    pub fn message(&self) -> String {
        self.to_string()
    }

    /// Stable machine-readable category for the JSON error object.
    #[must_use]
    pub fn category(&self) -> &'static str {
        match self {
            Self::InvalidInput(_)
            | Self::Typing(_)
            | Self::UnknownElement(_)
            | Self::UnknownType(_)
            | Self::Topology(_)
            | Self::UnknownCalculator(_) => "input",
            Self::Data(_) | Self::Json(_) => "internal",
            Self::Quantum(_) => "engine",
            Self::Core(error) => error.category(),
        }
    }
}
