//! External quantum engines: probes, input writers and output parsers
//! (Task 4.9).
//!
//! The engines are optional programs; nothing in the default build launches
//! one.  This module owns the parts that can be tested without a binary:
//!
//! * [`probe`] searches `PATH`, checks that the file is executable and runs
//!   `<binary> --version` (short output read).  It never panics; every failure
//!   is reported as [`EngineStatus::Missing`].  WASM has no process table, so
//!   it always reports "not available in the browser" (the recorded exception
//!   in the parity matrix).
//! * [`write_dftb_hsd`] and [`parse_dftb_output`] cover DFTB+ `GenFormat` HSD
//!   input and `detailed.out` energy, forces and stress.
//! * [`write_xtb_input`] and [`parse_xtb_output`] cover the xTB coordinate
//!   file and its `TOTAL ENERGY`, `TOTAL CHARGE` and force blocks.
//! * [`write_orca_input`] covers ORCA's simple-input format and makes every
//!   molecule whole across cell faces by translating bonded images, then
//!   shifts each connected component's centroid back into the cell.
//! * [`orca_catalogue`] lists the functionals and bases offered by the UI.
//! * [`try_run_engine`] is the opportunistic path: it probes first, skips a
//!   missing binary and runs DFTB+ or xTB when one is installed (Task 14.2.6).
//!   [`available_engines`] reports the probe result for every known engine.
//!
//! **Never required.**  No default build and no test needs an installed
//! engine; the opportunistic path reports [`EngineRunStatus::Skipped`] instead
//! of failing.
//!
//! **Units and documented assumptions.**  Energies leave this module in eV and
//! forces in eV/Angstrom, like every other `Calculator` boundary in the crate;
//! the parsers convert Hartree and Hartree/Bohr.  A stress tensor printed in
//! Hartree/Bohr^3 is converted to GPa, one printed in GPa is kept.  Writers
//! expand the input to P1, because every external engine expects the full
//! cell, and refuse the dummy `X` species.  DFTB+ forces are optional: a
//! `detailed.out` without a `Total Forces` block yields an empty vector, while
//! a missing `Total energy:` is an error.  The xTB parser validates
//! `TOTAL CHARGE` when present but does not return it because [`EngineResult`]
//! is the shared, engine-neutral record.  xTB selects its method on the command
//! line, not in the input file, so [`write_xtb_input`] records the method in the
//! XYZ comment for reproducibility.  DFTB+ uses the standard
//! `slakos/<Element1>-<Element2>.skf` naming under the `slakos/` prefix; the
//! user's chosen parameter set must live there (or the `Prefix` in the written
//! file must be edited).

mod available;
mod calculators;
mod dftb;
mod orca;
mod probe;
mod runner;
mod text;
mod xtb;

use serde::{Deserialize, Serialize};

pub use available::{
    EngineRunOptions, EngineRunStatus, PROBE_ENGINES, available_engines, try_run_engine,
};
pub use calculators::{DftbCalculator, XtbCalculator};
pub use dftb::{DftbOptions, parse_dftb_output, write_dftb_hsd};
pub use orca::{OrcaOptions, orca_catalogue, write_orca_input};
pub use probe::{probe, probe_binary};
pub use runner::{RunOutcome, RunSpec, RunStatus, run_in_scratch};
pub use xtb::{XtbOptions, parse_xtb_output, write_xtb_input};

/// Hartree in eV.
pub const HARTREE_TO_EV: f64 = 27.211_386_245_988;
/// Bohr in Angstrom.
pub const BOHR_TO_ANGSTROM: f64 = 0.529_177_210_903;
/// Hartree/Bohr in eV/Angstrom.
pub const HARTREE_PER_BOHR_TO_EV_PER_ANGSTROM: f64 = HARTREE_TO_EV / BOHR_TO_ANGSTROM;
/// Hartree/Bohr^3 in GPa.
pub const HARTREE_PER_BOHR3_TO_GPA: f64 = 29_421.010_7;

/// One supported external engine.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash, Serialize, Deserialize)]
pub enum EngineKind {
    /// DFTB+.
    #[serde(rename = "dftb+")]
    DftbPlus,
    /// xTB.
    #[serde(rename = "xtb")]
    Xtb,
    /// ORCA.
    #[serde(rename = "orca")]
    Orca,
}

impl EngineKind {
    /// Executable name searched on `PATH`.
    #[must_use]
    pub const fn binary(self) -> &'static str {
        match self {
            Self::DftbPlus => "dftb+",
            Self::Xtb => "xtb",
            Self::Orca => "orca",
        }
    }
}

/// Result of probing one engine.
#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
#[serde(rename_all = "snake_case")]
pub enum EngineStatus {
    /// The binary was found and answered `--version`.
    Available {
        /// Absolute or PATH-relative path that was executed.
        path: String,
        /// First non-empty line of the version output, capped at 160 chars.
        version: String,
    },
    /// The binary is not usable in this build or environment.
    Missing {
        /// Why the engine is missing; safe to show to a user.
        reason: String,
    },
}

/// Energy, forces and optional stress from one external run.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct EngineResult {
    /// Final energy in eV.
    pub energy: f64,
    /// Forces in eV/Angstrom, one vector per atom; empty when not printed.
    pub forces: Vec<[f64; 3]>,
    /// Stress tensor in GPa, when the engine printed one.
    pub stress: Option<[[f64; 3]; 3]>,
}
