//! Opportunistic external-engine runs (Task 14.2.6).
//!
//! The Task 4.9 calculators always try to launch a binary.  This module adds
//! the "run it when it exists" path a local CI or a batch workflow needs: probe
//! first, skip when the binary is missing, and report why.  Nothing here is
//! required by the default build or by the test suite; a machine without DFTB+
//! or xTB simply gets [`EngineRunStatus::Skipped`].
//!
//! ```text
//! match try_run_engine(EngineKind::Xtb, &structure, &options, &cancel)? {
//!     EngineRunStatus::Skipped { reason } => eprintln!("xTB not run: {reason}"),
//!     EngineRunStatus::Completed(result) => { /* energy, forces, stress */ }
//!     EngineRunStatus::Failed { reason } => eprintln!("xTB failed: {reason}"),
//! }
//! ```
//!
//! ORCA has an input writer and a parser for its output, but no runner: the
//! module returns [`Error::InvalidInput`] instead of pretending to launch it.
//! Every completed run leaves the units of [`EngineResult`]: eV, eV/Angstrom and
//! GPa, converted from the engine's own output.
//!
//! **DFTB+ data.**  xTB is self-contained.  DFTB+ additionally needs its
//! Slater-Koster files under `slakos/` in the working directory, and the
//! runner creates a fresh scratch folder for every run, so it does not copy
//! them.  A real DFTB+ run therefore succeeds only when the installed binary
//! finds the files it needs; otherwise it is reported as
//! [`EngineRunStatus::Failed`] with DFTB+'s own message.  This is a documented
//! limitation of the pre-existing DFTB+ runner (Task 4.9), not of the
//! opportunistic path.

use std::sync::Arc;
use std::sync::atomic::AtomicBool;
use std::time::Duration;

use serde::{Deserialize, Serialize};
use xtalforge_core::Structure;

use crate::Error;
use crate::external::{
    DftbCalculator, DftbOptions, EngineKind, EngineResult, XtbCalculator, XtbOptions, probe,
    probe_binary,
};

/// Engines [`crate::probe`] knows about.
pub const PROBE_ENGINES: [EngineKind; 3] =
    [EngineKind::DftbPlus, EngineKind::Xtb, EngineKind::Orca];

/// Settings for one opportunistic external-engine run.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct EngineRunOptions {
    /// DFTB+ input options.
    pub dftb: DftbOptions,
    /// xTB input options.
    pub xtb: XtbOptions,
    /// Executable name or path; `None` uses the engine's standard name on
    /// `PATH`.
    pub program: Option<String>,
    /// Wall-clock limit in seconds.
    pub timeout_secs: f64,
    /// Keep the scratch folder after the run.
    pub keep_scratch: bool,
}

impl Default for EngineRunOptions {
    fn default() -> Self {
        Self {
            dftb: DftbOptions::default(),
            xtb: XtbOptions::default(),
            program: None,
            timeout_secs: 600.0,
            keep_scratch: false,
        }
    }
}

impl EngineRunOptions {
    /// Checks the timeout and the engine options.
    ///
    /// # Errors
    ///
    /// Returns [`Error::InvalidInput`] for a non-positive or non-finite
    /// timeout, or invalid DFTB+/xTB options.
    pub fn validate(&self) -> Result<(), Error> {
        if !self.timeout_secs.is_finite() || self.timeout_secs <= 0.0 {
            return Err(Error::InvalidInput(format!(
                "external-engine timeout must be positive and finite, got {}",
                self.timeout_secs
            )));
        }
        self.dftb.validate()?;
        self.xtb.validate()?;
        Ok(())
    }
}

/// What an opportunistic run did.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
#[serde(rename_all = "snake_case")]
pub enum EngineRunStatus {
    /// The binary is not installed; nothing was launched.
    Skipped {
        /// Why the engine was skipped; safe to show to a user.
        reason: String,
    },
    /// The binary ran and its output parsed.
    Completed(EngineResult),
    /// The binary is installed but the run failed.
    Failed {
        /// What went wrong, as a human-readable message.
        reason: String,
    },
}

impl EngineRunStatus {
    /// Whether the run was skipped because the binary is missing.
    #[must_use]
    pub fn is_skipped(&self) -> bool {
        matches!(self, Self::Skipped { .. })
    }

    /// Whether the run completed and produced a result.
    #[must_use]
    pub fn is_completed(&self) -> bool {
        matches!(self, Self::Completed(_))
    }

    /// Whether the binary exists but the run failed.
    #[must_use]
    pub fn is_failed(&self) -> bool {
        matches!(self, Self::Failed { .. })
    }

    /// The parsed result of a completed run.
    #[must_use]
    pub fn result(&self) -> Option<&EngineResult> {
        match self {
            Self::Completed(result) => Some(result),
            _ => None,
        }
    }

    /// The skip or failure reason.
    #[must_use]
    pub fn reason(&self) -> Option<&str> {
        match self {
            Self::Skipped { reason } | Self::Failed { reason } => Some(reason),
            Self::Completed(_) => None,
        }
    }
}

/// Probes every engine this build knows about, in [`PROBE_ENGINES`] order.
///
/// No process is launched for a missing binary: the probe only searches
/// `PATH` and runs `<binary> --version` when the file exists and is
/// executable.
#[must_use]
pub fn available_engines() -> Vec<(EngineKind, crate::external::EngineStatus)> {
    PROBE_ENGINES
        .iter()
        .map(|kind| (*kind, probe(*kind)))
        .collect()
}

/// Runs DFTB+ or xTB when its binary exists, and reports a skip when it does
/// not.
///
/// `cancel` is shared with the calculator, so setting it kills the running
/// process and the caller gets [`EngineRunStatus::Failed`].
///
/// # Errors
///
/// Returns [`Error::InvalidInput`] for invalid options or for ORCA, which has
/// no runner in this crate.  A missing binary and a failed run are statuses,
/// not errors.
pub fn try_run_engine(
    kind: EngineKind,
    structure: &Structure,
    options: &EngineRunOptions,
    cancel: Arc<AtomicBool>,
) -> Result<EngineRunStatus, Error> {
    options.validate()?;
    let status = match options.program.as_deref() {
        Some(program) => probe_binary(program),
        None => probe(kind),
    };
    if let crate::external::EngineStatus::Missing { reason } = status {
        return Ok(EngineRunStatus::Skipped { reason });
    }
    let timeout = Duration::from_secs_f64(options.timeout_secs);
    let run = match kind {
        EngineKind::DftbPlus => {
            let mut calculator = DftbCalculator::new(options.dftb.clone())
                .with_timeout(timeout)
                .keep_scratch(options.keep_scratch)
                .with_cancel(cancel);
            if let Some(program) = &options.program {
                calculator = calculator.with_program(program.clone());
            }
            calculator.run(structure)
        }
        EngineKind::Xtb => {
            let mut calculator = XtbCalculator::new(options.xtb.clone())
                .with_timeout(timeout)
                .keep_scratch(options.keep_scratch)
                .with_cancel(cancel);
            if let Some(program) = &options.program {
                calculator = calculator.with_program(program.clone());
            }
            calculator.run(structure)
        }
        EngineKind::Orca => {
            return Err(Error::InvalidInput(
                "ORCA has an input writer but no runner in this crate; write the input and launch ORCA yourself"
                    .to_owned(),
            ));
        }
    };
    Ok(match run {
        Ok(result) => EngineRunStatus::Completed(result),
        Err(error) => EngineRunStatus::Failed {
            reason: error.message(),
        },
    })
}
