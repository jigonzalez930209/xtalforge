//! Running DFTB+ and xTB through the scratch-folder runner (Task 4.9).

use std::collections::BTreeMap;
use std::sync::Arc;
use std::sync::atomic::AtomicBool;
use std::time::Duration;

use xtalforge_core::Structure;

use crate::Error;
use crate::calculator::Calculator;

use super::{
    DftbOptions, EngineResult, RunSpec, RunStatus, XtbOptions, parse_dftb_output, parse_xtb_output,
    run_in_scratch, write_dftb_hsd, write_xtb_input,
};

const DEFAULT_TIMEOUT: Duration = Duration::from_secs(600);

/// A DFTB+ calculator that launches `dftb+` in a scratch folder.
#[derive(Debug, Clone)]
pub struct DftbCalculator {
    options: DftbOptions,
    program: String,
    timeout: Duration,
    keep_scratch: bool,
    cancel: Arc<AtomicBool>,
}

impl DftbCalculator {
    /// Creates a calculator using `dftb+` from `PATH`.
    #[must_use]
    pub fn new(options: DftbOptions) -> Self {
        Self {
            options,
            program: "dftb+".to_owned(),
            timeout: DEFAULT_TIMEOUT,
            keep_scratch: false,
            cancel: Arc::new(AtomicBool::new(false)),
        }
    }

    /// Uses a different executable name or path.
    #[must_use]
    pub fn with_program(mut self, program: impl Into<String>) -> Self {
        self.program = program.into();
        self
    }

    /// Sets the wall-clock timeout.
    #[must_use]
    pub fn with_timeout(mut self, timeout: Duration) -> Self {
        self.timeout = timeout;
        self
    }

    /// Keeps the scratch folder for inspection.
    #[must_use]
    pub fn keep_scratch(mut self, keep: bool) -> Self {
        self.keep_scratch = keep;
        self
    }

    /// Handle another thread can use to cancel the next run.
    #[must_use]
    pub fn cancel_handle(&self) -> Arc<AtomicBool> {
        Arc::clone(&self.cancel)
    }

    /// Uses a shared cancellation flag instead of the private one.
    ///
    /// A [`crate::CancelToken`] exposes its flag with
    /// [`crate::CancelToken::shared`].
    #[must_use]
    pub fn with_cancel(mut self, cancel: Arc<AtomicBool>) -> Self {
        self.cancel = cancel;
        self
    }

    /// Runs DFTB+ and parses `detailed.out`.
    ///
    /// # Errors
    ///
    /// Returns an error for invalid options, a failed launch, a non-zero exit,
    /// missing output or a parse failure.
    pub fn run(&self, structure: &Structure) -> Result<EngineResult, Error> {
        let input = write_dftb_hsd(structure, &self.options)?;
        let mut inputs = BTreeMap::new();
        inputs.insert("dftb_in.hsd".to_owned(), input);
        let mut spec = RunSpec::new(self.program.clone(), Vec::new());
        spec.input_files = inputs;
        spec.timeout = self.timeout;
        spec.keep_scratch = true;
        let outcome = run_in_scratch(&spec, &self.cancel)?;
        if outcome.status != RunStatus::Success {
            if !self.keep_scratch {
                let _ = std::fs::remove_dir_all(&outcome.scratch);
            }
            return Err(Error::InvalidInput(format!(
                "DFTB+ ended as {:?} (exit {:?}): {}",
                outcome.status,
                outcome.exit_code,
                outcome.stderr.trim()
            )));
        }
        let output = std::fs::read_to_string(outcome.scratch.join("detailed.out"));
        if !self.keep_scratch {
            let _ = std::fs::remove_dir_all(&outcome.scratch);
        }
        let output = output.map_err(|error| {
            Error::InvalidInput(format!(
                "DFTB+ did not write detailed.out ({}): {error}",
                outcome.stderr.trim()
            ))
        })?;
        parse_dftb_output(&output)
    }
}

impl Calculator for DftbCalculator {
    fn name(&self) -> &'static str {
        "dftb+"
    }

    fn energy(&self, structure: &Structure) -> Result<f64, Error> {
        Ok(self.run(structure)?.energy)
    }

    fn energy_forces(&self, structure: &Structure) -> Result<(f64, Vec<[f64; 3]>), Error> {
        let result = self.run(structure)?;
        let expected = structure.p1()?.n_atoms();
        if result.forces.len() != expected {
            return Err(Error::InvalidInput(format!(
                "DFTB+ returned {} force rows for {expected} atoms",
                result.forces.len()
            )));
        }
        Ok((result.energy, result.forces))
    }

    fn stress(&self, structure: &Structure) -> Result<Option<[[f64; 3]; 3]>, Error> {
        Ok(self.run(structure)?.stress)
    }
}

/// An xTB calculator that launches `xtb` in a scratch folder.
#[derive(Debug, Clone)]
pub struct XtbCalculator {
    options: XtbOptions,
    program: String,
    timeout: Duration,
    keep_scratch: bool,
    cancel: Arc<AtomicBool>,
}

impl XtbCalculator {
    /// Creates a calculator using `xtb` from `PATH`.
    #[must_use]
    pub fn new(options: XtbOptions) -> Self {
        Self {
            options,
            program: "xtb".to_owned(),
            timeout: DEFAULT_TIMEOUT,
            keep_scratch: false,
            cancel: Arc::new(AtomicBool::new(false)),
        }
    }

    /// Uses a different executable name or path.
    #[must_use]
    pub fn with_program(mut self, program: impl Into<String>) -> Self {
        self.program = program.into();
        self
    }

    /// Sets the wall-clock timeout.
    #[must_use]
    pub fn with_timeout(mut self, timeout: Duration) -> Self {
        self.timeout = timeout;
        self
    }

    /// Keeps the scratch folder for inspection.
    #[must_use]
    pub fn keep_scratch(mut self, keep: bool) -> Self {
        self.keep_scratch = keep;
        self
    }

    /// Handle another thread can use to cancel the next run.
    #[must_use]
    pub fn cancel_handle(&self) -> Arc<AtomicBool> {
        Arc::clone(&self.cancel)
    }

    /// Uses a shared cancellation flag instead of the private one.
    ///
    /// A [`crate::CancelToken`] exposes its flag with
    /// [`crate::CancelToken::shared`].
    #[must_use]
    pub fn with_cancel(mut self, cancel: Arc<AtomicBool>) -> Self {
        self.cancel = cancel;
        self
    }

    /// Runs xTB and parses its standard output.
    ///
    /// # Errors
    ///
    /// Returns an error for invalid options, an unsupported method, a failed
    /// launch, a non-zero exit or a parse failure.
    pub fn run(&self, structure: &Structure) -> Result<EngineResult, Error> {
        let input = write_xtb_input(structure, &self.options)?;
        let mut args = vec!["coord.xyz".to_owned()];
        args.extend(method_args(&self.options.method)?);
        args.push("--chrg".to_owned());
        args.push(format!("{:.6}", self.options.charge));
        args.push("--uhf".to_owned());
        args.push(self.options.uhf.to_string());
        let mut inputs = BTreeMap::new();
        inputs.insert("coord.xyz".to_owned(), input);
        let mut spec = RunSpec::new(self.program.clone(), args);
        spec.input_files = inputs;
        spec.timeout = self.timeout;
        spec.keep_scratch = self.keep_scratch;
        let outcome = run_in_scratch(&spec, &self.cancel)?;
        if outcome.status != RunStatus::Success {
            return Err(Error::InvalidInput(format!(
                "xTB ended as {:?} (exit {:?}): {}",
                outcome.status,
                outcome.exit_code,
                outcome.stderr.trim()
            )));
        }
        parse_xtb_output(&outcome.stdout)
    }
}

impl Calculator for XtbCalculator {
    fn name(&self) -> &'static str {
        "xtb"
    }

    fn energy(&self, structure: &Structure) -> Result<f64, Error> {
        Ok(self.run(structure)?.energy)
    }

    fn energy_forces(&self, structure: &Structure) -> Result<(f64, Vec<[f64; 3]>), Error> {
        let result = self.run(structure)?;
        let expected = structure.p1()?.n_atoms();
        if result.forces.len() != expected {
            return Err(Error::InvalidInput(format!(
                "xTB returned {} force rows for {expected} atoms",
                result.forces.len()
            )));
        }
        Ok((result.energy, result.forces))
    }
}

fn method_args(method: &str) -> Result<Vec<String>, Error> {
    match method.to_ascii_lowercase().as_str() {
        "gfn1" | "gfn1-xtb" => Ok(vec!["--gfn".to_owned(), "1".to_owned()]),
        "gfn2" | "gfn2-xtb" => Ok(vec!["--gfn".to_owned(), "2".to_owned()]),
        "gfn-ff" | "gfnff" => Ok(vec!["--gfnff".to_owned()]),
        other => Err(Error::InvalidInput(format!(
            "unsupported xTB method {other:?}; use GFN1-xTB, GFN2-xTB or GFN-FF"
        ))),
    }
}
