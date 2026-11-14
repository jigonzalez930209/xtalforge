//! Workflow `compute` steps: write, run, parse and checkpoint (Task 18.4.1).

use std::path::{Path, PathBuf};
use std::time::Duration;

use serde::{Deserialize, Serialize};
use xtalforge_core::Structure;

use crate::Error;
use crate::input_set::{InputFile, InputSetRequest, PresetOverrides, write_input_set};
use crate::result::CalculationResult;

use super::checkpoint::{
    JOB_CHECKPOINT_VERSION, JobCheckpoint, checkpoint_path, request_fingerprint,
};
use super::failure::detect_failure;
use super::inputs::InputFiles;
use super::parser::ParserRegistry;
use super::process::{EngineCommand, process_execution_available, run_process};
use super::util::validate_relative_path;
use super::{FailurePattern, JobState};

/// Output sentinel: parse the captured standard output instead of a file.
pub const STDOUT_OUTPUT: &str = "<stdout>";

/// A typed workflow `compute` step.
///
/// The preset and overrides are resolved by the crate's input-set registry;
/// the command, parser name and output file describe the engine run. The step
/// itself is serializable so a workflow can be stored and resumed.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct ComputeStep {
    /// Job name, also used for the scratch directory and checkpoint file.
    pub name: String,
    /// Input-set preset name.
    pub preset: String,
    /// Overrides merged into the preset.
    pub overrides: PresetOverrides,
    /// Engine command.
    pub command: EngineCommand,
    /// Parser name resolved by [`ParserRegistry`].
    pub parser: String,
    /// Output file read by the parser, or [`STDOUT_OUTPUT`].
    pub output: String,
}

impl ComputeStep {
    /// Builds a compute step.
    #[must_use]
    pub fn new(
        name: impl Into<String>,
        preset: impl Into<String>,
        command: EngineCommand,
        parser: impl Into<String>,
        output: impl Into<String>,
    ) -> Self {
        Self {
            name: name.into(),
            preset: preset.into(),
            overrides: PresetOverrides::new(),
            command,
            parser: parser.into(),
            output: output.into(),
        }
    }

    /// Builds the request through the input-set registry.
    ///
    /// # Errors
    ///
    /// Returns [`Error::UnknownPreset`] when the preset is not registered and
    /// any error raised by the preset writer.
    pub fn request(&self, structure: &Structure) -> Result<JobRequest, Error> {
        let request = InputSetRequest {
            preset: self.preset.clone(),
            overrides: self.overrides.clone(),
        };
        let summary = write_input_set(structure, &request)?;
        Ok(self.request_with_inputs(summary.files))
    }

    /// Builds the request from already materialized files.
    ///
    /// This is the integration point for engine writers that are not yet in
    /// the input-set registry and for callers that wrote the files themselves.
    #[must_use]
    pub fn request_with_inputs(&self, files: Vec<InputFile>) -> JobRequest {
        JobRequest {
            name: self.name.clone(),
            inputs: InputFiles::new(files),
            command: self.command.clone(),
            parser: self.parser.clone(),
            output: self.output.clone(),
        }
    }
}

/// One compute job ready to run.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct JobRequest {
    /// Job name.
    pub name: String,
    /// Input files written into the scratch directory.
    pub inputs: InputFiles,
    /// Engine command run with the scratch directory as working directory.
    pub command: EngineCommand,
    /// Parser name resolved by [`ParserRegistry`].
    pub parser: String,
    /// Output file read by the parser, or [`STDOUT_OUTPUT`].
    pub output: String,
}

impl JobRequest {
    /// Rejects unsafe names and an empty command or parser.
    ///
    /// # Errors
    ///
    /// Returns [`Error::Job`] for an invalid request.
    pub fn validate(&self) -> Result<(), Error> {
        validate_relative_path(&self.name, "job name")?;
        self.inputs.validate()?;
        self.command.validate()?;
        if self.parser.trim().is_empty() {
            return Err(Error::Job("a compute job needs a parser name".to_owned()));
        }
        if self.output != STDOUT_OUTPUT {
            validate_relative_path(&self.output, "output file")?;
        }
        Ok(())
    }
}

/// Runtime options of one compute run.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct RunOptions {
    /// Directory that holds the scratch folder and the checkpoint.
    pub work_dir: PathBuf,
    /// Write the inputs only; never launch the engine or touch a checkpoint.
    pub dry_run: bool,
    /// Reuse a completed checkpoint and restart a failed one.
    pub resume: bool,
    /// Wall-clock limit for the engine process.
    pub timeout: Duration,
    /// Maximum bytes kept from each captured log.
    pub max_log_bytes: usize,
}

impl RunOptions {
    /// Defaults: resume enabled, a one-hour timeout and a 1 MiB log cap.
    #[must_use]
    pub fn new(work_dir: impl Into<PathBuf>) -> Self {
        Self {
            work_dir: work_dir.into(),
            dry_run: false,
            resume: true,
            timeout: Duration::from_secs(3600),
            max_log_bytes: 1 << 20,
        }
    }
}

/// Typed outcome of a compute run.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct JobOutcome {
    /// Job name.
    pub name: String,
    /// Final state.
    pub state: JobState,
    /// Engine launches attempted (zero for a dry run).
    pub attempts: u32,
    /// Exit code when the platform reported one.
    pub exit_code: Option<i32>,
    /// Captured standard output.
    pub stdout: String,
    /// Captured standard error.
    pub stderr: String,
    /// Parsed result, when the run produced one.
    pub result: Option<CalculationResult>,
    /// Detected failure pattern with its corrective action.
    pub failure: Option<FailurePattern>,
    /// Whether this was a dry run.
    pub dry_run: bool,
    /// Whether a completed checkpoint supplied the outcome.
    pub resumed: bool,
    /// Scratch directory used for the run.
    pub scratch: PathBuf,
}

/// Writes the inputs, runs the engine and parses the result.
///
/// A non-zero engine exit is not a Rust error: the outcome is `failed` and
/// [`JobOutcome::failure`] carries the corrective action when a known pattern
/// matched. A launch failure, a missing output file, an unknown parser or
/// malformed output returns the corresponding [`Error`] after recording a
/// failed checkpoint.
///
/// # Errors
///
/// Returns [`Error::Job`] for an invalid request, options or checkpoint,
/// [`Error::Unavailable`] when the build has no `external` feature (the inputs
/// and a `prepared` checkpoint are written first), [`Error::Missing`] for a
/// missing output file and the parser error for malformed output.
pub fn run_job(
    request: &JobRequest,
    registry: &ParserRegistry,
    options: &RunOptions,
) -> Result<JobOutcome, Error> {
    request.validate()?;
    if options.timeout.is_zero() {
        return Err(Error::Job("the job timeout must be positive".to_owned()));
    }
    if options.max_log_bytes == 0 {
        return Err(Error::Job("the job log cap must be positive".to_owned()));
    }
    let scratch = options.work_dir.join(&request.name);
    std::fs::create_dir_all(&scratch).map_err(|error| {
        Error::Job(format!(
            "job scratch directory {} could not be created: {error}",
            scratch.display()
        ))
    })?;
    request.inputs.write_to(&scratch)?;

    if options.dry_run {
        return Ok(JobOutcome {
            name: request.name.clone(),
            state: JobState::Prepared,
            attempts: 0,
            exit_code: None,
            stdout: String::new(),
            stderr: String::new(),
            result: None,
            failure: None,
            dry_run: true,
            resumed: false,
            scratch,
        });
    }

    let checkpoint_file = checkpoint_path(&options.work_dir, &request.name)?;
    let fingerprint = request_fingerprint(request)?;
    let previous = if options.resume {
        JobCheckpoint::read_if_exists(&checkpoint_file)?
    } else {
        None
    };
    if let Some(checkpoint) = &previous {
        if checkpoint.checkpoint_version != JOB_CHECKPOINT_VERSION {
            return Err(Error::Job(format!(
                "job checkpoint {} uses version {}; expected {}",
                checkpoint_file.display(),
                checkpoint.checkpoint_version,
                JOB_CHECKPOINT_VERSION
            )));
        }
        if checkpoint.fingerprint != fingerprint {
            return Err(Error::Job(format!(
                "job checkpoint {} belongs to a different request; \
                 use a fresh work directory or disable resume",
                checkpoint_file.display()
            )));
        }
        if checkpoint.state == JobState::Completed {
            return Ok(JobOutcome {
                name: request.name.clone(),
                state: JobState::Completed,
                attempts: checkpoint.attempts,
                exit_code: checkpoint.exit_code,
                stdout: checkpoint.stdout.clone(),
                stderr: checkpoint.stderr.clone(),
                result: checkpoint.result.clone(),
                failure: checkpoint.failure.clone(),
                dry_run: false,
                resumed: true,
                scratch,
            });
        }
    }

    if !process_execution_available() {
        let attempts = previous
            .as_ref()
            .map_or(0, |checkpoint| checkpoint.attempts);
        let checkpoint =
            JobCheckpoint::new(&request.name, &fingerprint, JobState::Prepared, attempts);
        checkpoint.write(&checkpoint_file)?;
        return Err(Error::Unavailable(format!(
            "engine process execution requires the `external` feature; \
             the inputs and a prepared checkpoint were written to {}",
            scratch.display()
        )));
    }

    let attempts = if options.resume {
        previous
            .as_ref()
            .map_or(1, |checkpoint| checkpoint.attempts.saturating_add(1))
    } else {
        1
    };
    let mut running = JobCheckpoint::new(&request.name, &fingerprint, JobState::Running, attempts);
    running.write(&checkpoint_file)?;
    let process = match run_process(
        &request.command,
        &scratch,
        options.timeout,
        options.max_log_bytes,
    ) {
        Ok(process) => process,
        Err(error) => {
            running.state = JobState::Failed;
            running.error = Some(error.to_string());
            running.write(&checkpoint_file)?;
            return Err(error);
        }
    };

    let mut outcome = JobOutcome {
        name: request.name.clone(),
        state: JobState::Unknown,
        attempts,
        exit_code: process.exit_code,
        stdout: process.stdout,
        stderr: process.stderr,
        result: None,
        failure: None,
        dry_run: false,
        resumed: false,
        scratch,
    };
    if process.timed_out {
        outcome.state = JobState::Failed;
        outcome.failure = FailurePattern::from_code("walltime_exceeded");
    } else if !process.success {
        outcome.state = JobState::Failed;
        let combined = format!("{}\n{}", outcome.stdout, outcome.stderr);
        outcome.failure = detect_failure(&combined);
    } else {
        let output = if request.output == STDOUT_OUTPUT {
            Ok(outcome.stdout.clone())
        } else {
            let path = outcome.scratch.join(&request.output);
            std::fs::read_to_string(&path).map_err(|error| {
                if error.kind() == std::io::ErrorKind::NotFound {
                    Error::Missing(format!(
                        "engine output file {} was not produced",
                        path.display()
                    ))
                } else {
                    Error::Job(format!(
                        "engine output file {} could not be read: {error}",
                        path.display()
                    ))
                }
            })
        };
        let text = match output {
            Ok(text) => text,
            Err(error) => {
                outcome.state = JobState::Failed;
                save_checkpoint(
                    &checkpoint_file,
                    request,
                    &fingerprint,
                    &outcome,
                    Some(error.to_string()),
                )?;
                return Err(error);
            }
        };
        let combined = format!("{}\n{}\n{}", outcome.stdout, outcome.stderr, text);
        let failure = detect_failure(&combined);
        match registry.parse(&request.parser, &text) {
            Ok(result) => {
                outcome.result = Some(result);
                if let Some(pattern) = failure {
                    outcome.state = JobState::Failed;
                    outcome.failure = Some(pattern);
                } else {
                    outcome.state = JobState::Completed;
                }
            }
            Err(error) => {
                outcome.state = JobState::Failed;
                save_checkpoint(
                    &checkpoint_file,
                    request,
                    &fingerprint,
                    &outcome,
                    Some(error.to_string()),
                )?;
                return Err(error);
            }
        }
    }
    save_checkpoint(&checkpoint_file, request, &fingerprint, &outcome, None)?;
    Ok(outcome)
}

fn save_checkpoint(
    path: &Path,
    request: &JobRequest,
    fingerprint: &str,
    outcome: &JobOutcome,
    error: Option<String>,
) -> Result<(), Error> {
    let mut checkpoint =
        JobCheckpoint::new(&request.name, fingerprint, outcome.state, outcome.attempts);
    checkpoint.exit_code = outcome.exit_code;
    checkpoint.stdout = outcome.stdout.clone();
    checkpoint.stderr = outcome.stderr.clone();
    checkpoint.result = outcome.result.clone();
    checkpoint.failure = outcome.failure.clone();
    checkpoint.error = error;
    checkpoint.write(path)
}
