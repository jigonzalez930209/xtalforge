//! Job execution (Task 18.4).
//!
//! This module is the host-side half of a workflow `compute` step:
//!
//! * [`ComputeStep`] turns a structure into an [`InputFiles`] set through the
//!   crate's input-set registry and wraps it in a typed [`JobRequest`].
//!   [`run_job`] writes the files into a deterministic scratch directory, runs
//!   the configured engine command and hands the engine output to a parser
//!   registered in [`ParserRegistry`]. Process launch is only compiled in with
//!   the crate's `external` feature; without it a real run returns
//!   [`Error::Unavailable`] after the inputs have been written, while a dry run
//!   still writes the inputs and returns a `prepared` outcome.
//! * [`JobCheckpoint`] records the typed state after every stage, so a failed
//!   or interrupted run can be restarted and a completed run can be resumed
//!   without launching the engine again. A checkpoint belongs to one request
//!   fingerprint; a mismatch is refused instead of mixing results.
//! * [`scheduler_script`] writes deterministic SLURM and PBS submission
//!   scripts. [`submit`], [`poll`] and [`collect`] drive them behind a runtime
//!   probe and accept command overrides for fake schedulers in tests. A
//!   dry-run submission writes only the script.
//! * [`detect_failure`] maps the text of a failed (or misleadingly successful)
//!   run to a documented corrective action. The full catalogue is available
//!   through [`failure_catalogue`] and [`corrective_action`].
//!
//! Process execution and scheduler probes are host-only; WASM must use the
//! pure-Rust writing and parsing half. That exception is recorded in the
//! parity matrix by Task 18.5.1.
//!
//! # Known failure patterns
//!
//! Detection is case-insensitive and returns the first rule that matches, in
//! this order: walltime, out of memory, missing pseudopotential, `ZBRENT`,
//! SCF convergence. Each rule carries a stable code and a corrective action:
//!
//! | Code | Typical text | Corrective action |
//! |---|---|---|
//! | `walltime_exceeded` | `DUE TO TIME LIMIT` | resubmit with more walltime and restart from the saved state |
//! | `out_of_memory` | `oom-kill`, `Cannot allocate memory` | fewer ranks or more memory per node |
//! | `missing_pseudopotential` | `missing in the POTCAR`, `UPF file not found` | install or point at the pseudopotential library |
//! | `zbrent` | `ZBRENT` | smaller `POTIM`, restart from `CONTCAR`, gentler optimizer |
//! | `scf_not_converged` | `SCF not converged` | better mixing and restart from the saved density |

mod checkpoint;
mod compute;
mod failure;
mod inputs;
mod parser;
mod process;
mod scheduler;
mod script;
mod util;

use serde::{Deserialize, Serialize};

use crate::Error;

pub use checkpoint::{JOB_CHECKPOINT_VERSION, JobCheckpoint, checkpoint_path, request_fingerprint};
pub use compute::{ComputeStep, JobOutcome, JobRequest, RunOptions, STDOUT_OUTPUT, run_job};
pub use failure::{corrective_action, detect_failure, failure_catalogue};
pub use inputs::InputFiles;
pub use parser::{ParserFn, ParserRegistry};
pub use process::{EngineCommand, process_execution_available};
pub use scheduler::{
    CollectedJob, SchedulerConfig, SchedulerProbe, SubmittedJob, collect, poll, probe_scheduler,
    submit,
};

/// Supported schedulers.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
#[serde(rename_all = "snake_case")]
pub enum Scheduler {
    /// SLURM.
    Slurm,
    /// PBS/Torque.
    Pbs,
    /// A scheduler this build does not know; every script or probe request is
    /// refused with [`Error::Job`].
    Unknown,
}

impl Scheduler {
    /// File extension used for the submission script, or `None` for an unknown
    /// scheduler.
    #[must_use]
    pub const fn extension(self) -> Option<&'static str> {
        match self {
            Self::Slurm => Some("slurm"),
            Self::Pbs => Some("pbs"),
            Self::Unknown => None,
        }
    }

    /// Default submit and status commands, or `None` for an unknown scheduler.
    #[must_use]
    pub const fn commands(self) -> Option<(&'static str, &'static str)> {
        match self {
            Self::Slurm => Some(("sbatch", "squeue")),
            Self::Pbs => Some(("qsub", "qstat")),
            Self::Unknown => None,
        }
    }
}

/// Job options for one submission script.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct JobOptions {
    /// Scheduler.
    pub scheduler: Scheduler,
    /// Queue/partition.
    pub queue: String,
    /// Nodes.
    pub nodes: usize,
    /// Walltime (`HH:MM:SS` or `D-HH:MM:SS`).
    pub walltime: String,
    /// Command executed by the scheduler.
    pub command: String,
    /// Dry run only: write the script, never contact the scheduler.
    pub dry_run: bool,
    /// Job name used in the script and in the default file names.
    pub job_name: String,
    /// Scheduler standard-output path, relative to the submission directory.
    pub output: String,
    /// Scheduler standard-error path, relative to the submission directory.
    pub error: String,
}

impl JobOptions {
    /// Options with a single node, a one-hour walltime and `job.out`/`job.err`
    /// output paths.
    #[must_use]
    pub fn new(scheduler: Scheduler, command: impl Into<String>) -> Self {
        Self {
            scheduler,
            queue: "main".to_owned(),
            nodes: 1,
            walltime: "01:00:00".to_owned(),
            command: command.into(),
            dry_run: false,
            job_name: "xtalforge".to_owned(),
            output: "job.out".to_owned(),
            error: "job.err".to_owned(),
        }
    }
}

/// State of a compute job or scheduler submission.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
#[serde(rename_all = "snake_case")]
pub enum JobState {
    /// Inputs were written but nothing was executed (dry run or missing
    /// external feature).
    Prepared,
    /// Pending in a scheduler queue.
    Pending,
    /// Running.
    Running,
    /// Accepted by a scheduler but not yet seen in a status query.
    Submitted,
    /// Completed successfully.
    Completed,
    /// Failed; [`JobOutcome::failure`] may carry the corrective action.
    Failed,
    /// Cancelled by the scheduler or the user.
    Cancelled,
    /// The scheduler no longer knows the job, or the state was not understood.
    Unknown,
}

impl JobState {
    /// Whether no further transition is expected without a new submission.
    #[must_use]
    pub const fn is_terminal(self) -> bool {
        matches!(self, Self::Completed | Self::Failed | Self::Cancelled)
    }

    /// Whether the job completed successfully.
    #[must_use]
    pub const fn is_success(self) -> bool {
        matches!(self, Self::Completed)
    }
}

/// One scheduler job.
#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
pub struct SchedulerJob {
    /// Job id.
    pub id: String,
    /// State.
    pub state: JobState,
}

/// A known failure pattern with its documented corrective action.
#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
pub struct FailurePattern {
    /// Stable code (`scf_not_converged`, `zbrent`, `walltime_exceeded`,
    /// `out_of_memory`, `missing_pseudopotential`).
    pub code: String,
    /// Message fragment that matched, lower-case.
    pub pattern: String,
    /// Corrective action to show to the user.
    pub action: String,
}

impl FailurePattern {
    /// Builds the canonical pattern for a code, when the code is known.
    pub(crate) fn from_code(code: &str) -> Option<Self> {
        failure::pattern_for(code)
    }
}

/// Builds a deterministic submission script.
///
/// The command is emitted verbatim after the scheduler directives; it is the
/// user's shell line. Job name, queue, walltime and output paths must be single
/// tokens so they cannot inject extra directives.
///
/// # Errors
///
/// Returns [`Error::Job`] for an unknown scheduler, an empty command, a
/// non-positive node count or an invalid directive value.
pub fn scheduler_script(options: &JobOptions) -> Result<String, Error> {
    script::scheduler_script(options)
}

/// Writes a submission script into `work_dir` and returns its path.
///
/// The file is named `<job_name>.<slurm|pbs>` and is written atomically.
///
/// # Errors
///
/// Returns [`Error::Job`] for invalid options or an unwritable directory.
pub fn write_scheduler_script(
    options: &JobOptions,
    work_dir: impl AsRef<std::path::Path>,
) -> Result<std::path::PathBuf, Error> {
    script::write_scheduler_script(options, work_dir)
}
