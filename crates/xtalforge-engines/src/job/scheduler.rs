//! Scheduler submission, polling and collection behind a runtime probe
//! (Task 18.4.2).
//!
//! Submission scripts are always written. A dry run stops there (`prepared`
//! state, no scheduler call). A real submission probes the submit command,
//! runs it with the script path and parses the job id. Polling probes the
//! status command and maps the scheduler text to [`JobState`]; collection
//! reads the configured output and error files.

use std::path::{Path, PathBuf};

use serde::{Deserialize, Serialize};

use crate::Error;

use super::process::{probe_program, process_execution_available, run_command};
use super::script::write_scheduler_script;
use super::util::{validate_relative_path, validate_token};
use super::{JobOptions, JobState, Scheduler, SchedulerJob};

/// Maximum bytes returned from each collected scheduler log.
const MAX_COLLECTED_LOG_BYTES: usize = 1 << 20;

/// Everything needed to submit and poll one job.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct SchedulerConfig {
    /// Script and queue options.
    pub options: JobOptions,
    /// Submit command override (`sbatch`/`qsub` by default). Useful for fake
    /// schedulers in tests and for sites with a wrapper.
    pub submit_command: Option<String>,
    /// Status command override (`squeue`/`qstat` by default).
    pub status_command: Option<String>,
}

impl SchedulerConfig {
    /// A config that uses the default commands (or the
    /// `XTALFORGE_SBATCH`/`XTALFORGE_QSUB` environment overrides).
    #[must_use]
    pub fn new(options: JobOptions) -> Self {
        Self {
            options,
            submit_command: None,
            status_command: None,
        }
    }
}

/// Result of a runtime scheduler probe.
#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
pub struct SchedulerProbe {
    /// Scheduler that was probed.
    pub scheduler: Scheduler,
    /// Submit command that would be used.
    pub submit_command: String,
    /// Status command that would be used.
    pub status_command: String,
    /// Both commands were found and this build can launch them.
    pub available: bool,
    /// Human-readable reason, including the missing command when unavailable.
    pub reason: String,
}

/// A job accepted by (or prepared for) a scheduler.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct SubmittedJob {
    /// Job id and state.
    pub job: SchedulerJob,
    /// Submission script that was written.
    pub script: PathBuf,
    /// Whether the scheduler was not contacted.
    pub dry_run: bool,
}

/// Scheduler output collected after a poll.
#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
pub struct CollectedJob {
    /// Last observed job state.
    pub job: SchedulerJob,
    /// First 1 MiB of the standard-output file (empty when absent).
    pub stdout: String,
    /// First 1 MiB of the standard-error file (empty when absent).
    pub stderr: String,
}

/// Probes the submit and status commands without running a job.
///
/// The `XTALFORGE_SBATCH`, `XTALFORGE_SQUEUE`, `XTALFORGE_QSUB` and
/// `XTALFORGE_QSTAT` environment variables override the default commands.
#[must_use]
pub fn probe_scheduler(scheduler: Scheduler) -> SchedulerProbe {
    let Some((submit_command, status_command)) = probe_names(scheduler) else {
        return SchedulerProbe {
            scheduler,
            submit_command: String::new(),
            status_command: String::new(),
            available: false,
            reason: "unknown scheduler; expected slurm or pbs".to_owned(),
        };
    };
    if !process_execution_available() {
        return unavailable_probe(
            scheduler,
            submit_command,
            status_command,
            "process execution is disabled; rebuild with the `external` feature".to_owned(),
        );
    }
    if let Err(reason) = probe_program(&submit_command) {
        return unavailable_probe(
            scheduler,
            submit_command.clone(),
            status_command,
            format!("submit command {submit_command:?} is not available: {reason}"),
        );
    }
    if let Err(reason) = probe_program(&status_command) {
        return unavailable_probe(
            scheduler,
            submit_command,
            status_command.clone(),
            format!("status command {status_command:?} is not available: {reason}"),
        );
    }
    SchedulerProbe {
        scheduler,
        submit_command,
        status_command,
        available: true,
        reason: "scheduler commands are available".to_owned(),
    }
}

/// Writes the script and, unless this is a dry run, submits it.
///
/// # Errors
///
/// Returns [`Error::Job`] for invalid options or an unreadable job id,
/// [`Error::Unavailable`] for a missing scheduler command or a build without
/// the `external` feature, and [`Error::Job`] for a non-zero submit exit.
pub fn submit(config: &SchedulerConfig, work_dir: impl AsRef<Path>) -> Result<SubmittedJob, Error> {
    let script = write_scheduler_script(&config.options, work_dir)?;
    if config.options.dry_run {
        return Ok(SubmittedJob {
            job: SchedulerJob {
                id: "dry-run".to_owned(),
                state: JobState::Prepared,
            },
            script,
            dry_run: true,
        });
    }
    if !process_execution_available() {
        return Err(Error::Unavailable(
            "scheduler submission requires the `external` feature; the script was written but \
             not submitted"
                .to_owned(),
        ));
    }
    let program = resolve_command(config, true)?;
    probe_program(&program).map_err(|reason| {
        Error::Unavailable(format!(
            "scheduler submit command {program:?} is not available: {reason}"
        ))
    })?;
    let args = vec![script.to_string_lossy().into_owned()];
    let output = run_command(&program, &args)?;
    if !output.success {
        return Err(Error::Job(format!(
            "scheduler submission failed (exit {}): {}",
            exit_label(output.exit_code),
            output.stderr.trim()
        )));
    }
    let id = parse_job_id(config.options.scheduler, &output.stdout)?;
    Ok(SubmittedJob {
        job: SchedulerJob {
            id,
            state: JobState::Submitted,
        },
        script,
        dry_run: false,
    })
}

/// Queries the scheduler for the state of one job.
///
/// # Errors
///
/// Returns [`Error::Job`] for a dry-run config, an invalid id or a status
/// command that failed for an unexpected reason, and [`Error::Unavailable`]
/// for a missing command or a build without the `external` feature.
pub fn poll(config: &SchedulerConfig, job_id: &str) -> Result<SchedulerJob, Error> {
    validate_token(job_id, "job id")?;
    if config.options.dry_run {
        return Err(Error::Job(
            "a dry-run submission has no scheduler state".to_owned(),
        ));
    }
    if config.options.scheduler == Scheduler::Unknown {
        return Err(Error::Job(
            "unknown scheduler; expected slurm or pbs".to_owned(),
        ));
    }
    if !process_execution_available() {
        return Err(Error::Unavailable(
            "scheduler polling requires the `external` feature".to_owned(),
        ));
    }
    let program = resolve_command(config, false)?;
    probe_program(&program).map_err(|reason| {
        Error::Unavailable(format!(
            "scheduler status command {program:?} is not available: {reason}"
        ))
    })?;
    let args = match config.options.scheduler {
        Scheduler::Slurm => vec![
            "-h".to_owned(),
            "-j".to_owned(),
            job_id.to_owned(),
            "-o".to_owned(),
            "%T".to_owned(),
        ],
        Scheduler::Pbs => vec!["-f".to_owned(), job_id.to_owned()],
        Scheduler::Unknown => {
            return Err(Error::Job(
                "unknown scheduler; expected slurm or pbs".to_owned(),
            ));
        }
    };
    let output = run_command(&program, &args)?;
    if !output.success {
        let text = format!("{}\n{}", output.stdout, output.stderr);
        if text.trim().is_empty() || unknown_job(&text) {
            return Ok(SchedulerJob {
                id: job_id.to_owned(),
                state: JobState::Unknown,
            });
        }
        return Err(Error::Job(format!(
            "scheduler poll failed (exit {}): {}",
            exit_label(output.exit_code),
            text.trim()
        )));
    }
    Ok(SchedulerJob {
        id: job_id.to_owned(),
        state: parse_state(config.options.scheduler, &output.stdout),
    })
}

/// Reads the configured output and error files of a submitted job.
///
/// Missing files are reported as empty strings: not every scheduler creates an
/// error file for a job that never failed. Unreadable files are an error.
///
/// # Errors
///
/// Returns [`Error::Job`] for an unsafe path or an unreadable file.
pub fn collect(
    config: &SchedulerConfig,
    job: &SchedulerJob,
    work_dir: impl AsRef<Path>,
) -> Result<CollectedJob, Error> {
    validate_relative_path(&config.options.output, "output path")?;
    validate_relative_path(&config.options.error, "error path")?;
    let directory = work_dir.as_ref();
    let stdout = read_log(&directory.join(&config.options.output))?;
    let stderr = read_log(&directory.join(&config.options.error))?;
    Ok(CollectedJob {
        job: job.clone(),
        stdout: stdout.unwrap_or_default(),
        stderr: stderr.unwrap_or_default(),
    })
}

fn probe_names(scheduler: Scheduler) -> Option<(String, String)> {
    match scheduler {
        Scheduler::Slurm => Some((
            env_or("XTALFORGE_SBATCH", "sbatch"),
            env_or("XTALFORGE_SQUEUE", "squeue"),
        )),
        Scheduler::Pbs => Some((
            env_or("XTALFORGE_QSUB", "qsub"),
            env_or("XTALFORGE_QSTAT", "qstat"),
        )),
        Scheduler::Unknown => None,
    }
}

fn env_or(variable: &str, default: &str) -> String {
    std::env::var(variable)
        .ok()
        .filter(|value| !value.trim().is_empty())
        .unwrap_or_else(|| default.to_owned())
}

fn unavailable_probe(
    scheduler: Scheduler,
    submit_command: String,
    status_command: String,
    reason: String,
) -> SchedulerProbe {
    SchedulerProbe {
        scheduler,
        submit_command,
        status_command,
        available: false,
        reason,
    }
}

fn resolve_command(config: &SchedulerConfig, submit: bool) -> Result<String, Error> {
    let configured = if submit {
        config.submit_command.as_ref()
    } else {
        config.status_command.as_ref()
    };
    if let Some(program) = configured.filter(|program| !program.trim().is_empty()) {
        return Ok(program.clone());
    }
    let Some((default_submit, default_status)) = config.options.scheduler.commands() else {
        return Err(Error::Job(
            "unknown scheduler; expected slurm or pbs".to_owned(),
        ));
    };
    Ok(match config.options.scheduler {
        Scheduler::Slurm => {
            if submit {
                env_or("XTALFORGE_SBATCH", default_submit)
            } else {
                env_or("XTALFORGE_SQUEUE", default_status)
            }
        }
        Scheduler::Pbs => {
            if submit {
                env_or("XTALFORGE_QSUB", default_submit)
            } else {
                env_or("XTALFORGE_QSTAT", default_status)
            }
        }
        Scheduler::Unknown => {
            return Err(Error::Job(
                "unknown scheduler; expected slurm or pbs".to_owned(),
            ));
        }
    })
}

fn parse_job_id(scheduler: Scheduler, stdout: &str) -> Result<String, Error> {
    let tokens: Vec<&str> = stdout.split_whitespace().collect();
    let candidate = match scheduler {
        Scheduler::Slurm => tokens
            .iter()
            .position(|token| token.eq_ignore_ascii_case("job"))
            .and_then(|index| tokens.get(index.saturating_add(1)))
            .copied(),
        Scheduler::Pbs => tokens.first().copied(),
        Scheduler::Unknown => None,
    };
    let id = candidate
        .map(str::trim)
        .filter(|id| !id.is_empty())
        .ok_or_else(|| {
            Error::Job(format!(
                "scheduler did not report a job id (output {:?})",
                stdout.trim()
            ))
        })?;
    if id.chars().any(char::is_control) {
        return Err(Error::Job(
            "scheduler reported a job id with control characters".to_owned(),
        ));
    }
    Ok(id.to_owned())
}

fn parse_state(scheduler: Scheduler, text: &str) -> JobState {
    match scheduler {
        Scheduler::Slurm => parse_slurm_state(text),
        Scheduler::Pbs => parse_pbs_state(text),
        Scheduler::Unknown => JobState::Unknown,
    }
}

fn parse_slurm_state(text: &str) -> JobState {
    for token in text.split_whitespace() {
        let state = match token.to_ascii_uppercase().as_str() {
            "PENDING" | "CONFIGURING" | "REQUEUED" | "RESIZING" | "SUSPENDED" => {
                Some(JobState::Pending)
            }
            "RUNNING" | "COMPLETING" | "STAGE_OUT" => Some(JobState::Running),
            "COMPLETED" => Some(JobState::Completed),
            "FAILED" | "CANCELLED" | "TIMEOUT" | "OUT_OF_MEMORY" | "NODE_FAIL" | "PREEMPTED"
            | "BOOT_FAIL" | "DEADLINE" => Some(JobState::Failed),
            _ => None,
        };
        if let Some(state) = state {
            return state;
        }
    }
    JobState::Unknown
}

fn parse_pbs_state(text: &str) -> JobState {
    for line in text.lines() {
        let Some((key, value)) = line.split_once('=') else {
            continue;
        };
        if key.trim().eq_ignore_ascii_case("job_state")
            && let Some(state) = pbs_state(value.trim().chars().next())
        {
            return state;
        }
    }
    if !text.contains('=')
        && let Some(token) = text.split_whitespace().next()
        && let Some(state) = pbs_state(token.chars().next())
    {
        return state;
    }
    JobState::Unknown
}

fn pbs_state(value: Option<char>) -> Option<JobState> {
    match value.map(|character| character.to_ascii_uppercase()) {
        Some('Q' | 'H' | 'W' | 'T' | 'S') => Some(JobState::Pending),
        Some('R' | 'E' | 'B') => Some(JobState::Running),
        Some('C') => Some(JobState::Completed),
        Some('F') => Some(JobState::Failed),
        _ => None,
    }
}

fn unknown_job(text: &str) -> bool {
    let lower = text.to_ascii_lowercase();
    lower.contains("unknown job")
        || lower.contains("invalid job id")
        || lower.contains("invalid job")
        || lower.contains("no such job")
}

fn exit_label(code: Option<i32>) -> String {
    code.map_or_else(|| "signal".to_owned(), |value| value.to_string())
}

fn read_log(path: &Path) -> Result<Option<String>, Error> {
    match std::fs::read(path) {
        Ok(bytes) => {
            let bytes = bytes.get(..MAX_COLLECTED_LOG_BYTES).unwrap_or(&bytes);
            Ok(Some(String::from_utf8_lossy(bytes).into_owned()))
        }
        Err(error) if error.kind() == std::io::ErrorKind::NotFound => Ok(None),
        Err(error) => Err(Error::Job(format!(
            "job log {} could not be read: {error}",
            path.display()
        ))),
    }
}
