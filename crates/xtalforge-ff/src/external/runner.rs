//! Scratch-folder process runner with timeout, cancellation and logs.

use std::collections::BTreeMap;
use std::fs::{self, File};
use std::path::{Path, PathBuf};
use std::process::{Command, Stdio};
use std::sync::atomic::{AtomicBool, Ordering};
use std::time::{Duration, Instant, SystemTime, UNIX_EPOCH};

use serde::{Deserialize, Serialize};

use crate::Error;

const DEFAULT_MAX_LOG_BYTES: usize = 1 << 20;
const DEFAULT_POLL: Duration = Duration::from_millis(10);

/// One external program invocation.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct RunSpec {
    /// Executable name or absolute path.
    pub program: String,
    /// Command-line arguments.
    pub args: Vec<String>,
    /// Files written into the scratch folder before launch, relative names.
    pub input_files: BTreeMap<String, String>,
    /// Wall-clock limit.
    pub timeout: Duration,
    /// Keep the scratch folder after the run.
    pub keep_scratch: bool,
    /// Maximum bytes kept from each captured log.
    pub max_log_bytes: usize,
}

impl RunSpec {
    /// A run with a 60 s timeout and a 1 MiB log cap.
    #[must_use]
    pub fn new(program: impl Into<String>, args: Vec<String>) -> Self {
        Self {
            program: program.into(),
            args,
            input_files: BTreeMap::new(),
            timeout: Duration::from_secs(60),
            keep_scratch: false,
            max_log_bytes: DEFAULT_MAX_LOG_BYTES,
        }
    }

    /// Adds one relative input file.
    #[must_use]
    pub fn with_input(mut self, name: impl Into<String>, contents: impl Into<String>) -> Self {
        self.input_files.insert(name.into(), contents.into());
        self
    }
}

/// How an external run ended.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
#[serde(rename_all = "snake_case")]
pub enum RunStatus {
    /// Exit status zero.
    Success,
    /// The process exited non-zero.
    Failed,
    /// The wall-clock timeout expired and the process was killed.
    TimedOut,
    /// The cancellation flag was set and the process was killed.
    Cancelled,
}

/// Captured result of one external run.
#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct RunOutcome {
    /// How the process ended.
    pub status: RunStatus,
    /// Exit code when the platform reported one.
    pub exit_code: Option<i32>,
    /// Captured standard output.
    pub stdout: String,
    /// Captured standard error.
    pub stderr: String,
    /// Whether stdout was capped at `max_log_bytes`.
    pub stdout_truncated: bool,
    /// Whether stderr was capped at `max_log_bytes`.
    pub stderr_truncated: bool,
    /// Scratch folder used for the run (removed unless `keep_scratch`).
    pub scratch: PathBuf,
    /// Wall-clock duration.
    pub duration: Duration,
}

/// Runs `spec` in a private scratch folder.
///
/// Inputs are written before launch; stdout and stderr go to files so a chatty
/// engine cannot deadlock the poll loop. The child is killed when `cancel` is
/// set or the timeout expires, and the logs are still returned.
///
/// # Errors
///
/// Returns [`Error::InvalidInput`] for an unusable spec, an unsafe input name,
/// a filesystem failure, or a process that cannot be spawned.
pub fn run_in_scratch(spec: &RunSpec, cancel: &AtomicBool) -> Result<RunOutcome, Error> {
    validate(spec)?;
    let scratch = create_scratch()?;
    let result = run_in(&scratch, spec, cancel);
    if !spec.keep_scratch {
        let _ = fs::remove_dir_all(&scratch);
    }
    result
}

fn validate(spec: &RunSpec) -> Result<(), Error> {
    if spec.program.trim().is_empty() {
        return Err(Error::InvalidInput(
            "external program name must not be empty".to_owned(),
        ));
    }
    if spec.timeout.is_zero() {
        return Err(Error::InvalidInput(
            "external run timeout must be positive".to_owned(),
        ));
    }
    for name in spec.input_files.keys() {
        let path = Path::new(name);
        if path.is_absolute()
            || path.components().any(|component| {
                matches!(
                    component,
                    std::path::Component::ParentDir
                        | std::path::Component::RootDir
                        | std::path::Component::Prefix(_)
                )
            })
        {
            return Err(Error::InvalidInput(format!(
                "external input path {name:?} must be relative and stay in the scratch folder"
            )));
        }
    }
    Ok(())
}

fn run_in(scratch: &Path, spec: &RunSpec, cancel: &AtomicBool) -> Result<RunOutcome, Error> {
    for (name, contents) in &spec.input_files {
        let path = scratch.join(name);
        if let Some(parent) = path.parent() {
            fs::create_dir_all(parent).map_err(|error| {
                Error::InvalidInput(format!("external input {name:?}: {error}"))
            })?;
        }
        fs::write(&path, contents)
            .map_err(|error| Error::InvalidInput(format!("external input {name:?}: {error}")))?;
    }
    let stdout_path = scratch.join("stdout.log");
    let stderr_path = scratch.join("stderr.log");
    let mut child = loop {
        let stdout = File::create(&stdout_path)
            .map_err(|error| Error::InvalidInput(format!("external stdout log: {error}")))?;
        let stderr = File::create(&stderr_path)
            .map_err(|error| Error::InvalidInput(format!("external stderr log: {error}")))?;
        match Command::new(&spec.program)
            .args(&spec.args)
            .current_dir(scratch)
            .stdin(Stdio::null())
            .stdout(Stdio::from(stdout))
            .stderr(Stdio::from(stderr))
            .spawn()
        {
            Ok(child) => break child,
            Err(error) if error.raw_os_error() == Some(26) => {
                // ETXTBSY: another thread forked while the script was being
                // written. Retry briefly; it clears as soon as the child exits.
                std::thread::sleep(Duration::from_millis(10));
            }
            Err(error) => {
                return Err(Error::InvalidInput(format!(
                    "external program {:?} could not be launched: {error}",
                    spec.program
                )));
            }
        }
    };
    let start = Instant::now();
    let mut status = None;
    let mut timed_out = false;
    let mut cancelled = false;
    loop {
        if cancel.load(Ordering::Relaxed) {
            cancelled = true;
            let _ = child.kill();
            let _ = child.wait();
            break;
        }
        if start.elapsed() >= spec.timeout {
            timed_out = true;
            let _ = child.kill();
            let _ = child.wait();
            break;
        }
        match child.try_wait() {
            Ok(Some(value)) => {
                status = Some(value);
                break;
            }
            Ok(None) => std::thread::sleep(DEFAULT_POLL),
            Err(error) => {
                let _ = child.kill();
                let _ = child.wait();
                return Err(Error::InvalidInput(format!(
                    "external process wait failed: {error}"
                )));
            }
        }
    }
    let run_status = if cancelled {
        RunStatus::Cancelled
    } else if timed_out {
        RunStatus::TimedOut
    } else if status.is_some_and(|value| value.success()) {
        RunStatus::Success
    } else {
        RunStatus::Failed
    };
    let (stdout, stdout_truncated) = read_log(&stdout_path, spec.max_log_bytes)?;
    let (stderr, stderr_truncated) = read_log(&stderr_path, spec.max_log_bytes)?;
    Ok(RunOutcome {
        status: run_status,
        exit_code: status.and_then(|value| value.code()),
        stdout,
        stderr,
        stdout_truncated,
        stderr_truncated,
        scratch: scratch.to_path_buf(),
        duration: start.elapsed(),
    })
}

fn read_log(path: &Path, limit: usize) -> Result<(String, bool), Error> {
    let bytes = fs::read(path).map_err(|error| {
        Error::InvalidInput(format!("external log {}: {error}", path.display()))
    })?;
    let truncated = bytes.len() > limit;
    let bytes = bytes.get(..limit).unwrap_or(&bytes);
    Ok((String::from_utf8_lossy(bytes).into_owned(), truncated))
}

fn create_scratch() -> Result<PathBuf, Error> {
    let base = std::env::temp_dir();
    let stamp = SystemTime::now()
        .duration_since(UNIX_EPOCH)
        .map_or(0, |value| value.as_nanos());
    for attempt in 0..100 {
        let path = base.join(format!(
            "xtalforge-{}-{stamp}-{attempt}",
            std::process::id()
        ));
        match fs::create_dir(&path) {
            Ok(()) => return Ok(path),
            Err(error) if error.kind() == std::io::ErrorKind::AlreadyExists => continue,
            Err(error) => {
                return Err(Error::InvalidInput(format!(
                    "external scratch folder: {error}"
                )));
            }
        }
    }
    Err(Error::InvalidInput(
        "external scratch folder name could not be created".to_owned(),
    ))
}
