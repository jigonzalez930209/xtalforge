//! Typed resumable state for a compute job (Task 18.4.1).

use std::path::{Path, PathBuf};

use serde::{Deserialize, Serialize};

use crate::Error;
use crate::result::CalculationResult;

use super::compute::JobRequest;
use super::util::{validate_relative_path, write_atomic};
use super::{FailurePattern, JobState};

/// Version of the checkpoint layout.
pub const JOB_CHECKPOINT_VERSION: u32 = 1;

/// The resumable state of one compute job.
///
/// A checkpoint is written before the process starts (`running`), when the
/// build cannot launch processes (`prepared`), and after every terminal
/// outcome. Resuming a `completed` checkpoint with a matching fingerprint
/// returns the typed result without launching the engine; a failed or
/// interrupted checkpoint is restarted. A dry run never modifies a checkpoint.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct JobCheckpoint {
    /// Layout version; see [`JOB_CHECKPOINT_VERSION`].
    pub checkpoint_version: u32,
    /// Job name.
    pub name: String,
    /// Fingerprint of the request that produced this state.
    pub fingerprint: String,
    /// Last recorded state.
    pub state: JobState,
    /// Number of engine launches attempted.
    pub attempts: u32,
    /// Exit code when the platform reported one.
    pub exit_code: Option<i32>,
    /// Captured standard output, capped by the run options.
    pub stdout: String,
    /// Captured standard error, capped by the run options.
    pub stderr: String,
    /// Parsed result, when parsing succeeded.
    pub result: Option<CalculationResult>,
    /// Detected failure pattern, when one was recognised.
    pub failure: Option<FailurePattern>,
    /// Error message for a failure that has no engine pattern (parse error,
    /// unknown parser or launch failure).
    pub error: Option<String>,
}

impl JobCheckpoint {
    /// Reads and parses a checkpoint.
    ///
    /// # Errors
    ///
    /// Returns [`Error::Job`] for an unreadable or malformed file.
    pub fn read(path: impl AsRef<Path>) -> Result<Self, Error> {
        let path = path.as_ref();
        let text = std::fs::read_to_string(path).map_err(|error| {
            Error::Job(format!(
                "job checkpoint {} could not be read: {error}",
                path.display()
            ))
        })?;
        serde_json::from_str(&text).map_err(|error| {
            Error::Job(format!(
                "job checkpoint {} is malformed: {error}",
                path.display()
            ))
        })
    }

    /// Reads a checkpoint, returning `None` when the file does not exist.
    ///
    /// # Errors
    ///
    /// Returns [`Error::Job`] for an unreadable or malformed file.
    pub fn read_if_exists(path: impl AsRef<Path>) -> Result<Option<Self>, Error> {
        let path = path.as_ref();
        match std::fs::metadata(path) {
            Ok(_) => Self::read(path).map(Some),
            Err(error) if error.kind() == std::io::ErrorKind::NotFound => Ok(None),
            Err(error) => Err(Error::Job(format!(
                "job checkpoint {} could not be inspected: {error}",
                path.display()
            ))),
        }
    }

    /// Writes the checkpoint atomically.
    ///
    /// # Errors
    ///
    /// Returns [`Error::Job`] when the directory or file cannot be written.
    pub fn write(&self, path: impl AsRef<Path>) -> Result<(), Error> {
        let text = serde_json::to_string_pretty(self)
            .map_err(|error| Error::Job(format!("job checkpoint could not be encoded: {error}")))?;
        write_atomic(path.as_ref(), &text)
    }

    /// Fresh state for a job that has not run yet.
    pub(crate) fn new(name: &str, fingerprint: &str, state: JobState, attempts: u32) -> Self {
        Self {
            checkpoint_version: JOB_CHECKPOINT_VERSION,
            name: name.to_owned(),
            fingerprint: fingerprint.to_owned(),
            state,
            attempts,
            exit_code: None,
            stdout: String::new(),
            stderr: String::new(),
            result: None,
            failure: None,
            error: None,
        }
    }
}

/// Deterministic checkpoint path for a job in `work_dir`.
///
/// # Errors
///
/// Returns [`Error::Job`] for an unsafe job name.
pub fn checkpoint_path(work_dir: impl AsRef<Path>, name: &str) -> Result<PathBuf, Error> {
    validate_relative_path(name, "job name")?;
    Ok(work_dir.as_ref().join(format!("{name}.checkpoint.json")))
}

/// Stable fingerprint of a compute request.
///
/// The fingerprint covers the job name, every input file, the engine command,
/// the parser name and the output file. A checkpoint from another request is
/// refused instead of mixing results.
///
/// # Errors
///
/// Returns [`Error::Job`] when the request cannot be serialized.
pub fn request_fingerprint(request: &JobRequest) -> Result<String, Error> {
    let text = serde_json::to_string(request)
        .map_err(|error| Error::Job(format!("job request could not be fingerprinted: {error}")))?;
    let mut hash = 0xcbf2_9ce4_8422_2325_u64;
    for byte in text.as_bytes() {
        hash ^= u64::from(*byte);
        hash = hash.wrapping_mul(0x0000_0100_0000_01b3);
    }
    Ok(format!("fnv1a64:{hash:016x}"))
}
