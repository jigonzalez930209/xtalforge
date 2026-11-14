//! Process launch and runtime probes, compiled only with the `external`
//! feature.
//!
//! Without the feature the same functions exist and return
//! [`Error::Unavailable`], so callers do not need conditional compilation.

use std::path::{Path, PathBuf};
use std::time::Duration;

use serde::{Deserialize, Serialize};

use crate::Error;

#[cfg(feature = "external")]
const POLL_INTERVAL: Duration = Duration::from_millis(10);
#[cfg(feature = "external")]
const SPAWN_ATTEMPTS: usize = 10;

/// One engine invocation.
#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
pub struct EngineCommand {
    /// Executable name, absolute path or path relative to the job directory.
    pub program: String,
    /// Command-line arguments.
    pub args: Vec<String>,
}

impl EngineCommand {
    /// Builds an invocation.
    #[must_use]
    pub fn new(program: impl Into<String>, args: Vec<String>) -> Self {
        Self {
            program: program.into(),
            args,
        }
    }

    /// Rejects an empty program name.
    ///
    /// # Errors
    ///
    /// Returns [`Error::Job`] when the program name is empty.
    pub fn validate(&self) -> Result<(), Error> {
        if self.program.trim().is_empty() {
            return Err(Error::Job(
                "an engine command needs a program name".to_owned(),
            ));
        }
        if self.program.chars().any(char::is_control) {
            return Err(Error::Job(
                "an engine program name must not contain control characters".to_owned(),
            ));
        }
        Ok(())
    }
}

/// Whether this build can launch external processes and reach a scheduler.
///
/// Always `false` for the default pure-Rust build and for WASM.
#[must_use]
pub const fn process_execution_available() -> bool {
    cfg!(feature = "external")
}

/// Captured result of one engine process.
pub(crate) struct ProcessOutcome {
    /// Exit status zero.
    pub success: bool,
    /// The wall-clock limit expired and the process was killed.
    pub timed_out: bool,
    /// Exit code when the platform reported one.
    pub exit_code: Option<i32>,
    /// Captured standard output, capped at the requested size.
    pub stdout: String,
    /// Captured standard error, capped at the requested size.
    pub stderr: String,
}

/// Captured result of one scheduler command.
pub(crate) struct CommandOutcome {
    /// Exit status zero.
    pub success: bool,
    /// Exit code when the platform reported one.
    pub exit_code: Option<i32>,
    /// Captured standard output.
    pub stdout: String,
    /// Captured standard error.
    pub stderr: String,
}

/// Runs one engine in `cwd`, capturing logs and enforcing a wall-clock limit.
#[cfg(feature = "external")]
pub(crate) fn run_process(
    command: &EngineCommand,
    cwd: &Path,
    timeout: Duration,
    max_log_bytes: usize,
) -> Result<ProcessOutcome, Error> {
    use std::fs::File;
    use std::process::{Command, Stdio};

    let stdout_path = cwd.join(".xtalforge-stdout.log");
    let stderr_path = cwd.join(".xtalforge-stderr.log");
    let mut child = None;
    for attempt in 0..SPAWN_ATTEMPTS {
        let stdout = File::create(&stdout_path).map_err(|error| {
            Error::Job(format!(
                "engine stdout log {}: {error}",
                stdout_path.display()
            ))
        })?;
        let stderr = File::create(&stderr_path).map_err(|error| {
            Error::Job(format!(
                "engine stderr log {}: {error}",
                stderr_path.display()
            ))
        })?;
        match Command::new(&command.program)
            .args(&command.args)
            .current_dir(cwd)
            .stdin(Stdio::null())
            .stdout(Stdio::from(stdout))
            .stderr(Stdio::from(stderr))
            .spawn()
        {
            Ok(value) => {
                child = Some(value);
                break;
            }
            // ETXTBSY: another thread forked while the test script was being
            // written. It clears as soon as that child exits.
            Err(error) if error.raw_os_error() == Some(26) && attempt + 1 < SPAWN_ATTEMPTS => {
                std::thread::sleep(POLL_INTERVAL);
            }
            Err(error) => {
                return Err(Error::Unavailable(format!(
                    "engine program {:?} could not be launched: {error}",
                    command.program
                )));
            }
        }
    }
    let mut child = child.ok_or_else(|| {
        Error::Unavailable(format!(
            "engine program {:?} could not be launched",
            command.program
        ))
    })?;
    let start = std::time::Instant::now();
    let mut timed_out = false;
    let status = loop {
        if start.elapsed() >= timeout {
            timed_out = true;
            let _ = child.kill();
            let _ = child.wait();
            break None;
        }
        match child.try_wait() {
            Ok(Some(status)) => break Some(status),
            Ok(None) => std::thread::sleep(POLL_INTERVAL),
            Err(error) => {
                let _ = child.kill();
                let _ = child.wait();
                return Err(Error::Job(format!("engine process wait failed: {error}")));
            }
        }
    };
    Ok(ProcessOutcome {
        success: status.is_some_and(|status| status.success()),
        timed_out,
        exit_code: status.and_then(|status| status.code()),
        stdout: read_capped(&stdout_path, max_log_bytes)?,
        stderr: read_capped(&stderr_path, max_log_bytes)?,
    })
}

/// Runs one engine in `cwd`, capturing logs and enforcing a wall-clock limit.
#[cfg(not(feature = "external"))]
pub(crate) fn run_process(
    _command: &EngineCommand,
    _cwd: &Path,
    _timeout: Duration,
    _max_log_bytes: usize,
) -> Result<ProcessOutcome, Error> {
    Err(Error::Unavailable(
        "engine process execution requires the `external` feature".to_owned(),
    ))
}

/// Runs one scheduler command, capturing its output.
#[cfg(feature = "external")]
pub(crate) fn run_command(program: &str, args: &[String]) -> Result<CommandOutcome, Error> {
    use std::process::{Command, Stdio};

    let mut last_error = None;
    for attempt in 0..SPAWN_ATTEMPTS {
        match Command::new(program)
            .args(args)
            .stdin(Stdio::null())
            .output()
        {
            Ok(output) => {
                return Ok(CommandOutcome {
                    success: output.status.success(),
                    exit_code: output.status.code(),
                    stdout: String::from_utf8_lossy(&output.stdout).into_owned(),
                    stderr: String::from_utf8_lossy(&output.stderr).into_owned(),
                });
            }
            Err(error) if error.raw_os_error() == Some(26) && attempt + 1 < SPAWN_ATTEMPTS => {
                std::thread::sleep(POLL_INTERVAL);
                last_error = Some(error);
            }
            Err(error) => {
                return Err(Error::Unavailable(format!(
                    "scheduler program {program:?} could not be launched: {error}"
                )));
            }
        }
    }
    let reason = last_error.map_or_else(
        || "the spawn retry budget was exhausted".to_owned(),
        |error| error.to_string(),
    );
    Err(Error::Unavailable(format!(
        "scheduler program {program:?} could not be launched: {reason}"
    )))
}

/// Runs one scheduler command, capturing its output.
#[cfg(not(feature = "external"))]
pub(crate) fn run_command(_program: &str, _args: &[String]) -> Result<CommandOutcome, Error> {
    Err(Error::Unavailable(
        "scheduler commands require the `external` feature".to_owned(),
    ))
}

/// Resolves a program name against `PATH`, or checks an explicit path.
///
/// # Errors
///
/// Returns a human-readable reason when the program is missing or not
/// executable.
pub(crate) fn probe_program(program: &str) -> Result<PathBuf, String> {
    if program.trim().is_empty() {
        return Err("the command name is empty".to_owned());
    }
    let path = Path::new(program);
    let explicit = path.is_absolute() || path.components().count() > 1;
    if explicit {
        return if is_executable(path) {
            Ok(path.to_path_buf())
        } else {
            Err(format!("{program:?} is not an executable file"))
        };
    }
    let Some(paths) = std::env::var_os("PATH") else {
        return Err("PATH is not set".to_owned());
    };
    for directory in std::env::split_paths(&paths) {
        let candidate = directory.join(program);
        if is_executable(&candidate) {
            return Ok(candidate);
        }
    }
    Err(format!("{program:?} was not found on PATH"))
}

fn is_executable(path: &Path) -> bool {
    let Ok(metadata) = std::fs::metadata(path) else {
        return false;
    };
    if !metadata.is_file() {
        return false;
    }
    #[cfg(unix)]
    {
        use std::os::unix::fs::PermissionsExt;
        metadata.permissions().mode() & 0o111 != 0
    }
    #[cfg(not(unix))]
    {
        true
    }
}

#[cfg(feature = "external")]
fn read_capped(path: &Path, limit: usize) -> Result<String, Error> {
    let bytes = std::fs::read(path).map_err(|error| {
        Error::Job(format!(
            "engine log {} could not be read: {error}",
            path.display()
        ))
    })?;
    let bytes = bytes.get(..limit).unwrap_or(&bytes);
    Ok(String::from_utf8_lossy(bytes).into_owned())
}
