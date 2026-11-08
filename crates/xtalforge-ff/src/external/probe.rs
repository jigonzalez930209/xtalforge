//! Engine probes: PATH discovery, executable check and a short version run.

use super::{EngineKind, EngineStatus};

/// Probes an engine on `PATH`.
///
/// On WASM the result is always [`EngineStatus::Missing`] with a browser
/// reason; no process is ever launched.
#[cfg(not(target_arch = "wasm32"))]
#[must_use]
pub fn probe(kind: EngineKind) -> EngineStatus {
    probe_binary(kind.binary())
}

/// Probes an engine on `PATH`; WASM always reports it as missing.
#[cfg(target_arch = "wasm32")]
#[must_use]
pub fn probe(kind: EngineKind) -> EngineStatus {
    browser_missing(kind.binary())
}

/// Probes an executable by name on `PATH`.
///
/// This is the shared implementation of [`probe`] and is useful for a binary
/// the caller configured under a non-standard name.  On WASM it always
/// reports a missing engine.
#[cfg(not(target_arch = "wasm32"))]
#[must_use]
pub fn probe_binary(binary: &str) -> EngineStatus {
    if binary.trim().is_empty() {
        return EngineStatus::Missing {
            reason: "no binary name was given".to_owned(),
        };
    }
    let Some(path) = std::env::var_os("PATH") else {
        return EngineStatus::Missing {
            reason: "PATH is not set".to_owned(),
        };
    };
    let mut found = None;
    for directory in std::env::split_paths(&path) {
        let candidate = directory.join(binary);
        if is_executable(&candidate) {
            found = Some(candidate);
            break;
        }
        #[cfg(windows)]
        {
            let candidate = directory.join(format!("{binary}.exe"));
            if is_executable(&candidate) {
                found = Some(candidate);
                break;
            }
        }
    }
    let Some(path) = found else {
        return EngineStatus::Missing {
            reason: format!("{binary} was not found on PATH"),
        };
    };
    match std::process::Command::new(&path)
        .arg("--version")
        .stdin(std::process::Stdio::null())
        .output()
    {
        Ok(output) => {
            let mut text = String::from_utf8_lossy(&output.stdout).into_owned();
            if text.trim().is_empty() {
                text = String::from_utf8_lossy(&output.stderr).into_owned();
            }
            let version = text
                .lines()
                .map(str::trim)
                .find(|line| !line.is_empty())
                .unwrap_or("unknown")
                .chars()
                .take(160)
                .collect::<String>();
            EngineStatus::Available {
                path: path.to_string_lossy().into_owned(),
                version,
            }
        }
        Err(error) => EngineStatus::Missing {
            reason: format!("{} could not be launched: {error}", path.display()),
        },
    }
}

/// Probes an executable by name; WASM always reports it as missing.
#[cfg(target_arch = "wasm32")]
#[must_use]
pub fn probe_binary(binary: &str) -> EngineStatus {
    browser_missing(binary)
}

#[cfg(target_arch = "wasm32")]
fn browser_missing(binary: &str) -> EngineStatus {
    EngineStatus::Missing {
        reason: format!("{binary} is not available in the browser"),
    }
}

#[cfg(unix)]
fn is_executable(path: &std::path::Path) -> bool {
    use std::os::unix::fs::PermissionsExt;
    std::fs::metadata(path)
        .map(|metadata| metadata.is_file() && metadata.permissions().mode() & 0o111 != 0)
        .unwrap_or(false)
}

#[cfg(all(not(unix), not(target_arch = "wasm32")))]
fn is_executable(path: &std::path::Path) -> bool {
    std::fs::metadata(path)
        .map(|metadata| metadata.is_file())
        .unwrap_or(false)
}
