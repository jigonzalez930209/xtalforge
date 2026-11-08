//! Atomic file writes (Task 3.2).
//!
//! A writer targets a temporary file in the destination directory, syncs it
//! and renames it over the target. On any error the previous file is intact
//! and the temporary file is removed.

#[cfg(not(target_arch = "wasm32"))]
use std::io::Write;
use std::path::Path;
#[cfg(not(target_arch = "wasm32"))]
use std::path::PathBuf;

use xtalforge_core::Error;

/// Writes bytes atomically to `path`.
///
/// # Errors
///
/// Returns [`Error::InvalidInput`] when the directory, temporary file, sync or
/// rename fails.
#[cfg(not(target_arch = "wasm32"))]
pub fn write_bytes(path: &Path, bytes: &[u8]) -> Result<(), Error> {
    let temporary = temporary_path(path);
    let result = (|| {
        let mut file = std::fs::File::create(&temporary).map_err(|error| {
            Error::InvalidInput(format!("could not create {}: {error}", temporary.display()))
        })?;
        file.write_all(bytes).map_err(|error| {
            Error::InvalidInput(format!("could not write {}: {error}", temporary.display()))
        })?;
        file.sync_all().map_err(|error| {
            Error::InvalidInput(format!("could not sync {}: {error}", temporary.display()))
        })?;
        drop(file);
        std::fs::rename(&temporary, path).map_err(|error| {
            Error::InvalidInput(format!("could not replace {}: {error}", path.display()))
        })
    })();
    if result.is_err() {
        let _ = std::fs::remove_file(&temporary);
    }
    result
}

/// Writes text atomically to `path`.
///
/// # Errors
///
/// Returns the errors of [`write_bytes`].
#[cfg(not(target_arch = "wasm32"))]
pub fn write_text(path: &Path, text: &str) -> Result<(), Error> {
    write_bytes(path, text.as_bytes())
}

/// WASM has no filesystem; callers write the returned string or bytes.
///
/// # Errors
///
/// Always returns [`Error::InvalidInput`] on `wasm32`.
#[cfg(target_arch = "wasm32")]
pub fn write_bytes(_path: &Path, _bytes: &[u8]) -> Result<(), Error> {
    Err(Error::InvalidInput(
        "the browser has no filesystem; write the returned bytes instead".to_owned(),
    ))
}

/// WASM has no filesystem; callers write the returned string or bytes.
///
/// # Errors
///
/// Always returns [`Error::InvalidInput`] on `wasm32`.
#[cfg(target_arch = "wasm32")]
pub fn write_text(_path: &Path, _text: &str) -> Result<(), Error> {
    Err(Error::InvalidInput(
        "the browser has no filesystem; write the returned text instead".to_owned(),
    ))
}

#[cfg(not(target_arch = "wasm32"))]
fn temporary_path(path: &Path) -> PathBuf {
    let name = path
        .file_name()
        .and_then(|name| name.to_str())
        .unwrap_or("structure");
    let stamp = std::time::SystemTime::now()
        .duration_since(std::time::UNIX_EPOCH)
        .map_or(0, |duration| duration.as_nanos());
    path.with_file_name(format!(".{name}.{}.{}.tmp", std::process::id(), stamp))
}
