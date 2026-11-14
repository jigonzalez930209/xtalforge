//! Shared validation and atomic-file helpers for job execution.

use std::ffi::OsString;
use std::path::{Component, Path, PathBuf};

use crate::Error;

/// Accepts a relative path made only of normal components.
pub(crate) fn validate_relative_path(value: &str, what: &str) -> Result<(), Error> {
    let path = Path::new(value);
    if value.trim().is_empty()
        || path.is_absolute()
        || path
            .components()
            .any(|component| !matches!(component, Component::Normal(_)))
        || value.chars().any(char::is_control)
    {
        return Err(Error::Job(format!(
            "{what} {value:?} must be a relative path without `..`"
        )));
    }
    Ok(())
}

/// Accepts a non-empty single-token value (no whitespace or control chars).
pub(crate) fn validate_token(value: &str, what: &str) -> Result<(), Error> {
    if value.trim().is_empty()
        || value
            .chars()
            .any(|character| character.is_whitespace() || character.is_control())
    {
        return Err(Error::Job(format!(
            "{what} {value:?} must be a non-empty single token"
        )));
    }
    Ok(())
}

/// Writes `text` through a sibling temporary file and an atomic rename.
pub(crate) fn write_atomic(path: &Path, text: &str) -> Result<(), Error> {
    if let Some(parent) = path
        .parent()
        .filter(|parent| !parent.as_os_str().is_empty())
    {
        std::fs::create_dir_all(parent).map_err(|error| {
            Error::Job(format!(
                "job directory {} could not be created: {error}",
                parent.display()
            ))
        })?;
    }
    let mut temporary: OsString = path.as_os_str().to_owned();
    temporary.push(".tmp");
    let temporary = PathBuf::from(temporary);
    std::fs::write(&temporary, text.as_bytes()).map_err(|error| {
        Error::Job(format!(
            "job file {} could not be written: {error}",
            temporary.display()
        ))
    })?;
    std::fs::rename(&temporary, path).map_err(|error| {
        Error::Job(format!(
            "job file {} could not be replaced: {error}",
            path.display()
        ))
    })
}
