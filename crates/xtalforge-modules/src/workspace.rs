//! On-disk workspace of runs (Task 3.12).

use std::path::{Path, PathBuf};

use serde::{Deserialize, Serialize};
use xtalforge_core::Error;

/// One run record.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct RunRecord {
    /// Stable run id.
    pub id: String,
    /// Module or command that produced it.
    pub module: String,
    /// Input summary.
    pub input: String,
    /// Status.
    pub status: RunStatus,
    /// Output file names relative to the run directory.
    pub outputs: Vec<String>,
    /// Free-form parameters and results.
    pub data: serde_json::Value,
    /// Unix timestamp in seconds.
    pub created: u64,
}

/// Run status.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
pub enum RunStatus {
    /// Queued.
    Pending,
    /// Running.
    Running,
    /// Finished successfully.
    Done,
    /// Failed with a message in `data`.
    Failed,
}

/// A workspace directory with one subdirectory per run.
#[derive(Debug, Clone)]
pub struct Workspace {
    root: PathBuf,
}

impl Workspace {
    /// Opens or creates a workspace.
    ///
    /// # Errors
    ///
    /// Returns [`Error::InvalidInput`] when the directory cannot be created.
    #[cfg(not(target_arch = "wasm32"))]
    pub fn open(root: impl Into<PathBuf>) -> Result<Self, Error> {
        let root = root.into();
        std::fs::create_dir_all(&root).map_err(|error| {
            Error::InvalidInput(format!(
                "could not create workspace {}: {error}",
                root.display()
            ))
        })?;
        Ok(Self { root })
    }

    /// Workspace root.
    #[must_use]
    pub fn root(&self) -> &Path {
        &self.root
    }

    /// Directory of one run.
    #[must_use]
    pub fn run_dir(&self, id: &str) -> PathBuf {
        self.root.join("runs").join(id)
    }

    /// Writes a run record and creates its directory.
    ///
    /// # Errors
    ///
    /// Returns [`Error::InvalidInput`] when the directory or record cannot be
    /// written.
    #[cfg(not(target_arch = "wasm32"))]
    pub fn write_run(&self, record: &RunRecord) -> Result<PathBuf, Error> {
        let directory = self.run_dir(&record.id);
        std::fs::create_dir_all(&directory).map_err(|error| {
            Error::InvalidInput(format!(
                "could not create run {}: {error}",
                directory.display()
            ))
        })?;
        let path = directory.join("run.json");
        let text = serde_json::to_string_pretty(record)
            .map_err(|error| Error::InvalidInput(format!("could not serialize run: {error}")))?;
        std::fs::write(&path, text).map_err(|error| {
            Error::InvalidInput(format!("could not write {}: {error}", path.display()))
        })?;
        Ok(path)
    }

    /// Reads a run record.
    ///
    /// # Errors
    ///
    /// Returns [`Error::InvalidInput`] when the record is missing or invalid.
    #[cfg(not(target_arch = "wasm32"))]
    pub fn read_run(&self, id: &str) -> Result<RunRecord, Error> {
        let path = self.run_dir(id).join("run.json");
        let text = std::fs::read_to_string(&path).map_err(|error| {
            Error::InvalidInput(format!("could not read {}: {error}", path.display()))
        })?;
        serde_json::from_str(&text).map_err(|error| {
            Error::InvalidInput(format!("invalid run record {}: {error}", path.display()))
        })
    }

    /// Lists run ids in directory order.
    ///
    /// # Errors
    ///
    /// Returns [`Error::InvalidInput`] when the runs directory cannot be read
    /// (a missing directory is an empty list).
    #[cfg(not(target_arch = "wasm32"))]
    pub fn runs(&self) -> Result<Vec<String>, Error> {
        let directory = self.root.join("runs");
        let entries = match std::fs::read_dir(&directory) {
            Ok(entries) => entries,
            Err(error) if error.kind() == std::io::ErrorKind::NotFound => return Ok(Vec::new()),
            Err(error) => {
                return Err(Error::InvalidInput(format!(
                    "could not read {}: {error}",
                    directory.display()
                )));
            }
        };
        let mut ids = Vec::new();
        for entry in entries.flatten() {
            if entry.path().is_dir()
                && let Some(name) = entry.file_name().to_str()
            {
                ids.push(name.to_owned());
            }
        }
        ids.sort();
        Ok(ids)
    }
}
