//! The resumable checkpoint of a workflow run (Task 16.1.2).
//!
//! A checkpoint is written after every completed item and read on resume. It
//! records the workflow model version and a fingerprint of the document, so a
//! checkpoint from another workflow is refused instead of mixing results.

use std::path::Path;

use serde::{Deserialize, Serialize};

use super::WorkflowError;
use super::model::Workflow;
use super::runner::ItemReport;

/// Version of the checkpoint file layout.
pub const CHECKPOINT_VERSION: u32 = 1;

/// The resumable state of a run.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct WorkflowCheckpoint {
    /// Layout version; see [`CHECKPOINT_VERSION`].
    pub checkpoint_version: u32,
    /// Workflow model version the items were produced with.
    pub workflow_version: u32,
    /// Workflow name.
    pub name: String,
    /// Workflow fingerprint.
    pub fingerprint: String,
    /// Completed and failed items, in input order.
    pub items: Vec<ItemReport>,
}

impl WorkflowCheckpoint {
    /// Reads and parses a checkpoint file.
    ///
    /// # Errors
    ///
    /// Returns [`WorkflowError::Checkpoint`] for an unreadable or malformed
    /// file.
    pub fn read(path: impl AsRef<Path>) -> Result<Self, WorkflowError> {
        let path = path.as_ref();
        let text = std::fs::read_to_string(path).map_err(|error| WorkflowError::Checkpoint {
            path: path.display().to_string(),
            reason: error.to_string(),
        })?;
        serde_json::from_str(&text).map_err(|error| WorkflowError::Checkpoint {
            path: path.display().to_string(),
            reason: error.to_string(),
        })
    }

    /// Writes the checkpoint atomically (temporary file, then rename).
    ///
    /// # Errors
    ///
    /// Returns [`WorkflowError::Io`] or [`WorkflowError::Checkpoint`].
    pub fn write(&self, path: impl AsRef<Path>) -> Result<(), WorkflowError> {
        let path = path.as_ref();
        if let Some(parent) = path
            .parent()
            .filter(|parent| !parent.as_os_str().is_empty())
        {
            std::fs::create_dir_all(parent).map_err(|error| WorkflowError::Io {
                path: parent.display().to_string(),
                reason: error.to_string(),
            })?;
        }
        let text =
            serde_json::to_string_pretty(self).map_err(|error| WorkflowError::Checkpoint {
                path: path.display().to_string(),
                reason: error.to_string(),
            })?;
        let temporary = path.with_extension("tmp");
        std::fs::write(&temporary, text.as_bytes()).map_err(|error| WorkflowError::Io {
            path: temporary.display().to_string(),
            reason: error.to_string(),
        })?;
        std::fs::rename(&temporary, path).map_err(|error| WorkflowError::Io {
            path: path.display().to_string(),
            reason: error.to_string(),
        })
    }
}

/// Stable fingerprint of a workflow document.
///
/// # Errors
///
/// Returns [`WorkflowError::Parse`] when the document cannot be serialized.
pub fn fingerprint(workflow: &Workflow) -> Result<String, WorkflowError> {
    let text = serde_json::to_string(workflow).map_err(|error| WorkflowError::Parse {
        reason: error.to_string(),
    })?;
    let mut hash = 0xcbf2_9ce4_8422_2325_u64;
    for byte in text.as_bytes() {
        hash ^= u64::from(*byte);
        hash = hash.wrapping_mul(0x0000_0100_0000_01b3);
    }
    Ok(format!("fnv1a64:{hash:016x}"))
}
