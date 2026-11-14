//! The workflow document and its validation rules.
//!
//! The serde types are the wire format; [`Workflow::validate`] is the
//! authoritative rule set and mirrors the published JSON Schema. Every rule
//! reports the JSON path of the offending value so an agent can repair the
//! document without guessing.

use std::path::Path;

use serde::{Deserialize, Serialize};
use serde_json::Value;

use super::WorkflowError;

/// The `workflow_version` this build reads and writes.
pub const WORKFLOW_VERSION: u32 = 1;

/// One declarative workflow.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct Workflow {
    /// Model version; must be [`WORKFLOW_VERSION`].
    pub workflow_version: u32,
    /// Short name used in reports and checkpoints.
    pub name: String,
    /// Optional prose for humans and assistants.
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub description: Option<String>,
    /// Ordered pipeline; the first step is always a `read`.
    pub steps: Vec<Step>,
}

/// What a step does.
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash, Serialize, Deserialize)]
#[serde(rename_all = "snake_case")]
pub enum StepKind {
    /// Load one structure into a fresh document.
    Read,
    /// Apply a mutating verb from the catalogue.
    Transform,
    /// Apply a read-only analysis verb.
    Analyse,
    /// Apply a read-only force-field or numeric verb.
    Compute,
    /// Write the active document to a path.
    Write,
}

impl StepKind {
    /// Stable lower-case name used by JSON and the schema.
    #[must_use]
    pub fn as_str(self) -> &'static str {
        match self {
            Self::Read => "read",
            Self::Transform => "transform",
            Self::Analyse => "analyse",
            Self::Compute => "compute",
            Self::Write => "write",
        }
    }

    /// Whether the kind runs an agent verb.
    #[must_use]
    pub fn is_verb(self) -> bool {
        matches!(self, Self::Transform | Self::Analyse | Self::Compute)
    }
}

/// One step of a workflow.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct Step {
    /// Stable label used in results; defaults to `step-N` (1-based).
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub name: Option<String>,
    /// What the step does.
    pub kind: StepKind,
    /// Catalogue verb for `transform`, `analyse` and `compute`.
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub verb: Option<String>,
    /// Verb arguments, verbatim.
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub arguments: Option<serde_json::Map<String, Value>>,
    /// Path template for `read` and `write`.
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub path: Option<String>,
    /// Format name or extension for `read` and `write`.
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub format: Option<String>,
}

/// One step as validated, with its effective label and path.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct StepSummary {
    /// Effective label: the declared name or `step-N`.
    pub name: String,
    /// What the step does.
    pub kind: StepKind,
    /// Catalogue verb for verb steps.
    pub verb: Option<String>,
    /// Effective path template for `read` and `write`.
    pub path: Option<String>,
    /// Format name or extension for `read` and `write`.
    pub format: Option<String>,
}

/// A validated workflow, ready to describe or run.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct WorkflowSummary {
    /// Model version.
    pub workflow_version: u32,
    /// Workflow name.
    pub name: String,
    /// Optional prose.
    pub description: Option<String>,
    /// Steps in run order.
    pub steps: Vec<StepSummary>,
}

/// One machine-readable problem, used by validation answers and item reports.
#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
pub struct Problem {
    /// Stable code, for example `workflow_invalid` or `io`.
    pub code: String,
    /// Human-readable message.
    pub message: String,
    /// JSON path or filesystem path, when known.
    pub path: Option<String>,
    /// Step label, for runtime failures.
    pub step: Option<String>,
}

impl Workflow {
    /// Parses and validates a JSON document.
    ///
    /// # Errors
    ///
    /// Returns [`WorkflowError::Parse`] for malformed JSON,
    /// [`WorkflowError::Version`] for another model version and
    /// [`WorkflowError::Invalid`] for a rule violation.
    pub fn from_json(text: &str) -> Result<Self, WorkflowError> {
        let workflow: Self = serde_json::from_str(text).map_err(|error| WorkflowError::Parse {
            reason: error.to_string(),
        })?;
        workflow.validate()?;
        Ok(workflow)
    }

    /// Parses and validates a JSON value.
    ///
    /// # Errors
    ///
    /// Returns the errors of [`Self::from_json`].
    pub fn from_value(value: Value) -> Result<Self, WorkflowError> {
        let workflow: Self =
            serde_json::from_value(value).map_err(|error| WorkflowError::Parse {
                reason: error.to_string(),
            })?;
        workflow.validate()?;
        Ok(workflow)
    }

    /// Parses and validates a YAML document.
    ///
    /// Only the documented workflow subset is accepted (block mappings and
    /// sequences, one-line flow collections, quoted and plain scalars,
    /// comments, null/booleans/numbers). Anchors, aliases, tags, directives,
    /// block scalars and multiple documents are rejected. The document is
    /// converted to a JSON value and validated by the same rules as
    /// [`Self::from_json`]. A JSON document is valid YAML and is accepted
    /// through the JSON reader when it spans several lines.
    ///
    /// # Errors
    ///
    /// Returns [`WorkflowError::Parse`] with line and column context for
    /// syntax outside the subset, and the errors of [`Self::from_json`] for a
    /// value of the wrong shape, an unknown field, or a rule violation.
    pub fn from_yaml(text: &str) -> Result<Self, WorkflowError> {
        match super::yaml::parse(text) {
            Ok(value) => Self::from_value(value),
            Err(yaml_error) => {
                if looks_like_json(text) {
                    Self::from_json(text)
                } else {
                    Err(yaml_error)
                }
            }
        }
    }

    /// Parses and validates a document, detecting JSON versus YAML.
    ///
    /// JSON is tried first because it is the canonical format; a document that
    /// is not JSON is read with the YAML subset of [`Self::from_yaml`].
    ///
    /// # Errors
    ///
    /// Returns the error of the parser selected for the document.
    pub fn from_text(text: &str) -> Result<Self, WorkflowError> {
        match Self::from_json(text) {
            Ok(workflow) => Ok(workflow),
            Err(json_error) => Self::from_yaml(text).map_err(|yaml_error| {
                if looks_like_json(text) {
                    json_error
                } else {
                    yaml_error
                }
            }),
        }
    }

    /// Reads and validates a workflow file.
    ///
    /// `.yaml` and `.yml` (in any case) select the YAML reader, which also
    /// accepts JSON; every other extension, including `.json` and no
    /// extension, auto-detects JSON versus YAML like [`Self::from_text`].
    ///
    /// # Errors
    ///
    /// Returns [`WorkflowError::Io`] when the file cannot be read and the
    /// errors of [`Self::from_text`] otherwise.
    pub fn from_path(path: impl AsRef<Path>) -> Result<Self, WorkflowError> {
        let path = path.as_ref();
        let text = std::fs::read_to_string(path).map_err(|error| WorkflowError::Io {
            path: path.display().to_string(),
            reason: error.to_string(),
        })?;
        let extension = path
            .extension()
            .and_then(|extension| extension.to_str())
            .map(str::to_ascii_lowercase);
        match extension.as_deref() {
            Some("yaml" | "yml") => Self::from_yaml(&text),
            _ => Self::from_text(&text),
        }
    }

    /// Renders the workflow as indented JSON.
    ///
    /// # Errors
    ///
    /// Returns [`WorkflowError::Parse`] when a value cannot be serialized.
    pub fn to_json_pretty(&self) -> Result<String, WorkflowError> {
        serde_json::to_string_pretty(self).map_err(|error| WorkflowError::Parse {
            reason: error.to_string(),
        })
    }

    /// Checks every rule and returns the effective step labels and paths.
    ///
    /// # Errors
    ///
    /// Returns [`WorkflowError::Version`] or [`WorkflowError::Invalid`] with
    /// the JSON path of the offending value.
    pub fn validate(&self) -> Result<WorkflowSummary, WorkflowError> {
        super::rules::validate(self)
    }
}

fn looks_like_json(text: &str) -> bool {
    let text = text.strip_prefix('\u{feff}').unwrap_or(text);
    let trimmed = text.trim_start();
    trimmed.starts_with('{') || trimmed.starts_with('[')
}
