//! Declarative workflows (Tasks 16.1.1 to 16.1.3).
//!
//! A workflow is a versioned JSON document that describes one pipeline applied
//! to every input structure. A step is a `read`, a `transform`, an `analyse`, a
//! `compute` or a `write`; the middle ones name an agent verb from the shared
//! catalogue, so a workflow can do exactly what an assistant can do and no
//! second list of operations exists to drift.
//!
//! The document is validated against the JSON Schema in [`schema_json`]
//! (published at [`WORKFLOW_SCHEMA_ID`]) and against the live verb catalogue.
//! Bad JSON is a typed [`WorkflowError::Parse`], not a panic.
//!
//! # Example
//!
//! ```
//! use xtalforge_agent::workflow::Workflow;
//!
//! let text = r#"{
//!   "workflow_version": 1,
//!   "name": "inspect",
//!   "steps": [
//!     {"kind": "read", "path": "{{input}}"},
//!     {"kind": "analyse", "verb": "info"},
//!     {"kind": "write", "path": "out/{{stem}}.cif"}
//!   ]
//! }"#;
//! let workflow = Workflow::from_json(text).unwrap();
//! assert_eq!(workflow.validate().unwrap().steps.len(), 3);
//! ```
//!
//! The runner ([`run_workflow`]) expands an input pattern ([`expand_inputs`]),
//! processes one structure per session, isolates item failures, keeps the
//! report in input order and can resume from a JSON checkpoint.
//!
//! Only JSON is supported. YAML was considered (the roadmap names
//! `wf.yaml`); the pure-Rust YAML crates available for this workspace pull
//! dependencies that are not part of the locked set, so the build stays
//! JSON-only. JSON is a subset of YAML, so a JSON document stored with a
//! `.yaml` extension still loads; any other YAML syntax is reported as a parse
//! error that says so.

mod checkpoint;
mod exec;
mod inputs;
mod model;
mod rules;
mod runner;
mod schema;
mod template;
mod verbs;

pub use checkpoint::{CHECKPOINT_VERSION, WorkflowCheckpoint, fingerprint};
pub use inputs::expand_inputs;
pub use model::{
    Problem, Step, StepKind, StepSummary, WORKFLOW_VERSION, Workflow, WorkflowSummary,
};
pub use runner::{ItemReport, ItemStatus, RunOptions, StepReport, WorkflowReport, run_workflow};
pub use schema::{WORKFLOW_SCHEMA_ID, schema_json, schema_text, tool_schema};
pub use template::{TemplateContext, render_template};

pub(crate) use verbs::{run_answer, schema_answer, validate_answer};

/// An error raised while reading, validating or running a workflow.
#[derive(Debug, Clone, PartialEq, thiserror::Error)]
pub enum WorkflowError {
    /// The document is not valid JSON, or does not have the declared shape.
    #[error("workflow JSON: {reason}")]
    Parse {
        /// Parser message with the line and column.
        reason: String,
    },
    /// The document parsed but violates a rule of the model.
    #[error("invalid workflow at {path}: {reason}")]
    Invalid {
        /// JSON path of the offending value, for example `steps[1].verb`.
        path: String,
        /// What is wrong.
        reason: String,
    },
    /// The document declares a `workflow_version` this build does not read.
    #[error("unsupported workflow_version {found}; this build reads version {supported}")]
    Version {
        /// Version found in the document.
        found: u64,
        /// Version this build reads.
        supported: u32,
    },
    /// An input pattern could not be expanded.
    #[error("inputs {pattern}: {reason}")]
    Inputs {
        /// Pattern as the caller wrote it.
        pattern: String,
        /// What is wrong.
        reason: String,
    },
    /// A checkpoint is unreadable, inconsistent or belongs to another workflow.
    #[error("checkpoint {path}: {reason}")]
    Checkpoint {
        /// Checkpoint path.
        path: String,
        /// What is wrong.
        reason: String,
    },
    /// A file could not be read or written.
    #[error("io {path}: {reason}")]
    Io {
        /// Affected path.
        path: String,
        /// Operating-system message.
        reason: String,
    },
}

impl WorkflowError {
    /// Stable `snake_case` code for bindings and JSON errors.
    #[must_use]
    pub fn code(&self) -> &'static str {
        match self {
            Self::Parse { .. } => "workflow_parse",
            Self::Invalid { .. } => "workflow_invalid",
            Self::Version { .. } => "workflow_version",
            Self::Inputs { .. } => "inputs",
            Self::Checkpoint { .. } => "checkpoint",
            Self::Io { .. } => "io",
        }
    }

    /// Stable machine-readable category.
    #[must_use]
    pub fn category(&self) -> &'static str {
        match self {
            Self::Parse { .. }
            | Self::Invalid { .. }
            | Self::Version { .. }
            | Self::Inputs { .. } => "input",
            Self::Checkpoint { .. } | Self::Io { .. } => "io",
        }
    }

    /// JSON path or filesystem path the error refers to, when there is one.
    #[must_use]
    pub fn path(&self) -> Option<&str> {
        match self {
            Self::Invalid { path, .. } | Self::Checkpoint { path, .. } | Self::Io { path, .. } => {
                Some(path)
            }
            Self::Inputs { pattern, .. } => Some(pattern),
            Self::Parse { .. } | Self::Version { .. } => None,
        }
    }

    /// Renders the error as a workflow [`Problem`].
    #[must_use]
    pub fn problem(&self) -> Problem {
        Problem {
            code: self.code().to_owned(),
            message: self.to_string(),
            path: self.path().map(str::to_owned),
            step: None,
        }
    }
}
