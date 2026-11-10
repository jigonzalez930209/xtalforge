//! Agent layer for xtalforge: session, verbs, inspection, diagnostics,
//! change tracking and pictures (Phase 8).
//!
//! The layer is what an AI assistant drives. Every public operation is a
//! named verb in [`verbs::catalogue`], takes JSON arguments and returns a
//! typed, printable [`answers::Answer`]. The CLI, the Python package, the
//! WebAssembly package and the MCP server read the same catalogue, so there
//! is one list of capabilities and no second list to drift.
//!
//! Mutations never bypass [`xtalforge_commands::Document`]: the session
//! applies them through commands with [`xtalforge_commands::Origin::Agent`],
//! so a user can review or undo the assistant's work as a unit
//! ([`changes`]).
//!
//! # Example
//!
//! ```
//! use xtalforge_agent::{Session, answers::Answer};
//!
//! let mut session = Session::new();
//! let answer = session.run("capabilities", &serde_json::json!({})).unwrap();
//! assert!(matches!(answer, Answer::Capabilities { .. }));
//! ```
//!
//! Ported from Crystal-Builder `xtal/agent/` (MIT), split by verb group; the
//! reference is recorded in `docs/algorithms/agent.md`.

pub mod answers;
pub mod batch;
pub mod capabilities;
pub mod changes;
pub mod context;
pub mod diagnostics;
pub mod inspect;
pub mod mcp;
mod ops;
mod ops_analysis;
mod ops_ff;
mod ops_picture;
pub mod png;
pub mod render;
mod run;
pub mod session;
pub mod skill;
pub mod verbs;

pub use answers::Answer;
pub use batch::run_verbs;
pub use capabilities::{Capabilities, capability_report};
pub use changes::{ChangeBurst, changes_for};
pub use diagnostics::{Diagnostic, DiagnosticReport, Severity, diagnose};
pub use inspect::inspect_structure;
pub use mcp::{McpServer, serve_stdio};
pub use render::{RenderOptions, render_png, render_svg};
pub use session::{Session, SessionOptions};
pub use skill::{SkillDocument, default_directory};
pub use verbs::{ParameterKind, VerbInfo, catalogue, lookup};

/// Errors returned by the agent layer.
#[derive(Debug, Clone, PartialEq, thiserror::Error)]
pub enum Error {
    /// The verb name is not in the catalogue.
    #[error("unknown verb {0}")]
    UnknownVerb(String),
    /// A required argument is missing or has the wrong JSON type.
    #[error("invalid argument {name}: {reason}")]
    Argument {
        /// Argument name.
        name: String,
        /// What is wrong.
        reason: String,
    },
    /// The document name is not loaded in the session.
    #[error("unknown document {0}")]
    UnknownDocument(String),
    /// A verb needs a document and none is active.
    #[error("no active document; call `load` first")]
    NoDocument,
    /// File input or output failed.
    #[error("io: {0}")]
    Io(String),
    /// An error from the core crate.
    #[error("core: {0}")]
    Core(#[from] xtalforge_core::Error),
    /// An error from the analysis crate.
    #[error("analysis: {0}")]
    Analysis(String),
    /// An error from the force-field crate.
    #[error("force field: {0}")]
    ForceField(#[from] xtalforge_ff::Error),
    /// A JSON value could not be encoded or decoded.
    #[error("json: {0}")]
    Json(String),
}

impl Error {
    /// Stable `snake_case` code for this error, suitable for bindings and JSON.
    #[must_use]
    pub fn code(&self) -> &'static str {
        match self {
            Self::UnknownVerb(_) => "unknown_verb",
            Self::Argument { .. } => "invalid_argument",
            Self::UnknownDocument(_) => "unknown_document",
            Self::NoDocument => "no_document",
            Self::Io(_) => "io",
            Self::Core(error) => error.code(),
            Self::Analysis(_) => "analysis",
            Self::ForceField(error) => error.code(),
            Self::Json(_) => "json",
        }
    }
}
