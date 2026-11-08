//! Errors returned by the MOF crate.

/// Errors returned by topology loading, block loading and framework assembly.
#[derive(Debug, Clone, PartialEq, Eq, thiserror::Error)]
pub enum Error {
    /// An argument does not satisfy the documented preconditions.
    #[error("invalid input: {0}")]
    InvalidInput(String),
    /// The embedded topology or building-block data is missing or corrupt.
    #[error("data: {0}")]
    Data(String),
    /// A topology name is not present in the curated table.
    #[error("unknown topology {0}")]
    UnknownTopology(String),
    /// A building-block name is not present in the curated table.
    #[error("unknown building block {0}")]
    UnknownBlock(String),
    /// A JSON request could not be decoded or encoded.
    #[error("JSON: {0}")]
    Json(String),
    /// An error from the core crate.
    #[error("core: {0}")]
    Core(#[from] xtalforge_core::Error),
}

impl Error {
    /// Stable `snake_case` code for this error, suitable for bindings and JSON.
    ///
    /// The string is part of the public API (see `AGENTS.md`, section 5.2).
    #[must_use]
    pub fn code(&self) -> &'static str {
        match self {
            Self::InvalidInput(_) => "invalid_input",
            Self::Data(_) => "data",
            Self::UnknownTopology(_) => "unknown_topology",
            Self::UnknownBlock(_) => "unknown_block",
            Self::Json(_) => "json",
            Self::Core(error) => error.code(),
        }
    }

    /// Human-readable message for this error.
    #[must_use]
    pub fn message(&self) -> String {
        self.to_string()
    }
}
