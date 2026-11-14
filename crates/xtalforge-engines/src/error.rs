//! Errors of the simulation-engine bridge.

use xtalforge_core::Error as CoreError;

/// An input-set, parser or job failure.
#[derive(Debug, Clone, PartialEq, thiserror::Error)]
pub enum Error {
    /// A structure, lattice or numeric value is invalid.
    #[error(transparent)]
    Core(#[from] CoreError),
    /// The requested engine or preset is unknown.
    #[error("unknown engine preset {0}")]
    UnknownPreset(String),
    /// An input-set override is not accepted by the preset.
    #[error("unknown override {0}")]
    UnknownOverride(String),
    /// An output file is malformed or incomplete.
    #[error("{0}")]
    Parse(String),
    /// A required output file or field is missing.
    #[error("missing {0}")]
    Missing(String),
    /// The job or scheduler configuration is invalid.
    #[error("{0}")]
    Job(String),
    /// The external feature or program is unavailable.
    #[error("{0}")]
    Unavailable(String),
}

impl Error {
    /// Stable error code for the CLI, Python and WASM surfaces.
    #[must_use]
    pub fn code(&self) -> &'static str {
        match self {
            Self::Core(error) => error.code(),
            Self::UnknownPreset(_) => "unknown_preset",
            Self::UnknownOverride(_) => "unknown_override",
            Self::Parse(_) => "engine_parse",
            Self::Missing(_) => "engine_missing",
            Self::Job(_) => "engine_job",
            Self::Unavailable(_) => "engine_unavailable",
        }
    }
}
