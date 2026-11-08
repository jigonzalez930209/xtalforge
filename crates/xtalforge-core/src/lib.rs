//! Core crystal-structure types for xtalforge.
//!
//! This crate holds the element data, lattices, sites, structures and the
//! P1, supercell and slab operations developed in Phase 1. The foundation is
//! kept small so every domain can share a single error type per crate.

#![cfg_attr(
    test,
    allow(
        clippy::unwrap_used,
        clippy::expect_used,
        clippy::panic,
        clippy::indexing_slicing
    )
)]

pub mod data;

/// Errors returned by core operations.
#[derive(Debug, Clone, PartialEq, Eq, thiserror::Error)]
pub enum Error {
    /// An argument does not satisfy the documented preconditions.
    #[error("invalid input: {0}")]
    InvalidInput(String),
    /// A compressed data bundle could not be decompressed.
    #[error("bundle {name}: compressed data is invalid: {reason}")]
    DataCompression {
        /// Name of the bundle as declared in its manifest.
        name: String,
        /// Message from the decompressor.
        reason: String,
    },
    /// A data bundle failed its SHA-256 checksum.
    #[error("bundle {name}: SHA-256 checksum mismatch")]
    DataChecksum {
        /// Name of the bundle as declared in its manifest.
        name: String,
    },
    /// A data bundle declared a digest that is not 64 hexadecimal characters.
    #[error("bundle {name}: invalid SHA-256 digest")]
    DataDigest {
        /// Name of the bundle as declared in its manifest.
        name: String,
    },
    /// JSON decoding of a data bundle failed.
    #[error("bundle {name}: JSON decoding failed: {reason}")]
    DataJson {
        /// Name of the bundle as declared in its manifest.
        name: String,
        /// Message from the JSON decoder.
        reason: String,
    },
}

impl Error {
    /// Stable `snake_case` code for this error, suitable for bindings and JSON.
    ///
    /// The string is part of the public API (see `AGENTS.md`, section 5.2).
    #[must_use]
    pub fn code(&self) -> &'static str {
        match self {
            Self::InvalidInput(_) => "invalid_input",
            Self::DataCompression { .. } => "data_compression",
            Self::DataChecksum { .. } => "data_checksum",
            Self::DataDigest { .. } => "data_digest",
            Self::DataJson { .. } => "data_json",
        }
    }
}

#[cfg(test)]
mod tests {
    use super::Error;

    #[test]
    fn error_code_is_stable() {
        let error = Error::InvalidInput("negative volume".to_owned());
        assert_eq!(error.code(), "invalid_input");
        assert_eq!(error.to_string(), "invalid input: negative volume");
    }
}
