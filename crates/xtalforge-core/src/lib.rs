//! Core crystal-structure types for xtalforge.
//!
//! This crate holds the element data, lattices, sites, structures and the
//! P1, supercell and slab operations developed in Phase 1. The foundation is
//! kept small so every domain can share a single error type per crate.

pub mod cell;
pub mod data;
pub mod elements;
pub mod lattice;
pub mod limits;
pub mod p1;
pub mod site;
pub mod structure;
pub mod symop;

pub use cell::cell_volume;
pub use elements::{Element, Isotope, element_info};
pub use lattice::{Lattice, lattice_info, lattice_minimum_image};
pub use limits::{Limits, Profile, Verdict, limits_bounds, limits_check};
pub use p1::{P1Atom, P1Cell};
pub use site::Site;
pub use structure::{
    SCHEMA_VERSION, Structure, StructureBuilder, StructureDocument, StructureSummary,
    structure_info,
};
pub use symop::SymOp;

/// Version of the workspace, from the package metadata.
///
/// The value is the same on every surface: it comes from the workspace
/// `Cargo.toml` and is mirrored by the Python and WASM packages at release.
#[must_use]
pub fn version() -> &'static str {
    env!("CARGO_PKG_VERSION")
}

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
    /// A size limit would be exceeded; the verdict explains what would fit.
    #[error("{0}")]
    LimitExceeded(crate::limits::Verdict),
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
            Self::LimitExceeded(_) => "limit_exceeded",
        }
    }

    /// Human-readable message for this error.
    #[must_use]
    pub fn message(&self) -> String {
        self.to_string()
    }
}
