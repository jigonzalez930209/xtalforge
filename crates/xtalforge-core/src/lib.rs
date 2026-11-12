//! Core crystal-structure types for xtalforge.
//!
//! This crate holds the element data, lattices, sites, structures and the
//! P1, supercell and slab operations developed in Phase 1. The foundation is
//! kept small so every domain can share a single error type per crate.

pub mod bonding;
pub mod cell;
pub mod compare;
pub mod data;
pub mod deprecation;
pub mod describe;
pub mod edits;
pub mod elements;
pub mod eos;
pub mod groups;
pub mod hull;
pub mod lattice;
pub mod limits;
pub mod local;
pub mod measure;
pub mod neighbors;
pub mod p1;
pub mod parallel;
pub mod param_presets;
pub mod params;
pub mod polyhedron;
pub mod rings;
pub mod scalars;
pub mod selection;
pub mod site;
pub mod slab;
mod spatial;
pub mod structure;
pub mod supercell;
pub mod surface;
pub mod symop;
pub mod transforms;

pub use bonding::{Bond, BondGraph, BondRules, ExplicitBond, TopologyBond, perceive_bonds};
pub use cell::cell_volume;
pub use compare::{Alignment, interpolate, kabsch_alignment, periodic_rmsd};
pub use deprecation::warn_deprecated;
pub use describe::describe_atom;
pub use edits::{
    RattleOptions, RattleReport, SiteSort, SpeciesReplacement, SupercellChoice, add_sites, deform,
    move_site, optimal_supercell, rattle, remove_sites, replace_site, replace_species,
    scale_to_density, scale_to_volume, sort_sites, strain, tensor_to_voigt, voigt_to_tensor,
    wrap_sites,
};
pub use elements::{Element, Isotope, element_info};
pub use eos::{EosFit, EosKind, EosPoint, eos_energy, fit_eos};
pub use groups::{GroupKind, GroupMatch, find_groups, group_counts};
pub use hull::{ConvexHull, PhaseEntry, PhasePoint, convex_hull};
pub use lattice::{Lattice, lattice_info, lattice_minimum_image};
pub use limits::{Limits, Profile, Verdict, limits_bounds, limits_check};
pub use local::{
    AngleDistribution, CoordinationNeighbor, CoordinationShell, CoordinationStrategy,
    RadialDistribution, RdfOptions, ShapeMeasure, StructureFingerprint, VoronoiCell,
    VoronoiNeighbor, VoronoiReport, bond_angle_distribution, coordination, fingerprint,
    radial_distribution, shape_measures, voronoi,
};
pub use measure::{angle, centroid, distance, interplanar_angle, torsion};
pub use neighbors::{Neighbor, neighbors, neighbors_of};
pub use p1::{P1Atom, P1Cell};
pub use parallel::{ParallelOptions, set_global_threads};
pub use param_presets::{PRESET_NAMES, preset, presets};
pub use params::{ParamSet, ParamSpec, ParamType};
pub use polyhedron::{FULL_SPHERE, HalfSpace, Polyhedron};
pub use rings::{Ring, RingStep, is_aromatic, ring_counts, ring_membership, ring_paths, rings};
pub use scalars::{
    ScalarFn, ScalarInfo, ScalarRegistry, atom_charges, atom_masses, bond_lengths,
    coordination_numbers,
};
pub use selection::{AtomId, Selection, TrackingTable};
pub use site::Site;
pub use slab::{Slab, make_slab, plane_angle, surface_basis};
pub use structure::{
    SCHEMA_VERSION, Structure, StructureBuilder, StructureDocument, StructureSummary,
    StructureWarning, structure_info,
};
pub use surface::{
    AdsorptionKind, AdsorptionSite, SurfaceEnergy, WulffFacet, WulffShape, adsorption_sites,
    place_molecule, wulff,
};
pub use symop::SymOp;
pub use transforms::{Ellipsoid, Plane, displacement_ellipsoid};

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
    /// A symmetry search or space-group lookup failed.
    #[error("symmetry: {0}")]
    Symmetry(String),
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
            Self::Symmetry(_) => "symmetry",
        }
    }

    /// Human-readable message for this error.
    #[must_use]
    pub fn message(&self) -> String {
        self.to_string()
    }

    /// Stable machine-readable category for the JSON error object.
    ///
    /// Categories are `input`, `unsupported`, `numerical`, `limit`, `engine`
    /// and `internal`.
    #[must_use]
    pub fn category(&self) -> &'static str {
        match self {
            Self::InvalidInput(_) => "input",
            Self::DataCompression { .. }
            | Self::DataChecksum { .. }
            | Self::DataDigest { .. }
            | Self::DataJson { .. } => "internal",
            Self::LimitExceeded(_) => "limit",
            Self::Symmetry(_) => "numerical",
        }
    }
}
