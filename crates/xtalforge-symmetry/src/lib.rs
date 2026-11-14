//! Crystal symmetry for xtalforge.
//!
//! Space groups, Wyckoff positions, standard and primitive cells, Niggli,
//! Delaunay and Minkowski reductions, subgroups, point groups, Brillouin zones
//! and k-paths live here (Phase 2). `moyo` is an implementation detail: this
//! crate exposes its own types so the backend can be replaced without changing
//! the public API (Task 0.3.5).

pub mod brillouin;
pub mod convert;
pub mod dataset;
pub mod diagnostics;
pub mod kpath;
pub mod layer;
pub mod magnetic;
pub mod matcher;
mod matcher_canonical;
mod matcher_lattice;
pub mod molecular;
pub mod reduction;
pub mod spacegroup;
pub mod standardize;
pub mod subgroups;
pub mod workflow;
mod wyckoff_origin;

pub use brillouin::{BrillouinZone, brillouin_zone, brillouin_zone_for_structure};
pub use convert::AtomType;
pub use dataset::{SettingChoice, SymmetryDataset, SymmetryOptions, detect_symmetry};
pub use diagnostics::{
    AtomDeviation, DEFAULT_LADDER, SymmetryDiagnostics, ToleranceRange, ToleranceSample, diagnose,
    tolerance_ladder, tolerance_samples,
};
pub use kpath::{KPath, band_path, band_path_for_structure, classify_bravais};
pub use layer::{LayerDataset, detect_layer_symmetry};
pub use magnetic::{
    MagneticDataset, MagneticOperation, MagneticSymmetrization, MagneticSymmetryOptions,
    detect_magnetic_symmetry, detect_magnetic_symmetry_with_options, format_magmom, magmom_values,
    symmetrize_moments, symmetrize_structure_moments,
};
pub use matcher::{
    MatcherOptions, SpeciesComparator, StructureMatch, classify_structures, group_structures,
    interpolate_matched, match_structures,
};
pub use molecular::{MolecularPointGroup, molecular_point_group};
pub use reduction::{
    ReducedStructure, Reduction, lattices_equivalent, reduce_lattice, reduce_structure,
};
pub use spacegroup::{
    SpaceGroup, SpaceGroupFilter, WyckoffPosition, all_settings, catalogue, crystal_system,
    lattice_system, laue_class,
};
pub use standardize::{
    StandardizedStructure, change_setting, idealize, primitive_cell, standardize,
};
pub use subgroups::{
    Subgroup, SubgroupKind, identify_space_group, klassengleiche_subgroups, primitive_operations,
    translationengleiche_subgroups,
};
pub use workflow::{
    SetMode, SetReport, SetSummary, SymmetryResult, WyckoffAssignment, asymmetrize,
    detect_and_symmetrize, free_parameters, lower_to_subgroup, set_space_group, to_p1,
    wyckoff_assignments,
};
