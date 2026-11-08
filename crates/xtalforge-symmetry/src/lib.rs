//! Crystal symmetry for xtalforge.
//!
//! Space groups, Wyckoff positions, standard and primitive cells, Niggli,
//! Delaunay and Minkowski reductions, subgroups, point groups, Brillouin zones
//! and k-paths live here (Phase 2). `moyo` is an implementation detail: this
//! crate exposes its own types so the backend can be replaced without changing
//! the public API (Task 0.3.5).

pub mod convert;
pub mod dataset;
pub mod spacegroup;

pub use convert::AtomType;
pub use dataset::{SettingChoice, SymmetryDataset, SymmetryOptions, detect_symmetry};
pub use spacegroup::{SpaceGroup, WyckoffPosition, all_settings, crystal_system, lattice_system};
