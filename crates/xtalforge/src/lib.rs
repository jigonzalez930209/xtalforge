//! `xtalforge` builds, manipulates, analyses and exports crystal structures.
//!
//! This facade crate re-exports the domain crates under stable module names:
//!
//! - [`core`]: elements, lattices, sites, structures and cell operations.
//! - [`symmetry`]: space groups, Wyckoff positions, standard cells, subgroups,
//!   point groups and Brillouin zones.
//! - [`io`]: file-format readers and writers.
//!
//! Later phases add `ff`, `analysis`, `build`, `mof`, `powder` and `modules`
//! (see the roadmap in `docs/roadmap/`). The same operations are
//! exposed by the CLI, the Python package and the WebAssembly package; the
//! parity matrix in `docs/api/parity.md` tracks their status.
//!
//! # Features
//!
//! | Feature | Default | Effect |
//! |---|---|---|
//! | `parallel` | yes | Enables optional Rayon support in the domain crates. |
//! | `sci-form` | no | Reserved for quantum-chemistry methods (Task 4.10); currently a no-op. |
//! | `ml` | no | Reserved for machine-learned potentials (Task 8.4); currently a no-op. |
//! | `external` | no | Reserved for subprocess engines (Task 4.9); currently a no-op. |
//! | `mof-data` | no | Reserved for the embedded MOF bundle (Task 6.11.3); currently a no-op. |
//!
//! Disable default features with `default-features = false` to build without
//! Rayon; nothing in `parallel` changes results, only execution strategy.
//!
//! # Example
//!
//! ```
//! let error = xtalforge::core::Error::InvalidInput("negative volume".to_owned());
//! assert_eq!(error.code(), "invalid_input");
//! ```

pub use xtalforge_agent as agent;
pub use xtalforge_analysis as analysis;
pub use xtalforge_build as build;
pub use xtalforge_commands as commands;
pub use xtalforge_core as core;
pub use xtalforge_ff as ff;
pub use xtalforge_io as io;
pub use xtalforge_modules as modules;
pub use xtalforge_mof as mof;
pub use xtalforge_symmetry as symmetry;
pub use xtalforge_view as view;

/// Common imports for typical xtalforge programs.
pub mod prelude {
    pub use crate::core::{
        Element, Error, Lattice, Structure, cell_volume, element_info, lattice_info,
        lattice_minimum_image, structure_info, version,
    };
    pub use crate::symmetry::{SymmetryDataset, SymmetryOptions, detect_symmetry};
}
