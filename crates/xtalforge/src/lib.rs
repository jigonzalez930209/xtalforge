//! `xtalforge` builds, manipulates, analyses and exports crystal structures.
//!
//! This facade crate re-exports the domain crates under stable module names:
//!
//! - [`core`]: elements, lattices, sites, structures and cell operations.
//! - [`symmetry`]: space groups, Wyckoff positions, standard cells, subgroups,
//!   point groups and Brillouin zones.
//! - [`io`]: file-format readers and writers.
//!
//! Later phases add `ff`, `analysis`, `build`, `mof`, `powder`, `modules` and
//! `agent` (see the roadmap in `docs/roadmap/`). The same operations are
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

#![cfg_attr(
    test,
    allow(
        clippy::unwrap_used,
        clippy::expect_used,
        clippy::panic,
        clippy::indexing_slicing
    )
)]

pub use xtalforge_core as core;
pub use xtalforge_io as io;
pub use xtalforge_symmetry as symmetry;

/// Common imports for typical xtalforge programs.
pub mod prelude {
    pub use crate::core::Error;
}
