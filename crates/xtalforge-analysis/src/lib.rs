//! Analysis for xtalforge: PXRD, distance grids, porosity and nets (Phase 5).
//!
//! The first two modules implement the calculated powder pattern (Task 5.1)
//! and the periodic surface-distance field with the porosity numbers read off
//! it (Tasks 5.2 and 5.3). Every public operation returns a typed [`Error`]
//! and never panics on user input.

pub mod grid;
pub mod isosurface;
pub mod pxrd;

pub use grid::{DistanceGrid, PoreSizeBin, default_radii, distance_grid};
pub use isosurface::{TriangleMesh, isosurface};
pub use pxrd::{Pattern, PatternPoint, PxrdOptions, Radiation, Reflection, calculate_pxrd};
pub use xtalforge_core::Error;
