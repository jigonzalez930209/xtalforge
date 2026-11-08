//! Analysis for xtalforge: PXRD, distance grids, porosity and nets (Phase 5).
//!
//! The first two modules implement the calculated powder pattern (Task 5.1)
//! and the periodic surface-distance field with the porosity numbers read off
//! it (Tasks 5.2 and 5.3). Every public operation returns a typed [`Error`]
//! and never panics on user input.

pub mod fit;
pub mod grid;
pub mod indexing;
pub mod isosurface;
pub mod lebail;
pub mod merge;
pub mod nets;
pub mod overlay;
pub mod pareto;
pub mod peaks;
pub mod profiles;
pub mod pxrd;
pub mod rcsr;
pub mod rietveld;
pub mod zeo;

pub use fit::{FitOptions, FitResult, least_squares};
pub use grid::{DistanceGrid, PoreSizeBin, default_radii, distance_grid};
pub use indexing::{CubicIndex, IndexedPeak, index_cubic};
pub use isosurface::{TriangleMesh, isosurface};
pub use lebail::{LebailOptions, LebailResult, le_bail, pawley};
pub use merge::merge_nets;
pub use nets::{NetEdge, NetIndex, NetInvariants, QuotientGraph, interpenetrate, quotient_graph};
pub use overlay::{Overlay, OverlaySet, OverlayStyle, overlay_net};
pub use pareto::pareto_front;
pub use peaks::{Peak, PeakOptions, find_peaks};
pub use profiles::{Background, Profile, ProfileKind};
pub use pxrd::{Pattern, PatternPoint, PxrdOptions, Radiation, Reflection, calculate_pxrd};
pub use rcsr::{RCSR_LICENSE_NOTE, RcsrEntry, RcsrIndex};
pub use rietveld::{RietveldOptions, RietveldPhase, fit_phase};
pub use xtalforge_core::Error;
#[cfg(not(target_arch = "wasm32"))]
pub use zeo::run_zeo;
pub use zeo::{ZeoOptions, ZeoResult, parse_zeo_channels, parse_zeo_volume, write_zeo_inputs};
