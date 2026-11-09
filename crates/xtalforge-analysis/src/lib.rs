//! Analysis for xtalforge: PXRD, distance grids, porosity and nets (Phase 5).
//!
//! The first two modules implement the calculated powder pattern (Task 5.1)
//! and the periodic surface-distance field with the porosity numbers read off
//! it (Tasks 5.2 and 5.3). Every public operation returns a typed [`Error`]
//! and never panics on user input.

pub mod automatic;
pub mod cells;
pub mod fit;
pub mod grid;
pub mod indexing;
pub mod instrument;
pub mod isosurface;
pub mod lebail;
pub mod merge;
pub mod nets;
pub mod overlay;
pub mod pareto;
pub mod pawley;
pub mod peaks;
pub mod profiles;
pub mod pxrd;
pub mod rcsr;
pub mod refine;
pub mod rietveld;
pub mod rietx;
pub mod zeo;

pub use automatic::{AutomaticOptions, AutomaticResult, LeadingCell, run_automatic};
pub use cells::{
    CELL_PARAMETER_NAMES, TrigonalSetting, extinction_symbol, free_cell_parameters,
    free_cell_parameters_with_setting,
};
pub use fit::{FitOptions, FitResult, least_squares};
pub use grid::{DistanceGrid, PoreSizeBin, default_radii, distance_grid};
pub use indexing::{
    CellMetric, CubicIndex, FiguresOfMerit, IndexedCell, IndexedPeak, IndexingOptions,
    IndexingResult, LatticeSystem, figures_of_merit, index_cubic, index_pattern, match_peaks,
    peak_q, predict_reflections, two_theta_of_q,
};
pub use instrument::{Geometry, Instrument, doublet, doublet_checked};
pub use isosurface::{TriangleMesh, isosurface};
pub use lebail::{LebailOptions, LebailResult, le_bail, pawley};
pub use merge::merge_nets;
pub use nets::{NetEdge, NetIndex, NetInvariants, QuotientGraph, interpenetrate, quotient_graph};
pub use overlay::{Overlay, OverlaySet, OverlayStyle, overlay_net};
pub use pareto::pareto_front;
pub use pawley::{PawleyOptions, PawleyReflection, PawleyResult, pawley_refine};
pub use peaks::{Peak, PeakOptions, background_envelope, find_peaks, fit_peaks};
pub use profiles::{Background, Profile, ProfileKind};
pub use pxrd::{Pattern, PatternPoint, PxrdOptions, Radiation, Reflection, calculate_pxrd};
pub use rcsr::{RCSR_LICENSE_NOTE, RcsrEntry, RcsrIndex};
pub use refine::{RefineOptions, RefinedPhase, combined_objective, free_directions, refine_phase};
pub use rietveld::{RietveldOptions, RietveldPhase, fit_phase};
pub use rietx::{
    ACCEPTANCE_MAX_POSITION_ERROR, ACCEPTANCE_MAX_RWP, AcceptanceMetrics, RIETX_LICENSE_NOTE,
    RietxBridge, check_acceptance, rietx_inventory,
};
pub use xtalforge_core::Error;
#[cfg(not(target_arch = "wasm32"))]
pub use zeo::run_zeo;
pub use zeo::{ZeoOptions, ZeoResult, parse_zeo_channels, parse_zeo_volume, write_zeo_inputs};
