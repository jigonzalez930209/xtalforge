//! Analysis for xtalforge: PXRD, distance grids, porosity and nets (Phase 5).
//!
//! The first modules implement the calculated powder pattern (Task 5.1) and
//! the periodic surface-distance field with the porosity numbers read off it
//! (Tasks 5.2 and 5.3). Phase 14 adds the full diffraction model: the
//! instrument profile with Caglioti and TCH widths (Task 14.1.1), anisotropic
//! displacement parameters in every structure factor (Task 14.1.2), refinable
//! mixing in the Pawley, Le Bail and Rietveld engines (Task 14.1.3),
//! single-crystal structure factors with HKL export and electron diffraction
//! (Task 14.1.4), and whole-range peak refinement, persistent parameter sets,
//! Pareto sweeps, instrument presets and cancellation (Task 14.1.5).
//!
//! Phase 16 adds the JSON Lines dataset container (Task 16.3.1), dataset
//! deduplication with a fingerprint pre-filter (Task 16.3.2) and the
//! composition, RDF, Coulomb-matrix and SOAP-like descriptors (Task 16.3.3).
//!
//! Every public operation returns a typed [`Error`] and never panics on user
//! input.
//!
//! # Open notes (Task 14.1.5)
//!
//! * A "powder demo" is a documentation and playground item, not a crate
//!   entry point, so it is tracked by Task 14.6.1.
//! * "Indexing engine agreement" (the cubic and general indexers accepting the
//!   same cell on a shared pattern) is covered by the existing
//!   [`index_pattern`] / [`index_cubic`] cross-checks; a dedicated agreement
//!   verdict struct is left to Phase 22 because it would add a public result
//!   type that the frozen 1.x surface cannot deprecate first.

pub mod automatic;
pub mod bvs;
pub mod cells;
pub mod contraction;
mod contraction_metals;
mod contraction_simplify;
pub mod dataset;
pub mod dedup;
pub mod defects;
pub mod descriptors;
pub mod diffraction;
pub mod ewald;
pub mod extinction;
pub mod fit;
pub mod grid;
pub mod indexing;
pub mod instrument;
pub mod interfaces;
pub mod interpenetration;
pub mod isosurface;
pub mod lebail;
pub mod merge;
pub mod nets;
pub mod netsearch;
mod netsearch_query;
pub mod orderings;
pub mod overlay;
pub mod parameters;
pub mod pareto;
pub mod pawley;
pub mod peaks;
pub mod profile_pattern;
pub mod profiles;
pub mod progress;
pub mod pxrd;
pub mod rcsr;
pub mod refine;
pub mod rietveld;
pub mod rietx;
pub mod scattering;
pub mod sqs;
pub mod surfaces;
pub mod topology;
pub mod vertex_symbols;
pub mod zeo;

pub use automatic::{AutomaticOptions, AutomaticResult, LeadingCell, run_automatic};
pub use bvs::{
    BondValenceReport, BondValenceSite, BvsOptions, OxidationReport, OxidationSite,
    bond_valence_sums, oxidation_states,
};
pub use cells::{
    CELL_PARAMETER_NAMES, TrigonalSetting, extinction_symbol, free_cell_parameters,
    free_cell_parameters_with_setting,
};
pub use contraction::{
    ContractedEdge, ContractedNet, ContractedNode, ContractionOptions, NetNodeRule, NodeKind,
    contract as contract_structure,
};
pub use dataset::{DATASET_VERSION, Dataset, DatasetMetadata, DatasetRecord, PropertyValue};
pub use dedup::{DedupGroup, DedupOptions, DedupReport, DedupResult, deduplicate};
pub use defects::{DefectKind, PointDefect, interstitials, substitutions, vacancies};
pub use descriptors::{
    CompositionFeatures, CoulombMatrix, CoulombMatrixOptions, MAX_SOAP_FEATURES, MAX_SOAP_L,
    MAX_SOAP_NEIGHBORS, MAX_SOAP_RADIAL, RdfDescriptor, SoapLikeDescriptor, SoapLikeOptions,
    composition_features, coulomb_matrix, rdf_descriptor, soap_like_descriptor,
};
#[cfg(not(target_arch = "wasm32"))]
pub use diffraction::write_hkl;
pub use diffraction::{
    ElectronDiffractionOptions, ElectronPattern, ElectronReflection, StructureFactor,
    StructureFactorOptions, electron_diffraction, hkl_text, structure_factors,
};
pub use ewald::{EwaldOptions, ewald_energy};
pub use extinction::{
    ExtinctionClass, MAX_CANDIDATES as MAX_SPACE_GROUP_CANDIDATES, SpaceGroupMatch,
    extinction_classes, is_systematically_absent, space_group_candidates,
};
pub use fit::{FitOptions, FitResult, least_squares};
pub use grid::{
    DistanceGrid, Percolation, PoreSizeBin, PorosityOptions, PorosityReport, default_radii,
    distance_grid, porosity_report,
};
pub use indexing::{
    CellMetric, CubicIndex, FiguresOfMerit, IndexedCell, IndexedPeak, IndexingAgreement,
    IndexingOptions, IndexingResult, LatticeSystem, figures_of_merit, index_cubic, index_pattern,
    indexing_agreement, indexing_agreement_with_result, match_peaks, peak_q, predict_reflections,
    two_theta_of_q,
};
pub use instrument::{
    Caglioti, FcjAsymmetry, Geometry, INSTRUMENT_PRESET_COUNT, INSTRUMENT_PRESET_NAMES, Instrument,
    InstrumentPreset, InstrumentProfile, LineShape, Monochromator, RadiationKind, RadiationSource,
    Tch, doublet, doublet_checked, electron_wavelength, instrument_presets, preset_by_name,
    tch_preset_by_name,
};
pub use interfaces::{
    Heterostructure, InterfaceMatch, InterfaceOptions, lattice_matches, stack_slabs,
};
pub use interpenetration::{
    BlatovClass, InterpenetrationReport, classify as classify_interpenetration,
    classify_with_symmetry as classify_interpenetration_with_symmetry, cycle_voltages,
    saturation_index,
};
pub use isosurface::{TriangleMesh, isosurface};
pub use lebail::{
    LebailOptions, LebailProfileOptions, LebailProfileResult, LebailResult, le_bail,
    le_bail_profile, pawley,
};
pub use merge::merge_nets;
pub use nets::{NetEdge, NetIndex, NetInvariants, QuotientGraph, interpenetrate, quotient_graph};
pub use netsearch::{NetFacts, RcsrFactsIndex, plane_group_number, space_group_number};
pub use netsearch_query::NetQuery;
pub use orderings::{Ordering, OrderingOptions, orderings};
pub use overlay::{Overlay, OverlaySet, OverlayStyle, overlay_net};
pub use parameters::{PARAMETER_SET_VERSION, ParameterSet};
pub use pareto::{ParetoSample, ParetoSweep, knee_index, pareto_front, pareto_sweep};
pub use pawley::{
    PawleyOptions, PawleyProfileOptions, PawleyProfileResult, PawleyReflection, PawleyResult,
    pawley_refine, pawley_refine_profile,
};
pub use peaks::{
    Peak, PeakOptions, WholeRangeFit, WholeRangeOptions, background_envelope, find_peaks,
    fit_peaks, refine_whole_range, refine_whole_range_with_control,
};
pub use profile_pattern::{PxrdProfileOptions, calculate_pxrd_profile};
pub use profiles::{Background, Profile, ProfileKind};
pub use progress::{Cancellation, Progress, ProgressCallback, RunControl};
pub use pxrd::{Pattern, PatternPoint, PxrdOptions, Radiation, Reflection, calculate_pxrd};
pub use rcsr::{RCSR_LICENSE_NOTE, RcsrEntry, RcsrIndex};
pub use refine::{RefineOptions, RefinedPhase, combined_objective, free_directions, refine_phase};
pub use rietveld::{
    RietveldOptions, RietveldPhase, RietveldProfileOptions, RietveldProfileResult, fit_phase,
    fit_phase_profile,
};
pub use rietx::{
    ACCEPTANCE_MAX_POSITION_ERROR, ACCEPTANCE_MAX_RWP, AcceptanceMetrics, RIETX_LICENSE_NOTE,
    RietxBridge, check_acceptance, rietx_inventory,
};
pub use scattering::{electron_scattering_factor, neutron_scattering_length, xray_form_factor};
pub use sqs::{SqsOptions, SqsResult, SqsShell, special_quasirandom};
pub use surfaces::{
    MillerIndex, Termination, TerminationReport, minimal_slab, slab_terminations,
    slab_terminations_with_options, unique_miller_indices,
};
pub use topology::{DEPTH, MAX_RING, NetFingerprint, fingerprint};
pub use vertex_symbols::{VertexSymbol, vertex_symbol_text, vertex_symbols};
pub use xtalforge_core::Error;
#[cfg(not(target_arch = "wasm32"))]
pub use zeo::run_zeo;
pub use zeo::{ZeoOptions, ZeoResult, parse_zeo_channels, parse_zeo_volume, write_zeo_inputs};
