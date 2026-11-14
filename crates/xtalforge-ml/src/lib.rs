//! Machine-learning representations, datasets and inference (Phase 21).
//!
//! The crate is the preprocessing and inference layer for crystal ML:
//!
//! * [`crystal_graph`] — periodic graphs with image offsets and Cartesian
//!   edge vectors (Task 21.1.1).
//! * [`soap_power_spectrum`] — the exact DScribe-compatible SOAP power
//!   spectrum (Task 21.1.2); [`acsf()`] and [`mbtr()`] add the other exact
//!   descriptors (Task 21.1.3).
//! * [`mof_features()`] — revised autocorrelations and pore geometry
//!   (Task 21.1.4).
//! * [`split_dataset`] — leakage-safe random, composition, cluster and
//!   prototype splits (Task 21.2.2).
//! * [`import_dataset`] — extxyz, OPTIMADE/MP dumps and CoRE MOF/QMOF.
//! * [`PropertyModel`], [`Ensemble`] and the optional `OnnxModel` — native inference with
//!   uncertainty (Tasks 21.3).
//! * [`select_batch`] and [`active_learning_loop`] — active learning
//!   (Task 21.4).
//!
//! The `onnx` feature adds tract-based ONNX inference and the `parquet`
//! feature adds Arrow/Parquet containers; both are pure Rust and optional.

pub mod acsf;
pub mod active;
pub mod dataset;
pub mod ensemble;
pub mod graph;
pub mod import;
pub mod import_cif;
pub mod mbtr;
pub mod mof_features;
#[cfg(feature = "onnx")]
pub mod onnx;
pub mod property;
mod rng;
pub mod select;
pub mod soap;
pub mod split;

pub use acsf::{AcsfDescriptor, AcsfOptions, acsf};
pub use active::{
    ActiveLearningOptions, ActiveLearningReport, ActiveLearningRound, active_learning_loop,
    active_learning_workflow,
};
pub use dataset::{FeatureTable, feature_table};
#[cfg(feature = "parquet")]
pub use dataset::{read_parquet, write_parquet};
pub use ensemble::{Ensemble, EnsemblePrediction, Prediction, uncertainty};
pub use graph::{
    CrystalGraph, GraphArrays, GraphOptions, canonical_p1, crystal_graph, crystal_graph_arrays,
    graph_from_positions,
};
pub use import::{
    ImportFormat, ImportOptions, ImportReport, import_dataset, import_extxyz, import_json_records,
};
pub use import_cif::{import_cif_directory, structure_from_arrays};
pub use mbtr::{MbtrDescriptor, MbtrOptions, mbtr};
pub use mof_features::{
    MofFeatures, MofFeaturesOptions, mof_features, pore_geometry_vector, revised_autocorrelations,
};
#[cfg(feature = "onnx")]
pub use onnx::{OnnxCalculator, OnnxModel, OnnxOptions};
pub use property::{
    Predictor, PropertyModel, WorkflowPrediction, append_predictions, predict_dataset,
};
pub use select::{
    SelectionOptions, SelectionStrategy, batch_uncertainty, farthest_point, select_batch,
};
pub use soap::{SoapDescriptor, SoapOptions, soap_feature_count, soap_power_spectrum};
pub use split::{
    SplitAssignment, SplitOptions, SplitStrategy, SplitSummary, split_dataset, split_summary,
};
