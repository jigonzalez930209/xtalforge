//! Simulation-engine bridge (Phase 18).
//!
//! The crate turns a structure into reproducible inputs for the major
//! simulation codes, reads their outputs back into typed results, analyses
//! electronic structure and submits jobs. The external programs stay optional:
//! writing and parsing text is pure Rust, and only the scheduler/process paths
//! are behind the `external` feature.

pub mod electronic;
pub mod error;
pub mod input_set;
pub mod job;
pub mod parse;
pub mod result;

pub use electronic::*;
pub use error::Error;
pub use input_set::*;
pub use job::*;
pub use parse::*;
pub use result::{
    BaderEntry, BandData, BandProjection, CalculationResult, Convergence, DdecEntry, DosChannel,
    DosData, ElectronicStep, ProjectionData, Snapshot, Timing, Warning, bader_charges, bader_table,
    ddec_charges, ddec_table,
};
