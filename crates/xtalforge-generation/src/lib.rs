//! Structure generation and global search for xtalforge (Phase 20).
//!
//! The crate turns a composition and a space group into a crystal, places
//! rigid molecules on Wyckoff positions, builds layer and rod structures,
//! searches the energy landscape with the native force field, and enumerates
//! hypothetical MOFs. Every random decision comes from an explicit seed
//! (AGENTS.md rule 6); the same options and seed give the same structures.
//!
//! * [`generate_crystal`] — Wyckoff-constrained random crystals (Task 20.1.1).
//! * [`generate_molecular_crystal`] — rigid molecules on general and special
//!   positions (Task 20.1.2).
//! * [`generate_layer_crystal`] and [`generate_rod_crystal`] — low-dimensional
//!   structures (Task 20.1.3).
//! * [`random_search`] and [`evolutionary_search`] — AIRSS-style relaxation
//!   and evolutionary global search (Tasks 20.2.1 and 20.2.2), with
//!   checkpointed workflow documents (Task 20.2.3).
//! * [`prototype_library`], [`identify_prototype`] and
//!   [`substitution_candidates`] — the AFLOW prototype subset (Task 20.3).
//! * [`generate_mof_dataset`], [`global_orientation_fit`] and
//!   [`functionalization_sweep`] — hypothetical MOFs at scale (Task 20.4).

pub mod evolutionary;
mod evolutionary_check;
pub mod functionalize;
pub mod geometry;
pub mod lowdim;
pub mod mofgen;
pub mod molecular;
pub mod prototype;
pub mod random;
pub mod rng;
pub mod rod;
pub mod substitution;
pub mod workflow;
pub mod wyckoff;

/// Converts a force-field error into the crate's core error type.
pub(crate) fn ff_error(error: xtalforge_ff::Error) -> xtalforge_core::Error {
    xtalforge_core::Error::InvalidInput(format!("force field: {}", error.message()))
}

pub use evolutionary::{
    EvolutionaryOptions, EvolutionarySearch, GenerationReport, evolutionary_search,
    evolutionary_search_with,
};
pub use functionalize::{FunctionalizationRecord, FunctionalizationSweep, functionalization_sweep};
pub use lowdim::{
    LAYER_VACUUM, LayerCrystalOptions, LowDimCrystal, LowDimOptions, generate_layer_crystal,
    layer_groups, layer_groups_catalogue,
};
pub use mofgen::{
    MofDatasetOptions, MofDatasetReport, MofGenerationRecord, OrientationFit, generate_mof_dataset,
    global_orientation_fit,
};
pub use molecular::{MolecularCrystal, MolecularCrystalOptions, generate_molecular_crystal};
pub use prototype::{
    Prototype, PrototypeIdentification, PrototypeLattice, PrototypeOptions, PrototypeSite,
    build_prototype, identify_prototype, prototype_library, prototype_names,
};
pub use random::{
    RandomSearchOptions, RandomSearchReport, SearchCalculator, SearchHit, SearchMinimizer,
    random_search, random_search_with,
};
pub use rod::{RodCrystalOptions, RodGroup, generate_rod_crystal, rod_groups};
pub use substitution::{
    IonicRadius, SubstitutionCandidate, SubstitutionOptions, SubstitutionReport,
    apply_substitution, group_by_element, ionic_radii, substituted_structures,
    substitution_candidates, substitution_report,
};
pub use workflow::{
    CHECKPOINT_VERSION, SearchCheckpoint, SearchChunk, SearchKind, read_checkpoint,
    run_checkpointed, run_evolutionary_chunk, run_evolutionary_chunk_with, run_random_chunk,
    run_random_chunk_with, search_provenance, search_workflow_document, search_workflow_json,
    write_checkpoint,
};
pub use wyckoff::{
    CrystalOptions, GeneratedCrystal, GenerationStatistics, crystal_json, generate_crystal,
    generate_crystal_range, generation_statistics,
};
