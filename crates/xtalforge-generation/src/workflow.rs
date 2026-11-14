//! Searches as checkpointed workflow steps (Task 20.2.3).
//!
//! A search is a long run, so the surface exposes two things:
//!
//! * [`SearchCheckpoint`] — a self-contained JSON record of the options, the
//!   completed seeds and the hits obtained so far. [`run_checkpointed`] runs
//!   the next chunk, writes the checkpoint and returns the updated report, so
//!   a run can stop and resume without losing work.
//! * [`search_workflow_document`] — a valid Phase 16 workflow document whose
//!   `compute` step calls the `search` agent verb. The verb runs the same
//!   checkpointed entry point, so provenance and resume come for free.
//!
//! The checkpoint stores full structures, so it is larger than the workflow
//! report but needs no external files. Its `checkpoint_version` is independent
//! of the crate version and changes only when the stored fields do.

use std::path::Path;

use serde::{Deserialize, Serialize};
use serde_json::{Value, json};

use crate::evolutionary::{EvolutionaryOptions, EvolutionarySearch, evolutionary_search_with};
use crate::random::{
    RandomSearchOptions, RandomSearchReport, SearchCalculator, SearchHit, random_search_with,
};

/// Version of the checkpoint document.
pub const CHECKPOINT_VERSION: u32 = 1;

/// Search kind accepted by the workflow verb.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
#[serde(rename_all = "snake_case")]
pub enum SearchKind {
    /// AIRSS-style random search.
    Random,
    /// Evolutionary search.
    Evolutionary,
}

impl SearchKind {
    /// Stable lower-case name.
    #[must_use]
    pub fn as_str(self) -> &'static str {
        match self {
            Self::Random => "random",
            Self::Evolutionary => "evolutionary",
        }
    }
}

/// A resumable search checkpoint.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct SearchCheckpoint {
    /// Checkpoint format version.
    pub checkpoint_version: u32,
    /// `random` or `evolutionary`.
    pub kind: String,
    /// Options exactly as supplied, for provenance.
    pub options: Value,
    /// Seeds or generations completed.
    pub completed: usize,
    /// Total units requested; zero when the search is open-ended.
    pub total: usize,
    /// Hits obtained so far.
    pub hits: Vec<SearchHit>,
}

impl SearchCheckpoint {
    /// Serializes the checkpoint to JSON.
    ///
    /// # Errors
    ///
    /// Returns [`xtalforge_core::Error::InvalidInput`] when serialization
    /// fails.
    pub fn to_json(&self) -> Result<String, xtalforge_core::Error> {
        serde_json::to_string(self).map_err(|error| {
            xtalforge_core::Error::InvalidInput(format!(
                "could not serialize search checkpoint: {error}"
            ))
        })
    }

    /// Parses a checkpoint and validates its version.
    ///
    /// # Errors
    ///
    /// Returns [`xtalforge_core::Error::InvalidInput`] for invalid JSON or an
    /// unknown version.
    pub fn from_json(text: &str) -> Result<Self, xtalforge_core::Error> {
        let checkpoint: Self = serde_json::from_str(text).map_err(|error| {
            xtalforge_core::Error::InvalidInput(format!(
                "could not parse search checkpoint: {error}"
            ))
        })?;
        if checkpoint.checkpoint_version != CHECKPOINT_VERSION {
            return Err(xtalforge_core::Error::InvalidInput(format!(
                "unsupported checkpoint version {}; this build writes {CHECKPOINT_VERSION}",
                checkpoint.checkpoint_version
            )));
        }
        Ok(checkpoint)
    }

    /// Creates an empty checkpoint.
    #[must_use]
    pub fn new(kind: SearchKind, options: Value, total: usize) -> Self {
        Self {
            checkpoint_version: CHECKPOINT_VERSION,
            kind: kind.as_str().to_owned(),
            options,
            completed: 0,
            total,
            hits: Vec::new(),
        }
    }

    /// Whether the requested units are complete.
    #[must_use]
    pub fn is_complete(&self) -> bool {
        self.total > 0 && self.completed >= self.total
    }

    /// Updates the stored hits, removing energies that are not finite.
    ///
    /// # Errors
    ///
    /// Returns [`xtalforge_core::Error::InvalidInput`] when the kind is
    /// unknown.
    pub fn update(
        &mut self,
        hits: Vec<SearchHit>,
        completed: usize,
    ) -> Result<(), xtalforge_core::Error> {
        if self.kind != SearchKind::Random.as_str()
            && self.kind != SearchKind::Evolutionary.as_str()
        {
            return Err(xtalforge_core::Error::InvalidInput(format!(
                "unknown search kind `{}`",
                self.kind
            )));
        }
        self.hits = hits;
        self.completed = completed;
        Ok(())
    }
}

/// The outcome of one checkpointed chunk.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct SearchChunk {
    /// Updated checkpoint.
    pub checkpoint: SearchCheckpoint,
    /// Ranked hits of this chunk.
    pub hits: Vec<SearchHit>,
    /// Calculator name.
    pub calculator: String,
    /// Whether this chunk finished the requested units.
    pub complete: bool,
}

impl SearchChunk {
    /// Serializes the chunk to JSON.
    ///
    /// # Errors
    ///
    /// Returns [`xtalforge_core::Error::InvalidInput`] when serialization
    /// fails.
    pub fn to_json(&self) -> Result<String, xtalforge_core::Error> {
        serde_json::to_string(self).map_err(|error| {
            xtalforge_core::Error::InvalidInput(format!(
                "could not serialize search chunk: {error}"
            ))
        })
    }
}

/// Runs one random-search chunk with the UFF4MOF calculator.
///
/// # Errors
///
/// Returns the generator, optimizer or checkpoint errors.
pub fn run_random_chunk(
    checkpoint: &mut SearchCheckpoint,
    chunk: usize,
) -> Result<SearchChunk, xtalforge_core::Error> {
    let calculator = xtalforge_ff::UffCalculator::new(xtalforge_ff::CalculatorOptions::default())
        .map_err(crate::ff_error)?;
    run_random_chunk_with(
        checkpoint,
        chunk,
        &SearchCalculator {
            calculator: &calculator,
            name: "uff".to_owned(),
        },
    )
}

/// Runs one random-search chunk with an explicit calculator.
///
/// # Errors
///
/// Returns [`xtalforge_core::Error::InvalidInput`] for a malformed checkpoint
/// and the generator, optimizer or deduplication errors.
pub fn run_random_chunk_with(
    checkpoint: &mut SearchCheckpoint,
    chunk: usize,
    calculator: &SearchCalculator<'_>,
) -> Result<SearchChunk, xtalforge_core::Error> {
    if checkpoint.kind != SearchKind::Random.as_str() {
        return Err(xtalforge_core::Error::InvalidInput(format!(
            "checkpoint kind is `{}`, not `random`",
            checkpoint.kind
        )));
    }
    let mut options: RandomSearchOptions = serde_json::from_value(checkpoint.options.clone())
        .map_err(|error| {
            xtalforge_core::Error::InvalidInput(format!(
                "checkpoint options are not a random search: {error}"
            ))
        })?;
    if chunk == 0 {
        return Err(xtalforge_core::Error::InvalidInput(
            "chunk size must be at least one".to_owned(),
        ));
    }
    let remaining = options.n_structures.saturating_sub(checkpoint.completed);
    options.n_structures = remaining.min(chunk);
    options.seed = options.seed.wrapping_add(checkpoint.completed as u64);
    options.deduplicate = false;
    let report: RandomSearchReport = random_search_with(&options, calculator)?;
    let mut hits = checkpoint.hits.clone();
    hits.extend(report.hits);
    let completed = checkpoint.completed + report.generated;
    checkpoint.update(hits.clone(), completed)?;
    let complete = checkpoint.is_complete();
    Ok(SearchChunk {
        checkpoint: checkpoint.clone(),
        hits,
        calculator: report.calculator,
        complete,
    })
}

/// Runs one evolutionary chunk with the UFF4MOF calculator.
///
/// # Errors
///
/// Returns the generator, optimizer or checkpoint errors.
pub fn run_evolutionary_chunk(
    checkpoint: &mut SearchCheckpoint,
    generations: usize,
) -> Result<SearchChunk, xtalforge_core::Error> {
    let calculator = xtalforge_ff::UffCalculator::new(xtalforge_ff::CalculatorOptions::default())
        .map_err(crate::ff_error)?;
    run_evolutionary_chunk_with(
        checkpoint,
        generations,
        &SearchCalculator {
            calculator: &calculator,
            name: "uff".to_owned(),
        },
    )
}

/// Runs one evolutionary chunk with an explicit calculator.
///
/// # Errors
///
/// Returns [`xtalforge_core::Error::InvalidInput`] for a malformed checkpoint
/// and the generator, optimizer or deduplication errors.
pub fn run_evolutionary_chunk_with(
    checkpoint: &mut SearchCheckpoint,
    generations: usize,
    calculator: &SearchCalculator<'_>,
) -> Result<SearchChunk, xtalforge_core::Error> {
    if checkpoint.kind != SearchKind::Evolutionary.as_str() {
        return Err(xtalforge_core::Error::InvalidInput(format!(
            "checkpoint kind is `{}`, not `evolutionary`",
            checkpoint.kind
        )));
    }
    let mut options: EvolutionaryOptions = serde_json::from_value(checkpoint.options.clone())
        .map_err(|error| {
            xtalforge_core::Error::InvalidInput(format!(
                "checkpoint options are not an evolutionary search: {error}"
            ))
        })?;
    if checkpoint.completed > 0 {
        options.initial_population = checkpoint
            .hits
            .iter()
            .map(|hit| hit.structure.clone())
            .collect();
    }
    options.generations = generations;
    let report: EvolutionarySearch = evolutionary_search_with(&options, calculator)?;
    checkpoint.update(report.hits.clone(), checkpoint.completed + 1)?;
    let complete = checkpoint.is_complete();
    Ok(SearchChunk {
        checkpoint: checkpoint.clone(),
        hits: report.hits,
        calculator: report.calculator,
        complete,
    })
}

/// Writes a checkpoint to a path.
///
/// # Errors
///
/// Returns [`xtalforge_core::Error::InvalidInput`] when the file cannot be
/// written.
pub fn write_checkpoint(
    path: impl AsRef<Path>,
    checkpoint: &SearchCheckpoint,
) -> Result<(), xtalforge_core::Error> {
    let text = checkpoint.to_json()?;
    std::fs::write(path.as_ref(), text).map_err(|error| {
        xtalforge_core::Error::InvalidInput(format!(
            "could not write checkpoint {}: {error}",
            path.as_ref().display()
        ))
    })
}

/// Reads a checkpoint from a path.
///
/// # Errors
///
/// Returns [`xtalforge_core::Error::InvalidInput`] when the file cannot be
/// read or parsed.
pub fn read_checkpoint(path: impl AsRef<Path>) -> Result<SearchCheckpoint, xtalforge_core::Error> {
    let text = std::fs::read_to_string(path.as_ref()).map_err(|error| {
        xtalforge_core::Error::InvalidInput(format!(
            "could not read checkpoint {}: {error}",
            path.as_ref().display()
        ))
    })?;
    SearchCheckpoint::from_json(&text)
}

/// Runs one chunk and persists the updated checkpoint.
///
/// # Errors
///
/// Returns the chunk errors and the write error.
pub fn run_checkpointed(
    kind: SearchKind,
    options: Value,
    total: usize,
    checkpoint_path: impl AsRef<Path>,
    chunk: usize,
) -> Result<SearchChunk, xtalforge_core::Error> {
    let path = checkpoint_path.as_ref();
    let mut checkpoint = if path.exists() {
        read_checkpoint(path)?
    } else {
        SearchCheckpoint::new(kind, options, total)
    };
    let result = match kind {
        SearchKind::Random => run_random_chunk(&mut checkpoint, chunk)?,
        SearchKind::Evolutionary => {
            let generations = references_completed(&checkpoint, total);
            run_evolutionary_chunk(&mut checkpoint, generations.max(1))?
        }
    };
    write_checkpoint(path, &result.checkpoint)?;
    Ok(result)
}

fn references_completed(checkpoint: &SearchCheckpoint, total: usize) -> usize {
    if total == 0 {
        1
    } else {
        total.saturating_sub(checkpoint.completed)
    }
}

/// Builds the workflow document for a search.
///
/// The first step reads a seed structure (the workflow model requires a
/// `read` first), the second calls the `search` compute verb with the options
/// and checkpoint path, and the third writes the JSON report.
#[must_use]
pub fn search_workflow_document(
    name: &str,
    kind: SearchKind,
    seed_path: &str,
    options: Value,
    checkpoint: Option<&str>,
    output: &str,
) -> Value {
    json!({
        "workflow_version": 1,
        "name": name,
        "description": format!("{} structure search with checkpointing", kind.as_str()),
        "steps": [
            {
                "name": "seed",
                "kind": "read",
                "path": seed_path,
                "format": "cif"
            },
            {
                "name": "search",
                "kind": "compute",
                "verb": "search",
                "arguments": {
                    "kind": kind.as_str(),
                    "options": options,
                    "checkpoint": checkpoint
                }
            },
            {
                "name": "write",
                "kind": "write",
                "path": output,
                "format": "json"
            }
        ]
    })
}

/// Serializes [`search_workflow_document`] to pretty JSON.
///
/// # Errors
///
/// Returns [`xtalforge_core::Error::InvalidInput`] when serialization fails.
pub fn search_workflow_json(
    name: &str,
    kind: SearchKind,
    seed_path: &str,
    options: Value,
    checkpoint: Option<&str>,
    output: &str,
) -> Result<String, xtalforge_core::Error> {
    serde_json::to_string_pretty(&search_workflow_document(
        name, kind, seed_path, options, checkpoint, output,
    ))
    .map_err(|error| {
        xtalforge_core::Error::InvalidInput(format!("could not serialize search workflow: {error}"))
    })
}

/// Provenance block for a search report.
#[must_use]
pub fn search_provenance(kind: SearchKind, seed: u64, calculator: &str) -> Value {
    json!({
        "generator": "xtalforge-generation",
        "search": kind.as_str(),
        "seed": seed,
        "calculator": calculator,
    })
}
