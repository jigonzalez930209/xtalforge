#![allow(clippy::indexing_slicing)] // bounded indices from validated loops and fixed-size vectors

//! Active-learning loops (Task 21.4.2).
//!
//! [`active_learning_loop`] runs generate → predict → select → compute →
//! append rounds over an existing dataset. The caller supplies a featurizer,
//! a predictor, a retraining closure and a compute closure; the loop records
//! every round, the selected identifiers, the predicted uncertainty and the
//! computed labels.
//!
//! [`active_learning_workflow`] emits the same pipeline as a Phase 16 workflow
//! document. The document names the `search`, `ml_predict`, `ml_select` and
//! `ml_append` verbs; running it needs the corresponding agent verbs and a
//! filesystem for the dataset.

use serde::{Deserialize, Serialize};
use serde_json::Value;
use xtalforge_analysis::dataset::{Dataset, PropertyValue};
use xtalforge_core::{Error, Structure};

use crate::property::Predictor;
use crate::select::{SelectionOptions, select_batch};

/// Options for [`active_learning_loop`].
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct ActiveLearningOptions {
    /// Number of rounds.
    pub rounds: usize,
    /// Selection options.
    pub selection: SelectionOptions,
    /// Property name written back to the dataset.
    pub target: String,
    /// Stop when the batch has no candidate.
    pub stop_when_empty: bool,
}

impl Default for ActiveLearningOptions {
    fn default() -> Self {
        Self {
            rounds: 4,
            selection: SelectionOptions::default(),
            target: "computed".to_owned(),
            stop_when_empty: true,
        }
    }
}

/// One active-learning round.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct ActiveLearningRound {
    /// Round number (1-based).
    pub round: usize,
    /// Number of records before the round.
    pub n_before: usize,
    /// Selected record identifiers.
    pub selected: Vec<String>,
    /// Predictive uncertainty of each selected record.
    pub uncertainty: Vec<f64>,
    /// Labels computed for the selected records.
    pub labels: Vec<f64>,
    /// Number of records after the round.
    pub n_after: usize,
}

/// Result of an active-learning run.
#[derive(Debug, Clone, PartialEq)]
pub struct ActiveLearningReport {
    /// Rounds in order.
    pub rounds: Vec<ActiveLearningRound>,
    /// Final dataset.
    pub dataset: Dataset,
    /// Final model name.
    pub model: String,
}

impl ActiveLearningReport {
    /// Serializes the report without the dataset structures.
    ///
    /// # Errors
    ///
    /// Returns [`Error::InvalidInput`] when serialization fails.
    pub fn to_json(&self) -> Result<String, Error> {
        let rounds = serde_json::to_value(&self.rounds)
            .map_err(|error| Error::InvalidInput(format!("active learning: {error}")))?;
        serde_json::to_string(&serde_json::json!({
            "rounds": rounds,
            "model": self.model,
            "n_records": self.dataset.len(),
        }))
        .map_err(|error| Error::InvalidInput(format!("active learning: {error}")))
    }
}

/// Runs an active-learning loop.
///
/// `featurize` maps a structure to a descriptor, `uncertainty` returns the
/// predictive uncertainty for a descriptor, `compute` returns the label for a
/// structure, and `retrain` may return an updated predictor after each round.
/// If `retrain` returns `None` the current predictor stays.
///
/// # Errors
///
/// Returns the featurizer, predictor, compute or retrain errors.
pub fn active_learning_loop<F, U, C, R>(
    dataset: &Dataset,
    mut model: Box<dyn Predictor>,
    mut featurize: F,
    mut uncertainty: U,
    mut compute: C,
    mut retrain: R,
    options: &ActiveLearningOptions,
) -> Result<ActiveLearningReport, Error>
where
    F: FnMut(&Structure) -> Result<Vec<f64>, Error>,
    U: FnMut(&dyn Predictor, &[f64]) -> Result<f64, Error>,
    C: FnMut(&Structure) -> Result<f64, Error>,
    R: FnMut(&Dataset) -> Result<Option<Box<dyn Predictor>>, Error>,
{
    if options.rounds == 0 {
        return Err(Error::InvalidInput(
            "active learning needs at least one round".to_owned(),
        ));
    }
    let mut current = clone_dataset(dataset)?;
    let mut rounds = Vec::new();
    for round in 1..=options.rounds {
        let n_before = current.len();
        let mut features = Vec::with_capacity(n_before);
        let mut uncertainties = Vec::with_capacity(n_before);
        for record in current.iter() {
            let descriptor = featurize(&record.structure)?;
            let value = uncertainty(model.as_ref(), &descriptor)?;
            features.push(descriptor);
            uncertainties.push(value);
        }
        let selected = select_batch(&features, &uncertainties, &options.selection)?;
        if selected.is_empty() && options.stop_when_empty {
            break;
        }
        let mut ids = Vec::with_capacity(selected.len());
        let mut selected_uncertainty = Vec::with_capacity(selected.len());
        let mut labels = Vec::with_capacity(selected.len());
        let mut additions = Vec::with_capacity(selected.len());
        for index in &selected {
            let Some(record) = current.get(*index) else {
                continue;
            };
            let label = compute(&record.structure)?;
            ids.push(record.id.clone());
            selected_uncertainty.push(uncertainties[*index]);
            labels.push(label);
            let mut clone = record.clone();
            clone.id = format!("{}-al{round}-{}", clone.id, additions.len());
            additions.push((clone, label));
        }
        for (record, label) in additions {
            let mut record = record;
            record = record.with_property(options.target.clone(), PropertyValue::Number(label))?;
            current.push(record)?;
        }
        rounds.push(ActiveLearningRound {
            round,
            n_before,
            selected: ids,
            uncertainty: selected_uncertainty,
            labels,
            n_after: current.len(),
        });
        if let Some(updated) = retrain(&current)? {
            model = updated;
        }
    }
    Ok(ActiveLearningReport {
        rounds,
        dataset: current,
        model: model.name().to_owned(),
    })
}

fn clone_dataset(dataset: &Dataset) -> Result<Dataset, Error> {
    let mut output = Dataset::new(dataset.metadata().name.clone());
    *output.metadata_mut() = dataset.metadata().clone();
    for record in dataset.iter() {
        output.push(record.clone())?;
    }
    Ok(output)
}

/// Emits the active-learning workflow document.
#[must_use]
pub fn active_learning_workflow(
    name: &str,
    dataset_path: &str,
    model_path: &str,
    batch: usize,
    rounds: usize,
    checkpoint: Option<&str>,
) -> Value {
    serde_json::json!({
        "workflow_version": 1,
        "name": name,
        "description": "generate, predict, select, compute and append in rounds",
        "steps": [
            {"name": "read", "kind": "read", "path": dataset_path, "format": "json"},
            {"name": "search", "kind": "compute", "verb": "search", "arguments": {
                "kind": "random",
                "chunk": 1,
                "checkpoint": checkpoint,
                "options": {"n_structures": batch}
            }},
            {"name": "predict", "kind": "compute", "verb": "ml_predict", "arguments": {
                "model": model_path
            }},
            {"name": "select", "kind": "compute", "verb": "ml_select", "arguments": {
                "batch": batch, "strategy": "uncertainty", "rounds": rounds
            }},
            {"name": "append", "kind": "transform", "verb": "ml_append", "arguments": {
                "target": "computed"
            }},
            {"name": "write", "kind": "write", "path": "active-learning.json", "format": "json"}
        ]
    })
}
