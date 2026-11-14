//! Phase 21 exports: ML representations, datasets and inference.

use wasm_bindgen::prelude::*;

use crate::common::{parse_any, to_js_error};

fn json(value: &impl serde::Serialize) -> Result<JsValue, JsValue> {
    let text = serde_json::to_string(value)
        .map_err(|error| JsValue::from_str(&format!("phase 21 JSON: {error}")))?;
    js_sys::JSON::parse(&text).map_err(|_| JsValue::from_str("phase 21 result is not JSON"))
}

/// Builds a periodic crystal graph.
#[wasm_bindgen]
pub fn crystal_graph_info(
    text: &str,
    cutoff: f64,
    max_neighbors: Option<usize>,
) -> Result<JsValue, JsValue> {
    let structure = parse_any(text).map_err(to_js_error)?;
    let graph = xtalforge::ml::crystal_graph_arrays(
        &structure,
        &xtalforge::ml::GraphOptions {
            cutoff,
            max_neighbors,
        },
    )
    .map_err(to_js_error)?;
    json(&graph)
}

/// Calculates the exact SOAP power spectrum.
#[wasm_bindgen]
pub fn soap_descriptor_info(
    text: &str,
    cutoff: f64,
    n_max: usize,
    l_max: usize,
    sigma: f64,
    average: bool,
) -> Result<JsValue, JsValue> {
    let structure = parse_any(text).map_err(to_js_error)?;
    let descriptor = xtalforge::ml::soap_power_spectrum(
        &structure,
        &xtalforge::ml::SoapOptions {
            cutoff,
            n_max,
            l_max,
            sigma,
            average,
        },
    )
    .map_err(to_js_error)?;
    json(&descriptor)
}

/// Calculates ACSF features.
#[wasm_bindgen]
pub fn acsf_descriptor_info(text: &str, cutoff: f64, average: bool) -> Result<JsValue, JsValue> {
    let structure = parse_any(text).map_err(to_js_error)?;
    let descriptor = xtalforge::ml::acsf(
        &structure,
        &xtalforge::ml::AcsfOptions {
            cutoff,
            average,
            ..xtalforge::ml::AcsfOptions::default()
        },
    )
    .map_err(to_js_error)?;
    json(&descriptor)
}

/// Calculates an MBTR descriptor.
#[wasm_bindgen]
pub fn mbtr_descriptor_info(
    text: &str,
    cutoff: f64,
    grid: usize,
    sigma: f64,
    angular: bool,
) -> Result<JsValue, JsValue> {
    let structure = parse_any(text).map_err(to_js_error)?;
    let descriptor = xtalforge::ml::mbtr(
        &structure,
        &xtalforge::ml::MbtrOptions {
            cutoff,
            n_grid: grid,
            sigma,
            angular,
            ..xtalforge::ml::MbtrOptions::default()
        },
    )
    .map_err(to_js_error)?;
    json(&descriptor)
}

/// Calculates MOF-specific features.
#[wasm_bindgen]
pub fn mof_ml_features_info(text: &str, max_path: usize, pores: bool) -> Result<JsValue, JsValue> {
    let structure = parse_any(text).map_err(to_js_error)?;
    let features = xtalforge::ml::mof_features(
        &structure,
        &xtalforge::ml::MofFeaturesOptions {
            max_path,
            pores,
            ..xtalforge::ml::MofFeaturesOptions::default()
        },
    )
    .map_err(to_js_error)?;
    json(&features)
}

/// Splits a JSON Lines dataset.
#[wasm_bindgen]
pub fn ml_split_info(
    dataset_jsonl: &str,
    strategy: &str,
    train: f64,
    validation: f64,
    test: f64,
    seed: f64,
) -> Result<JsValue, JsValue> {
    let dataset = xtalforge::analysis::dataset::Dataset::from_json_lines(dataset_jsonl)
        .map_err(to_js_error)?;
    let strategy = match strategy {
        "random" => xtalforge::ml::SplitStrategy::Random,
        "composition" => xtalforge::ml::SplitStrategy::Composition,
        "cluster" => xtalforge::ml::SplitStrategy::Cluster {
            matcher: Default::default(),
        },
        "prototype" => xtalforge::ml::SplitStrategy::Prototype,
        other => {
            return Err(JsValue::from_str(&format!(
                "unknown split strategy `{other}`"
            )));
        }
    };
    let assignment = xtalforge::ml::split_dataset(
        &dataset,
        &xtalforge::ml::SplitOptions {
            train,
            validation,
            test,
            seed: seed as u64,
            strategy,
        },
    )
    .map_err(to_js_error)?;
    json(&assignment)
}

/// Imports extxyz, JSON or OPTIMADE text.
#[wasm_bindgen]
pub fn ml_import_info(text: &str, format: &str, name: &str) -> Result<JsValue, JsValue> {
    let format = match format {
        "extxyz" => xtalforge::ml::ImportFormat::Extxyz,
        "json" => xtalforge::ml::ImportFormat::Json,
        "optimade" => xtalforge::ml::ImportFormat::Optimade,
        other => {
            return Err(JsValue::from_str(&format!(
                "unknown import format `{other}`"
            )));
        }
    };
    let report =
        xtalforge::ml::import_dataset(text, &xtalforge::ml::ImportOptions::new(format, name))
            .map_err(to_js_error)?;
    let dataset = report.dataset.to_json_lines().map_err(to_js_error)?;
    json(&serde_json::json!({
        "n_records": report.n_records,
        "skipped": report.skipped,
        "dataset": dataset,
    }))
}

/// Predicts a property with a linear model.
#[wasm_bindgen]
pub fn ml_predict_info(model_json: &str, features_json: &str) -> Result<JsValue, JsValue> {
    let model = xtalforge::ml::PropertyModel::from_json(model_json).map_err(to_js_error)?;
    let features: Vec<f64> = serde_json::from_str(features_json)
        .map_err(|error| JsValue::from_str(&format!("feature JSON: {error}")))?;
    let value = xtalforge::ml::Predictor::predict(&model, &features).map_err(to_js_error)?;
    json(&serde_json::json!({ "name": model.name, "value": value }))
}

/// Selects a batch from a feature table.
#[wasm_bindgen]
pub fn ml_select_info(
    features_json: &str,
    uncertainties_json: &str,
    batch: usize,
    strategy: &str,
    diversity_radius: f64,
) -> Result<JsValue, JsValue> {
    let features: Vec<Vec<f64>> = serde_json::from_str(features_json)
        .map_err(|error| JsValue::from_str(&format!("feature JSON: {error}")))?;
    let uncertainties: Vec<f64> = serde_json::from_str(uncertainties_json)
        .map_err(|error| JsValue::from_str(&format!("uncertainty JSON: {error}")))?;
    let strategy = match strategy {
        "uncertainty" => xtalforge::ml::SelectionStrategy::Uncertainty,
        "farthest_point" => xtalforge::ml::SelectionStrategy::FarthestPoint,
        "diversity" => xtalforge::ml::SelectionStrategy::Diversity,
        other => {
            return Err(JsValue::from_str(&format!(
                "unknown selection strategy `{other}`"
            )));
        }
    };
    let selected = xtalforge::ml::select_batch(
        &features,
        &uncertainties,
        &xtalforge::ml::SelectionOptions {
            batch,
            strategy,
            diversity_radius,
        },
    )
    .map_err(to_js_error)?;
    json(&selected)
}

/// Emits an active-learning workflow document.
#[wasm_bindgen]
pub fn ml_active_workflow_info(
    name: &str,
    dataset: &str,
    model: &str,
    batch: usize,
    rounds: usize,
    checkpoint: Option<String>,
) -> Result<JsValue, JsValue> {
    let document = xtalforge::ml::active_learning_workflow(
        name,
        dataset,
        model,
        batch,
        rounds,
        checkpoint.as_deref(),
    );
    json(&document)
}
