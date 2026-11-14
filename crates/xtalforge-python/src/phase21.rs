//! Phase 21 bindings: representations, datasets and inference.

use pyo3::prelude::*;

use xtalforge::analysis::dataset::{Dataset, DatasetRecord, PropertyValue};
use xtalforge::core::StructureDocument;
use xtalforge::ml::{
    AcsfOptions, GraphOptions, MbtrOptions, MofFeaturesOptions, Predictor, PropertyModel,
    SelectionOptions, SelectionStrategy, SplitOptions, SplitStrategy, acsf,
    active_learning_workflow, append_predictions, crystal_graph_arrays, import_dataset, mbtr,
    mof_features, predict_dataset, select_batch, soap_power_spectrum, split_dataset,
};

fn json_value(py: Python<'_>, value: &impl serde::Serialize) -> PyResult<Py<PyAny>> {
    let text = serde_json::to_string(value)
        .map_err(|error| crate::invalid_input(format!("phase 21 JSON: {error}")))?;
    let module = py.import("json")?;
    module
        .call_method1("loads", (text,))
        .map(|value| value.unbind())
}

fn parse_dataset(text: &str) -> PyResult<Dataset> {
    if let Ok(dataset) = Dataset::from_json_lines(text) {
        return Ok(dataset);
    }
    let mut dataset = Dataset::new("python");
    for line in text.lines() {
        if line.trim().is_empty() {
            continue;
        }
        let value: serde_json::Value = serde_json::from_str(line)
            .map_err(|error| crate::invalid_input(format!("dataset line: {error}")))?;
        let id = value
            .get("id")
            .and_then(serde_json::Value::as_str)
            .unwrap_or("record")
            .to_owned();
        let Some(structure_value) = value.get("structure") else {
            continue;
        };
        let structure = if let Some(document) = structure_value
            .get("json")
            .and_then(serde_json::Value::as_str)
        {
            StructureDocument::from_json(document)
                .map_err(crate::map_error)?
                .structure
        } else {
            serde_json::from_value(structure_value.clone())
                .map_err(|error| crate::invalid_input(format!("structure: {error}")))?
        };
        let mut record = DatasetRecord::new(id, structure).map_err(crate::map_error)?;
        if let Some(object) = value
            .get("properties")
            .and_then(serde_json::Value::as_object)
        {
            for (key, value) in object {
                let property = match value {
                    serde_json::Value::Number(number) => {
                        PropertyValue::Number(number.as_f64().unwrap_or(0.0))
                    }
                    serde_json::Value::String(text) => PropertyValue::Text(text.clone()),
                    serde_json::Value::Bool(flag) => PropertyValue::Boolean(*flag),
                    serde_json::Value::Array(values) => {
                        let numbers: Vec<f64> = values
                            .iter()
                            .filter_map(serde_json::Value::as_f64)
                            .collect();
                        PropertyValue::Vector(numbers)
                    }
                    _ => continue,
                };
                record = record
                    .with_property(key.clone(), property)
                    .map_err(crate::map_error)?;
            }
        }
        dataset.push(record).map_err(crate::map_error)?;
    }
    Ok(dataset)
}

fn strategy(name: &str) -> PyResult<SelectionStrategy> {
    match name {
        "uncertainty" => Ok(SelectionStrategy::Uncertainty),
        "farthest_point" => Ok(SelectionStrategy::FarthestPoint),
        "diversity" => Ok(SelectionStrategy::Diversity),
        other => Err(crate::invalid_input(format!(
            "unknown selection strategy `{other}`"
        ))),
    }
}

/// Builds a periodic crystal graph and returns typed arrays.
#[pyfunction(name = "crystal_graph_arrays")]
#[pyo3(signature = (text, cutoff = 5.0, max_neighbors = None, format_hint = None))]
pub(crate) fn crystal_graph_arrays_py(
    py: Python<'_>,
    text: &str,
    cutoff: f64,
    max_neighbors: Option<usize>,
    format_hint: Option<&str>,
) -> PyResult<Py<PyAny>> {
    let structure = crate::structures::load_structure(text, format_hint)?;
    let graph = crystal_graph_arrays(
        &structure,
        &GraphOptions {
            cutoff,
            max_neighbors,
        },
    )
    .map_err(crate::map_error)?;
    json_value(py, &graph)
}

/// Calculates the exact SOAP power spectrum.
#[allow(clippy::too_many_arguments)] // one Python keyword per option
#[pyfunction]
#[pyo3(signature = (text, cutoff = 5.0, n_max = 4, l_max = 4, sigma = 0.5, average = false, format_hint = None))]
pub(crate) fn soap_descriptor(
    py: Python<'_>,
    text: &str,
    cutoff: f64,
    n_max: usize,
    l_max: usize,
    sigma: f64,
    average: bool,
    format_hint: Option<&str>,
) -> PyResult<Py<PyAny>> {
    let structure = crate::structures::load_structure(text, format_hint)?;
    let descriptor = soap_power_spectrum(
        &structure,
        &xtalforge::ml::SoapOptions {
            cutoff,
            n_max,
            l_max,
            sigma,
            average,
        },
    )
    .map_err(crate::map_error)?;
    json_value(py, &descriptor)
}

/// Calculates ACSF features.
#[pyfunction]
#[pyo3(signature = (text, cutoff = 5.0, average = false, format_hint = None))]
pub(crate) fn acsf_descriptor(
    py: Python<'_>,
    text: &str,
    cutoff: f64,
    average: bool,
    format_hint: Option<&str>,
) -> PyResult<Py<PyAny>> {
    let structure = crate::structures::load_structure(text, format_hint)?;
    let descriptor = acsf(
        &structure,
        &AcsfOptions {
            cutoff,
            average,
            ..AcsfOptions::default()
        },
    )
    .map_err(crate::map_error)?;
    json_value(py, &descriptor)
}

/// Calculates an MBTR descriptor.
#[pyfunction]
#[pyo3(signature = (text, cutoff = 5.0, grid = 64, sigma = 0.1, angular = true, format_hint = None))]
pub(crate) fn mbtr_descriptor(
    py: Python<'_>,
    text: &str,
    cutoff: f64,
    grid: usize,
    sigma: f64,
    angular: bool,
    format_hint: Option<&str>,
) -> PyResult<Py<PyAny>> {
    let structure = crate::structures::load_structure(text, format_hint)?;
    let descriptor = mbtr(
        &structure,
        &MbtrOptions {
            cutoff,
            n_grid: grid,
            sigma,
            angular,
            ..MbtrOptions::default()
        },
    )
    .map_err(crate::map_error)?;
    json_value(py, &descriptor)
}

/// Calculates MOF-specific RAC and pore features.
#[pyfunction]
#[pyo3(signature = (text, max_path = 3, pores = true, format_hint = None))]
pub(crate) fn mof_ml_features(
    py: Python<'_>,
    text: &str,
    max_path: usize,
    pores: bool,
    format_hint: Option<&str>,
) -> PyResult<Py<PyAny>> {
    let structure = crate::structures::load_structure(text, format_hint)?;
    let features = mof_features(
        &structure,
        &MofFeaturesOptions {
            max_path,
            pores,
            ..MofFeaturesOptions::default()
        },
    )
    .map_err(crate::map_error)?;
    json_value(py, &features)
}

/// Splits a JSON Lines dataset.
#[pyfunction]
#[pyo3(signature = (text, strategy = "random", train = 0.6, validation = 0.2, test = 0.2, seed = 0))]
pub(crate) fn ml_split(
    py: Python<'_>,
    text: &str,
    strategy: &str,
    train: f64,
    validation: f64,
    test: f64,
    seed: u64,
) -> PyResult<Py<PyAny>> {
    let dataset = parse_dataset(text)?;
    let strategy = match strategy {
        "random" => SplitStrategy::Random,
        "composition" => SplitStrategy::Composition,
        "cluster" => SplitStrategy::Cluster {
            matcher: Default::default(),
        },
        "prototype" => SplitStrategy::Prototype,
        other => {
            return Err(crate::invalid_input(format!(
                "unknown split strategy `{other}`"
            )));
        }
    };
    let assignment = split_dataset(
        &dataset,
        &SplitOptions {
            train,
            validation,
            test,
            seed,
            strategy,
        },
    )
    .map_err(crate::map_error)?;
    json_value(py, &assignment)
}

/// Imports extxyz, JSON or OPTIMADE text into a JSON Lines dataset.
#[pyfunction]
#[pyo3(signature = (text, format = "extxyz", name = "imported"))]
pub(crate) fn ml_import(
    py: Python<'_>,
    text: &str,
    format: &str,
    name: &str,
) -> PyResult<Py<PyAny>> {
    let format = match format {
        "extxyz" => xtalforge::ml::ImportFormat::Extxyz,
        "json" => xtalforge::ml::ImportFormat::Json,
        "optimade" => xtalforge::ml::ImportFormat::Optimade,
        other => {
            return Err(crate::invalid_input(format!(
                "unknown import format `{other}`"
            )));
        }
    };
    let report = import_dataset(text, &xtalforge::ml::ImportOptions::new(format, name))
        .map_err(crate::map_error)?;
    let jsonl = report.dataset.to_json_lines().map_err(crate::map_error)?;
    let value = serde_json::json!({
        "n_records": report.n_records,
        "skipped": report.skipped,
        "dataset": jsonl,
    });
    json_value(py, &value)
}

/// Predicts a property with a linear model.
#[pyfunction]
#[pyo3(signature = (model_json, features))]
pub(crate) fn ml_predict(
    py: Python<'_>,
    model_json: &str,
    features: Vec<f64>,
) -> PyResult<Py<PyAny>> {
    let model = PropertyModel::from_json(model_json).map_err(crate::map_error)?;
    let value = model.predict(&features).map_err(crate::map_error)?;
    json_value(
        py,
        &serde_json::json!({ "name": model.name, "value": value }),
    )
}

/// Predicts every record in a JSON Lines dataset.
#[pyfunction]
#[pyo3(signature = (dataset_jsonl, model_json, target = "predicted", feature = "soap"))]
pub(crate) fn ml_predict_dataset(
    py: Python<'_>,
    dataset_jsonl: &str,
    model_json: &str,
    target: &str,
    feature: &str,
) -> PyResult<Py<PyAny>> {
    let dataset = parse_dataset(dataset_jsonl)?;
    let model = PropertyModel::from_json(model_json).map_err(crate::map_error)?;
    let feature = feature.to_owned();
    let prediction = predict_dataset(&dataset, &model, |structure| {
        feature_vector(structure, &feature).map_err(xtalforge::core::Error::InvalidInput)
    })
    .map_err(crate::map_error)?;
    let updated = append_predictions(&dataset, &prediction, target).map_err(crate::map_error)?;
    let jsonl = updated.to_json_lines().map_err(crate::map_error)?;
    json_value(
        py,
        &serde_json::json!({ "dataset": jsonl, "prediction": prediction }),
    )
}

fn feature_vector(structure: &xtalforge::core::Structure, kind: &str) -> Result<Vec<f64>, String> {
    match kind {
        "soap" => soap_power_spectrum(structure, &Default::default())
            .map_err(|error| error.message())?
            .rows
            .first()
            .cloned()
            .ok_or_else(|| "SOAP produced no row".to_owned()),
        "acsf" => acsf(
            structure,
            &AcsfOptions {
                average: true,
                ..AcsfOptions::default()
            },
        )
        .map_err(|error| error.message())?
        .rows
        .first()
        .cloned()
        .ok_or_else(|| "ACSF produced no row".to_owned()),
        "mbtr" => Ok(mbtr(structure, &MbtrOptions::default())
            .map_err(|error| error.message())?
            .features),
        other => Err(format!("unknown feature kind `{other}`")),
    }
}

/// Selects a batch from a feature table.
#[pyfunction]
#[pyo3(signature = (features, uncertainties, batch = 8, strategy = "uncertainty", diversity_radius = 0.0))]
pub(crate) fn ml_select(
    py: Python<'_>,
    features: Vec<Vec<f64>>,
    uncertainties: Vec<f64>,
    batch: usize,
    strategy: &str,
    diversity_radius: f64,
) -> PyResult<Py<PyAny>> {
    let strategy = strategy_fn(strategy)?;
    let selected = select_batch(
        &features,
        &uncertainties,
        &SelectionOptions {
            batch,
            strategy,
            diversity_radius,
        },
    )
    .map_err(crate::map_error)?;
    json_value(py, &selected)
}

fn strategy_fn(name: &str) -> PyResult<SelectionStrategy> {
    strategy(name)
}

/// Emits an active-learning workflow document.
#[pyfunction]
#[pyo3(signature = (name = "active-learning", dataset = "dataset.jsonl", model = "model.json", batch = 8, rounds = 4, checkpoint = None))]
pub(crate) fn ml_active_workflow(
    py: Python<'_>,
    name: &str,
    dataset: &str,
    model: &str,
    batch: usize,
    rounds: usize,
    checkpoint: Option<&str>,
) -> PyResult<Py<PyAny>> {
    let document = active_learning_workflow(name, dataset, model, batch, rounds, checkpoint);
    json_value(py, &document)
}

/// Registers the Phase 21 functions.
pub(crate) fn register(module: &Bound<'_, PyModule>) -> PyResult<()> {
    module.add_function(wrap_pyfunction!(crystal_graph_arrays_py, module)?)?;
    module.add_function(wrap_pyfunction!(soap_descriptor, module)?)?;
    module.add_function(wrap_pyfunction!(acsf_descriptor, module)?)?;
    module.add_function(wrap_pyfunction!(mbtr_descriptor, module)?)?;
    module.add_function(wrap_pyfunction!(mof_ml_features, module)?)?;
    module.add_function(wrap_pyfunction!(ml_split, module)?)?;
    module.add_function(wrap_pyfunction!(ml_import, module)?)?;
    module.add_function(wrap_pyfunction!(ml_predict, module)?)?;
    module.add_function(wrap_pyfunction!(ml_predict_dataset, module)?)?;
    module.add_function(wrap_pyfunction!(ml_select, module)?)?;
    module.add_function(wrap_pyfunction!(ml_active_workflow, module)?)?;
    Ok(())
}
