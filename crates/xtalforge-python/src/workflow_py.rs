//! Workflow and provenance bindings (Phase 16).

use std::path::PathBuf;

use pyo3::prelude::*;

use xtalforge::agent::workflow::{RunOptions, Workflow, run_workflow, schema_text};

/// Returns the published workflow JSON Schema.
#[pyfunction]
pub(crate) fn workflow_schema() -> PyResult<String> {
    schema_text().map_err(|error| crate::invalid_input(error.to_string()))
}

/// Parses and validates a workflow, returning its summary as JSON.
#[pyfunction]
pub(crate) fn workflow_validate(text: &str) -> PyResult<String> {
    let workflow =
        Workflow::from_json(text).map_err(|error| crate::invalid_input(error.to_string()))?;
    let summary = workflow
        .validate()
        .map_err(|error| crate::invalid_input(error.to_string()))?;
    serde_json::to_string(&summary).map_err(|error| crate::invalid_input(error.to_string()))
}

/// Runs a workflow over filesystem inputs and returns the report as JSON.
#[pyfunction]
#[pyo3(signature = (text, inputs, jobs = 1, checkpoint = None))]
pub(crate) fn workflow_run(
    text: &str,
    inputs: Vec<String>,
    jobs: usize,
    checkpoint: Option<String>,
) -> PyResult<String> {
    let workflow =
        Workflow::from_json(text).map_err(|error| crate::invalid_input(error.to_string()))?;
    let options = RunOptions {
        inputs: inputs.into_iter().map(PathBuf::from).collect(),
        jobs,
        checkpoint: checkpoint.map(PathBuf::from),
        resume: true,
    };
    let report = run_workflow(&workflow, &options)
        .map_err(|error| crate::invalid_input(error.to_string()))?;
    serde_json::to_string(&report).map_err(|error| crate::invalid_input(error.to_string()))
}

/// Returns the embedded provenance record as JSON, when present.
#[pyfunction]
pub(crate) fn provenance_json(text: &str) -> PyResult<Option<String>> {
    xtalforge::io::provenance(text)
        .map_err(crate::map_error)?
        .map(|record| record.to_json().map_err(crate::map_error))
        .transpose()
}

pub(crate) fn register(module: &Bound<'_, PyModule>) -> PyResult<()> {
    module.add_function(wrap_pyfunction!(workflow_schema, module)?)?;
    module.add_function(wrap_pyfunction!(workflow_validate, module)?)?;
    module.add_function(wrap_pyfunction!(workflow_run, module)?)?;
    module.add_function(wrap_pyfunction!(provenance_json, module)?)?;
    Ok(())
}
