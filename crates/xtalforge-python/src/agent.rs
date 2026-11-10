//! Python bindings for the agent layer (Task 8.1).

use pyo3::prelude::*;

use crate::map_agent_error;

/// Returns the runtime capability report as JSON.
#[pyfunction]
pub(crate) fn agent_capabilities() -> PyResult<String> {
    let report = xtalforge::agent::capability_report();
    serde_json::to_string(&report)
        .map_err(|error| pyo3::exceptions::PyRuntimeError::new_err(error.to_string()))
}

/// Returns the verb catalogue as JSON.
#[pyfunction]
pub(crate) fn agent_verbs() -> PyResult<String> {
    serde_json::to_string(xtalforge::agent::catalogue())
        .map_err(|error| pyo3::exceptions::PyRuntimeError::new_err(error.to_string()))
}

/// Runs a JSON script against a fresh session and returns the answers.
#[pyfunction]
pub(crate) fn agent_run(script_json: &str) -> PyResult<String> {
    let answers = xtalforge::agent::run_verbs(script_json).map_err(map_agent_error)?;
    serde_json::to_string(&answers)
        .map_err(|error| pyo3::exceptions::PyRuntimeError::new_err(error.to_string()))
}

/// A stateful session an assistant drives.
#[pyclass]
pub(crate) struct AgentSession {
    session: xtalforge::agent::Session,
}

#[pymethods]
impl AgentSession {
    /// Creates an empty session.
    #[new]
    fn new() -> Self {
        Self {
            session: xtalforge::agent::Session::new(),
        }
    }

    /// Runs one verb with a JSON object and returns the answer as JSON.
    fn run(&mut self, verb: &str, arguments_json: &str) -> PyResult<String> {
        let arguments: serde_json::Value = serde_json::from_str(arguments_json)
            .map_err(|error| pyo3::exceptions::PyValueError::new_err(error.to_string()))?;
        let answer = self
            .session
            .run(verb, &arguments)
            .map_err(map_agent_error)?;
        serde_json::to_string(&answer)
            .map_err(|error| pyo3::exceptions::PyRuntimeError::new_err(error.to_string()))
    }

    /// Names of the open documents as JSON.
    fn documents_json(&self) -> PyResult<String> {
        let summaries = self.session.summaries();
        serde_json::to_string(&summaries)
            .map_err(|error| pyo3::exceptions::PyRuntimeError::new_err(error.to_string()))
    }

    /// One document as structure JSON.
    fn structure_json(&self, name: Option<&str>) -> PyResult<String> {
        let structure = self.session.structure(name).map_err(map_agent_error)?;
        serde_json::to_string(structure)
            .map_err(|error| pyo3::exceptions::PyRuntimeError::new_err(error.to_string()))
    }

    /// Name of the active document.
    fn active_name(&self) -> Option<String> {
        self.session.active_name().map(str::to_owned)
    }
}

pub(crate) fn register(module: &Bound<'_, PyModule>) -> PyResult<()> {
    module.add_function(wrap_pyfunction!(agent_capabilities, module)?)?;
    module.add_function(wrap_pyfunction!(agent_verbs, module)?)?;
    module.add_function(wrap_pyfunction!(agent_run, module)?)?;
    module.add_class::<AgentSession>()?;
    Ok(())
}
