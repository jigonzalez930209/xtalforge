//! Python-side plugins (Task 8.3.2).
//!
//! Python callables can be registered as modules, calculators and formats.
//! The contract is deliberately JSON in, JSON out so a plugin cannot depend on
//! a Rust type or on the internal object layout:
//!
//! - module: `callable(parameters_json) -> result_json`
//! - calculator: `callable(structure_json) -> {"energy": ..., "forces": [...]}`
//! - format reader: `callable(text) -> structure_json`
//! - format writer: `callable(structure_json) -> text`
//!
//! The plugins live in a [`Plugins`] instance, not in a global registry, so two
//! sessions or two tests do not interfere. A plugin that must be visible to
//! the Rust paths (the optimizer, the CLI, MCP) is a Rust plugin crate; see
//! `examples/plugin`.

use std::collections::BTreeMap;

use pyo3::exceptions::{PyTypeError, PyValueError};
use pyo3::prelude::*;
use serde_json::json;

use crate::InputError;

struct FormatPlugin {
    extensions: Vec<String>,
    reader: Py<PyAny>,
    writer: Option<Py<PyAny>>,
}

/// A registry of Python callables.
#[pyclass]
pub struct Plugins {
    modules: BTreeMap<String, Py<PyAny>>,
    calculators: BTreeMap<String, Py<PyAny>>,
    formats: BTreeMap<String, FormatPlugin>,
}

fn check_name(name: &str) -> PyResult<()> {
    if name.trim().is_empty() {
        return Err(InputError::new_err((
            "invalid_input",
            "a plugin needs a non-empty name",
        )));
    }
    Ok(())
}

fn call_json(py: Python<'_>, callable: &Py<PyAny>, argument: &str) -> PyResult<String> {
    let result = callable.bind(py).call1((argument,))?;
    result
        .extract::<String>()
        .map_err(|_| PyTypeError::new_err("a plugin must return a JSON string"))
}

#[pymethods]
impl Plugins {
    /// Creates an empty plugin registry.
    #[new]
    fn new() -> Self {
        Self {
            modules: BTreeMap::new(),
            calculators: BTreeMap::new(),
            formats: BTreeMap::new(),
        }
    }

    /// Registers a module callable.
    fn register_module(&mut self, name: &str, callable: Py<PyAny>) -> PyResult<()> {
        check_name(name)?;
        if self.modules.contains_key(name) {
            return Err(InputError::new_err((
                "invalid_input",
                format!("module {name} is already registered"),
            )));
        }
        self.modules.insert(name.to_owned(), callable);
        Ok(())
    }

    /// Registers a calculator callable.
    fn register_calculator(&mut self, name: &str, callable: Py<PyAny>) -> PyResult<()> {
        check_name(name)?;
        if self.calculators.contains_key(name) {
            return Err(InputError::new_err((
                "invalid_input",
                format!("calculator {name} is already registered"),
            )));
        }
        self.calculators.insert(name.to_owned(), callable);
        Ok(())
    }

    /// Registers a format reader and an optional writer.
    fn register_format(
        &mut self,
        name: &str,
        extensions: Vec<String>,
        reader: Py<PyAny>,
        writer: Option<Py<PyAny>>,
    ) -> PyResult<()> {
        check_name(name)?;
        if extensions.is_empty() {
            return Err(InputError::new_err((
                "invalid_input",
                format!("format {name} needs at least one extension"),
            )));
        }
        if self.formats.contains_key(name) {
            return Err(InputError::new_err((
                "invalid_input",
                format!("format {name} is already registered"),
            )));
        }
        self.formats.insert(
            name.to_owned(),
            FormatPlugin {
                extensions,
                reader,
                writer,
            },
        );
        Ok(())
    }

    /// Registered names and format extensions as JSON.
    fn names_json(&self) -> String {
        json!({
            "modules": self.modules.keys().collect::<Vec<_>>(),
            "calculators": self.calculators.keys().collect::<Vec<_>>(),
            "formats": self.formats.iter().map(|(name, plugin)| json!({
                "name": name,
                "extensions": plugin.extensions,
                "read": true,
                "write": plugin.writer.is_some(),
            })).collect::<Vec<_>>(),
        })
        .to_string()
    }

    /// Runs a registered module with a JSON parameter object.
    fn run_module(&self, py: Python<'_>, name: &str, parameters_json: &str) -> PyResult<String> {
        let callable = self
            .modules
            .get(name)
            .ok_or_else(|| InputError::new_err(("unknown_module", name.to_owned())))?;
        call_json(py, callable, parameters_json)
    }

    /// Runs a registered calculator on structure JSON.
    fn calculator_energy(
        &self,
        py: Python<'_>,
        name: &str,
        structure_json: &str,
    ) -> PyResult<String> {
        let callable = self
            .calculators
            .get(name)
            .ok_or_else(|| InputError::new_err(("unknown_calculator", name.to_owned())))?;
        call_json(py, callable, structure_json)
    }

    /// Reads text with a registered format.
    fn read_format(&self, py: Python<'_>, name: &str, text: &str) -> PyResult<String> {
        let plugin = self
            .formats
            .get(name)
            .ok_or_else(|| InputError::new_err(("unknown_format", name.to_owned())))?;
        call_json(py, &plugin.reader, text)
    }

    /// Writes structure JSON with a registered format.
    fn write_format(&self, py: Python<'_>, name: &str, structure_json: &str) -> PyResult<String> {
        let plugin = self
            .formats
            .get(name)
            .ok_or_else(|| InputError::new_err(("unknown_format", name.to_owned())))?;
        let writer = plugin
            .writer
            .as_ref()
            .ok_or_else(|| PyValueError::new_err(format!("format {name} is read-only")))?;
        call_json(py, writer, structure_json)
    }
}

pub(crate) fn register(module: &Bound<'_, PyModule>) -> PyResult<()> {
    module.add_class::<Plugins>()?;
    Ok(())
}
