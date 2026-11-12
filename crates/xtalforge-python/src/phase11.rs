//! Python bindings for the Phase 11 usability helpers: one-call IO, SMILES,
//! sketches, fragments, substitution and guest placement.

use pyo3::prelude::*;
use pyo3::types::{PyDict, PyList};

use xtalforge::build::{
    EmbedOptions, embed_smiles, find_fragment, fragments, from_smiles, parse_smiles, place_guest,
    substitute,
};

use crate::structures::{load_structure, structure_payload};

fn load(text: &str, hint: Option<&str>) -> PyResult<xtalforge::core::Structure> {
    load_structure(text, hint)
}

/// Reads a structure from a filesystem path.
#[pyfunction]
pub(crate) fn read<'py>(py: Python<'py>, path: &str) -> PyResult<Bound<'py, PyDict>> {
    let result = xtalforge::read(path).map_err(crate::map_error)?;
    let dict = structure_payload(py, &result.structure)?;
    dict.set_item(
        "warnings",
        crate::phase9_analysis::json_value(&result.warnings)?,
    )?;
    Ok(dict)
}

/// Writes a structure text to a filesystem path.
#[pyfunction]
#[pyo3(signature = (path, text, format_hint = ""))]
pub(crate) fn write<'py>(
    py: Python<'py>,
    path: &str,
    text: &str,
    format_hint: &str,
) -> PyResult<Bound<'py, PyDict>> {
    let structure = load(text, Some(format_hint))?;
    xtalforge::write(path, &structure).map_err(crate::map_error)?;
    let dict = PyDict::new(py);
    dict.set_item("path", path)?;
    dict.set_item(
        "bytes",
        std::fs::metadata(path).map(|meta| meta.len()).unwrap_or(0),
    )?;
    dict.set_item("summary", structure_payload(py, &structure)?)?;
    Ok(dict)
}

/// Parses a SMILES string and embeds it in 3D.
#[pyfunction]
#[pyo3(signature = (smiles, seed = 7, iterations = 1500))]
pub(crate) fn molecule_from_smiles<'py>(
    py: Python<'py>,
    smiles: &str,
    seed: u64,
    iterations: usize,
) -> PyResult<Bound<'py, PyDict>> {
    let graph = parse_smiles(smiles).map_err(crate::map_error)?;
    let structure = embed_smiles(
        &graph,
        &EmbedOptions {
            seed,
            iterations,
            ..EmbedOptions::default()
        },
    )
    .map_err(crate::map_error)?;
    structure_payload(py, &structure)
}

/// Parses a SMILES and writes it back as the canonical depth-first string.
#[pyfunction]
pub(crate) fn sketch_roundtrip(smiles: &str) -> PyResult<String> {
    let sketch = from_smiles(smiles).map_err(crate::map_error)?;
    sketch.to_smiles().map_err(crate::map_error)
}

/// The embedded fragment library.
#[pyfunction]
pub(crate) fn fragment_library() -> PyResult<Py<PyAny>> {
    crate::phase9_analysis::json_value(&fragments().map_err(crate::map_error)?)
}

/// Embeds one fragment by name.
#[pyfunction]
#[pyo3(signature = (name, seed = 7, iterations = 1500))]
pub(crate) fn fragment_structure<'py>(
    py: Python<'py>,
    name: &str,
    seed: u64,
    iterations: usize,
) -> PyResult<Bound<'py, PyDict>> {
    let fragment = find_fragment(name)
        .map_err(crate::map_error)?
        .ok_or_else(|| crate::invalid_input(format!("unknown fragment {name:?}")))?;
    let graph = parse_smiles(&fragment.smiles).map_err(crate::map_error)?;
    let structure = embed_smiles(
        &graph,
        &EmbedOptions {
            seed,
            iterations,
            ..EmbedOptions::default()
        },
    )
    .map_err(crate::map_error)?;
    let dict = structure_payload(py, &structure)?;
    dict.set_item("fragment", fragment.name)?;
    Ok(dict)
}

/// Replaces one asymmetric site with another element after a clearance check.
#[pyfunction]
#[pyo3(signature = (text, index, element, format_hint = None))]
pub(crate) fn substitute_structure<'py>(
    py: Python<'py>,
    text: &str,
    index: usize,
    element: &str,
    format_hint: Option<&str>,
) -> PyResult<Bound<'py, PyDict>> {
    let structure =
        substitute(&load(text, format_hint)?, index, element).map_err(crate::map_error)?;
    structure_payload(py, &structure)
}

/// Places a guest molecule at a Cartesian offset after a clearance check.
#[pyfunction]
#[pyo3(signature = (host, guest, offset = (0.0, 0.0, 0.0), clearance = 1.0, format_hint = None))]
pub(crate) fn place_guest_structure<'py>(
    py: Python<'py>,
    host: &str,
    guest: Vec<(String, Vec<f64>)>,
    offset: (f64, f64, f64),
    clearance: f64,
    format_hint: Option<&str>,
) -> PyResult<Bound<'py, PyDict>> {
    let mut points = Vec::with_capacity(guest.len());
    for (element, point) in &guest {
        let [x, y, z] = point.as_slice() else {
            return Err(crate::invalid_input("guest points need three coordinates"));
        };
        points.push((element.clone(), [*x, *y, *z]));
    }
    let structure = place_guest(
        &load(host, format_hint)?,
        &points,
        [offset.0, offset.1, offset.2],
        clearance,
    )
    .map_err(crate::map_error)?;
    structure_payload(py, &structure)
}

/// Lists the fragment names as a Python list.
#[pyfunction]
pub(crate) fn fragment_names<'py>(py: Python<'py>) -> PyResult<Bound<'py, PyList>> {
    let names: Vec<String> = fragments()
        .map_err(crate::map_error)?
        .into_iter()
        .map(|fragment| fragment.name)
        .collect();
    PyList::new(py, names)
}

/// Emits a Python `DeprecationWarning` pointing at the replacement.
#[pyfunction]
pub(crate) fn deprecation_warning(py: Python<'_>, name: &str, replacement: &str) -> PyResult<()> {
    let message =
        format!("{name} is deprecated; use {replacement} instead; see docs/guide/deprecation.md");
    let warnings = py.import("warnings")?;
    let category = py.import("builtins")?.getattr("DeprecationWarning")?;
    warnings.call_method1("warn", (message, category))?;
    Ok(())
}
