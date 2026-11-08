//! Project, parameter, document and workspace bindings (Tasks 3.7, 3.8,
//! 3.12 and 3.16).

use pyo3::prelude::*;
use pyo3::types::{PyBytes, PyDict, PyList};

use xtalforge::commands::{
    Command, Document as CoreDocument, Interpenetrate, MoveSites, Origin, SetCharges, SetElement,
};
use xtalforge::modules::{RunRecord, Workspace, archive_records, read_archived_records};

use crate::structures::{load_structure, structure_payload};

/// Returns the JSON schema of a built-in parameter preset.
#[pyfunction]
pub(crate) fn param_schema<'py>(py: Python<'py>, name: &str) -> PyResult<Bound<'py, PyDict>> {
    let set = xtalforge::core::preset(name).map_err(crate::map_error)?;
    let value = set.schema();
    let text = serde_json::to_string(&value)
        .map_err(|error| crate::invalid_input(format!("could not serialize schema: {error}")))?;
    let value: serde_json::Value = serde_json::from_str(&text)
        .map_err(|error| crate::invalid_input(format!("could not parse schema: {error}")))?;
    let json = py.import("json")?;
    let object = json.call_method1("loads", (text.clone(),))?;
    let _ = value;
    let dict = object.cast::<PyDict>()?;
    Ok(dict.clone())
}

/// Packs structure text into `.xtalproj` bytes.
#[pyfunction]
#[pyo3(signature = (text, format_hint = None))]
pub(crate) fn project_pack(
    py: Python<'_>,
    text: &str,
    format_hint: Option<&str>,
) -> PyResult<Py<PyBytes>> {
    let structure = load_structure(text, format_hint)?;
    let graph = structure.bond_graph().map_err(crate::map_error)?;
    let bonds = serde_json::to_value(
        graph
            .bonds
            .iter()
            .map(|bond| {
                serde_json::json!({
                    "i": bond.i,
                    "j": bond.j,
                    "image": bond.image,
                    "distance": bond.distance,
                    "order": bond.order,
                    "explicit": bond.explicit,
                })
            })
            .collect::<Vec<_>>(),
    )
    .map_err(|error| crate::invalid_input(format!("could not serialize bonds: {error}")))?;
    let parts = xtalforge::io::project::ProjectParts::from_structure(
        &structure,
        None,
        Some(bonds),
        None,
        None,
        None,
    )
    .map_err(crate::map_error)?;
    let bytes = xtalforge::io::project::write_project(&parts).map_err(crate::map_error)?;
    Ok(PyBytes::new(py, &bytes).unbind())
}

/// Reads a `.xtalproj` archive and returns its structure payload.
#[pyfunction]
pub(crate) fn project_read<'py>(py: Python<'py>, bytes: &[u8]) -> PyResult<Bound<'py, PyDict>> {
    let parts = xtalforge::io::project::read_project(bytes).map_err(crate::map_error)?;
    let structure = parts.structure().map_err(crate::map_error)?;
    let dict = structure_payload(py, &structure)?;
    dict.set_item("format_version", xtalforge::io::project::FORMAT_VERSION)?;
    dict.set_item("has_bonds", parts.bonds.is_some())?;
    dict.set_item("has_sites", parts.sites.is_some())?;
    dict.set_item("has_view", parts.view.is_some())?;
    dict.set_item("has_session", parts.session.is_some())?;
    Ok(dict)
}

/// An undoable structure document.
#[pyclass(name = "Document", module = "xtalforge")]
pub(crate) struct Document {
    inner: CoreDocument,
}

#[pymethods]
impl Document {
    #[new]
    #[pyo3(signature = (text, format_hint = None))]
    fn new(text: &str, format_hint: Option<&str>) -> PyResult<Self> {
        Ok(Self {
            inner: CoreDocument::new(load_structure(text, format_hint)?),
        })
    }

    /// Moves sites by a fractional delta.
    fn move_sites(&mut self, sites: Vec<usize>, delta: [f64; 3]) -> PyResult<()> {
        let targets: Vec<(usize, [f64; 3])> = sites
            .into_iter()
            .map(|site| {
                let atom =
                    self.inner.structure().sites().get(site).ok_or_else(|| {
                        crate::invalid_input(format!("site {site} is out of range"))
                    })?;
                let frac = atom.frac();
                Ok((
                    site,
                    [frac[0] + delta[0], frac[1] + delta[1], frac[2] + delta[2]],
                ))
            })
            .collect::<PyResult<_>>()?;
        self.apply(Box::new(MoveSites::new(targets)))
    }

    /// Changes the element of one site.
    fn set_element(&mut self, site: usize, element: &str) -> PyResult<()> {
        let command = SetElement::new(vec![site], element).map_err(crate::map_error)?;
        self.apply(Box::new(command))
    }

    /// Sets the charge of every asymmetric site.
    fn set_charges(&mut self, charges: Vec<f64>) -> PyResult<()> {
        let values = charges.into_iter().map(Some).collect();
        self.apply(Box::new(SetCharges::new(values)))
    }

    /// Interpenetrates the structure at a fractional shift.
    fn interpenetrate(&mut self, shift: [f64; 3]) -> PyResult<()> {
        self.apply(Box::new(Interpenetrate::new(shift)))
    }

    fn undo(&mut self) -> PyResult<bool> {
        Ok(self.inner.undo())
    }

    fn redo(&mut self) -> PyResult<bool> {
        Ok(self.inner.redo())
    }

    #[getter]
    fn can_undo(&self) -> bool {
        self.inner.can_undo()
    }

    #[getter]
    fn can_redo(&self) -> bool {
        self.inner.can_redo()
    }

    #[getter]
    fn dirty(&self) -> bool {
        self.inner.dirty()
    }

    #[getter]
    fn history(&self) -> Vec<String> {
        self.inner
            .history()
            .iter()
            .map(|record| record.description.clone())
            .collect()
    }

    fn to_cif(&self) -> PyResult<String> {
        xtalforge::io::write_structure(self.inner.structure(), "cif").map_err(crate::map_error)
    }

    fn to_json(&self) -> PyResult<String> {
        self.inner
            .structure()
            .to_document()
            .to_json()
            .map_err(crate::map_error)
    }

    fn summary<'py>(&self, py: Python<'py>) -> PyResult<Bound<'py, PyDict>> {
        structure_payload(py, self.inner.structure())
    }
}

impl Document {
    fn apply(&mut self, command: Box<dyn Command>) -> PyResult<()> {
        self.inner
            .apply(command, Origin::Script)
            .map_err(crate::map_error)
    }
}

/// Opens an on-disk workspace.
#[pyfunction]
pub(crate) fn workspace_open(path: &str) -> PyResult<String> {
    let workspace = Workspace::open(path).map_err(crate::map_error)?;
    Ok(workspace.root().to_string_lossy().into_owned())
}

/// Writes a run record (JSON) into a workspace.
#[pyfunction]
pub(crate) fn workspace_write_run(path: &str, record_json: &str) -> PyResult<String> {
    let record: RunRecord = serde_json::from_str(record_json)
        .map_err(|error| crate::invalid_input(format!("invalid run record: {error}")))?;
    let workspace = Workspace::open(path).map_err(crate::map_error)?;
    Ok(workspace
        .write_run(&record)
        .map_err(crate::map_error)?
        .to_string_lossy()
        .into_owned())
}

/// Lists run ids in a workspace.
#[pyfunction]
pub(crate) fn workspace_runs(path: &str) -> PyResult<Vec<String>> {
    let workspace = Workspace::open(path).map_err(crate::map_error)?;
    workspace.runs().map_err(crate::map_error)
}

/// Reads one run record as JSON.
#[pyfunction]
pub(crate) fn workspace_read_run(path: &str, id: &str) -> PyResult<String> {
    let workspace = Workspace::open(path).map_err(crate::map_error)?;
    let record = workspace.read_run(id).map_err(crate::map_error)?;
    serde_json::to_string_pretty(&record)
        .map_err(|error| crate::invalid_input(format!("could not serialize run: {error}")))
}

/// Archives run records (JSON list) into an in-memory ZIP.
#[pyfunction]
pub(crate) fn workspace_archive(py: Python<'_>, records_json: &str) -> PyResult<Py<PyBytes>> {
    let records: Vec<RunRecord> = serde_json::from_str(records_json)
        .map_err(|error| crate::invalid_input(format!("invalid run records: {error}")))?;
    let bytes = archive_records(&records).map_err(crate::map_error)?;
    Ok(PyBytes::new(py, &bytes).unbind())
}

/// Reads run records from an in-memory ZIP.
#[pyfunction]
pub(crate) fn workspace_read_archive<'py>(
    py: Python<'py>,
    bytes: &[u8],
) -> PyResult<Bound<'py, PyList>> {
    let records = read_archived_records(bytes).map_err(crate::map_error)?;
    let items = PyList::empty(py);
    for record in records {
        items.append(record.id)?;
    }
    Ok(items)
}
