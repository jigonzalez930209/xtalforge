//! MOF catalogue and builder bindings (Tasks 6.11-6.17).

use pyo3::prelude::*;
use pyo3::types::{PyDict, PyList};

use xtalforge::mof::{MofRequest, blocks, build_mof as build, topologies};

/// Lists the topology catalogue.
#[pyfunction]
pub(crate) fn mof_topologies<'py>(py: Python<'py>) -> PyResult<Bound<'py, PyList>> {
    let items = PyList::empty(py);
    for topology in topologies().map_err(crate::map_mof_error)? {
        let dict = PyDict::new(py);
        dict.set_item("name", topology.name)?;
        dict.set_item("nodes", topology.nodes.len())?;
        dict.set_item("edges", topology.edges.len())?;
        items.append(dict)?;
    }
    Ok(items)
}

/// Lists the building blocks.
#[pyfunction]
pub(crate) fn mof_blocks<'py>(py: Python<'py>) -> PyResult<Bound<'py, PyList>> {
    let items = PyList::empty(py);
    for block in blocks().map_err(crate::map_mof_error)? {
        let dict = PyDict::new(py);
        dict.set_item("name", block.name)?;
        dict.set_item("kind", format!("{:?}", block.kind))?;
        dict.set_item("atoms", block.elements.len())?;
        dict.set_item("connections", block.connection_points.len())?;
        items.append(dict)?;
    }
    Ok(items)
}

/// Builds a framework from a topology, node and linker.
#[pyfunction]
#[pyo3(signature = (topology, node, linker))]
pub(crate) fn build_mof<'py>(
    py: Python<'py>,
    topology: &str,
    node: &str,
    linker: &str,
) -> PyResult<Bound<'py, PyDict>> {
    let request = MofRequest::new(topology, node, linker);
    let built = build(&request).map_err(crate::map_mof_error)?;
    let dict = crate::structures::structure_payload(py, &built.structure)?;
    dict.set_item("warnings", built.warnings)?;
    Ok(dict)
}
