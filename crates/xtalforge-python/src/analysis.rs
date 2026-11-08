//! Bonding, ring, group, neighbor and measurement bindings (Task 1.17.2).

use pyo3::prelude::*;
use pyo3::types::{PyDict, PyList};

use crate::structures::load_structure;

/// Perceives distance-based bonds on the P1 cell.
///
/// Each dictionary has ``i``, ``j``, ``image``, ``distance``, ``order`` and
/// ``explicit``.
#[pyfunction]
#[pyo3(signature = (text, format_hint = None))]
pub(crate) fn bonds<'py>(
    py: Python<'py>,
    text: &str,
    format_hint: Option<&str>,
) -> PyResult<Bound<'py, PyList>> {
    let structure = load_structure(text, format_hint)?;
    let graph = structure.bond_graph().map_err(crate::map_error)?;
    let mut items = Vec::with_capacity(graph.bonds.len());
    for bond in &graph.bonds {
        let dict = PyDict::new(py);
        dict.set_item("i", bond.i)?;
        dict.set_item("j", bond.j)?;
        dict.set_item("image", bond.image.to_vec())?;
        dict.set_item("distance", bond.distance)?;
        dict.set_item("order", bond.order)?;
        dict.set_item("explicit", bond.explicit)?;
        items.push(dict);
    }
    PyList::new(py, items)
}

/// Finds chordless primitive rings up to ``max_size`` atoms.
///
/// Returns atom indices in the P1 cell, one sorted index list per ring.
#[pyfunction]
#[pyo3(signature = (text, max_size, format_hint = None))]
pub(crate) fn rings(
    text: &str,
    max_size: i64,
    format_hint: Option<&str>,
) -> PyResult<Vec<Vec<usize>>> {
    let structure = load_structure(text, format_hint)?;
    let size = usize::try_from(max_size).map_err(|_| {
        crate::invalid_input(format!(
            "ring size limit must be a non-negative integer, got {max_size}"
        ))
    })?;
    Ok(xtalforge::core::rings(&structure, size)
        .map_err(crate::map_error)?
        .into_iter()
        .map(|ring| ring.atoms)
        .collect())
}

/// Finds functional groups on the perceived bond graph.
///
/// Each dictionary has ``kind`` and ``atoms``.
#[pyfunction]
#[pyo3(signature = (text, format_hint = None))]
pub(crate) fn groups<'py>(
    py: Python<'py>,
    text: &str,
    format_hint: Option<&str>,
) -> PyResult<Bound<'py, PyList>> {
    let structure = load_structure(text, format_hint)?;
    let matches = xtalforge::core::find_groups(&structure).map_err(crate::map_error)?;
    let mut items = Vec::with_capacity(matches.len());
    for group in matches {
        let dict = PyDict::new(py);
        dict.set_item("kind", group.kind.name())?;
        dict.set_item("atoms", group.atoms)?;
        items.push(dict);
    }
    PyList::new(py, items)
}

/// Lists all periodic neighbors within ``cutoff`` Angstrom.
///
/// Each dictionary has ``i``, ``j``, ``image`` and ``distance``.
#[pyfunction]
#[pyo3(signature = (text, cutoff, format_hint = None))]
pub(crate) fn neighbors<'py>(
    py: Python<'py>,
    text: &str,
    cutoff: f64,
    format_hint: Option<&str>,
) -> PyResult<Bound<'py, PyList>> {
    let structure = load_structure(text, format_hint)?;
    let found = xtalforge::core::neighbors(&structure, cutoff).map_err(crate::map_error)?;
    let mut items = Vec::with_capacity(found.len());
    for neighbor in found {
        let dict = PyDict::new(py);
        dict.set_item("i", neighbor.i)?;
        dict.set_item("j", neighbor.j)?;
        dict.set_item("image", neighbor.image.to_vec())?;
        dict.set_item("distance", neighbor.distance)?;
        items.push(dict);
    }
    PyList::new(py, items)
}

/// Minimum-image distance between two P1 atoms in Angstrom.
#[pyfunction]
#[pyo3(signature = (text, i, j, format_hint = None))]
pub(crate) fn distance(text: &str, i: i64, j: i64, format_hint: Option<&str>) -> PyResult<f64> {
    let structure = load_structure(text, format_hint)?;
    xtalforge::core::distance(&structure, atom_index(i)?, atom_index(j)?).map_err(crate::map_error)
}

/// Angle ``i-j-k`` in degrees, using minimum images around ``j``.
#[pyfunction]
#[pyo3(signature = (text, i, j, k, format_hint = None))]
pub(crate) fn angle(
    text: &str,
    i: i64,
    j: i64,
    k: i64,
    format_hint: Option<&str>,
) -> PyResult<f64> {
    let structure = load_structure(text, format_hint)?;
    xtalforge::core::angle(&structure, atom_index(i)?, atom_index(j)?, atom_index(k)?)
        .map_err(crate::map_error)
}

/// Torsion ``i-j-k-l`` in degrees.
#[pyfunction]
#[pyo3(signature = (text, i, j, k, l, format_hint = None))]
pub(crate) fn torsion(
    text: &str,
    i: i64,
    j: i64,
    k: i64,
    l: i64,
    format_hint: Option<&str>,
) -> PyResult<f64> {
    let structure = load_structure(text, format_hint)?;
    xtalforge::core::torsion(
        &structure,
        atom_index(i)?,
        atom_index(j)?,
        atom_index(k)?,
        atom_index(l)?,
    )
    .map_err(crate::map_error)
}

fn atom_index(value: i64) -> PyResult<usize> {
    usize::try_from(value).map_err(|_| {
        crate::invalid_input(format!(
            "atom index must be a non-negative integer, got {value}"
        ))
    })
}
