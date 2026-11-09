//! Bonding, ring, group, neighbor and measurement bindings (Task 1.17.2).

use pyo3::prelude::*;
use pyo3::types::{PyDict, PyList};

use crate::structures::load_structure;

type RingPath = Vec<Vec<(usize, (i32, i32, i32))>>;

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

/// Ring paths with the lattice image of every atom, in ring order.
///
/// Each path is a list of ``(atom, (i, j, k))`` pairs; add the image to the
/// atom's fractional position to lay the ring out in Cartesian space.
#[pyfunction]
#[pyo3(signature = (text, max_size, format_hint = None))]
pub(crate) fn ring_paths(
    text: &str,
    max_size: i64,
    format_hint: Option<&str>,
) -> PyResult<RingPath> {
    let structure = load_structure(text, format_hint)?;
    let size = usize::try_from(max_size).map_err(|_| {
        crate::invalid_input(format!(
            "ring size limit must be a non-negative integer, got {max_size}"
        ))
    })?;
    Ok(xtalforge::core::ring_paths(&structure, size)
        .map_err(crate::map_error)?
        .into_iter()
        .map(|path| {
            path.into_iter()
                .map(|step| (step.atom, (step.image[0], step.image[1], step.image[2])))
                .collect()
        })
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

/// Calculates a PXRD pattern.
#[pyfunction]
#[pyo3(signature = (text, radiation = "cu-ka".to_string(), max_two_theta = 80.0, step = 0.02, format_hint = None))]
#[allow(clippy::too_many_arguments)] // one Python argument per pattern option
pub(crate) fn pxrd<'py>(
    py: Python<'py>,
    text: &str,
    radiation: String,
    max_two_theta: f64,
    step: f64,
    format_hint: Option<&str>,
) -> PyResult<Bound<'py, PyDict>> {
    use xtalforge::analysis::{PxrdOptions, Radiation, calculate_pxrd};
    let structure = load_structure(text, format_hint)?;
    let radiation = match radiation.as_str() {
        "cu-ka" | "cuka" => Radiation::CuKa,
        "cu-ka1" | "cuka1" => Radiation::CuKa1,
        "mo-ka" | "moka" => Radiation::MoKa,
        "co-ka" | "coka" => Radiation::CoKa,
        "fe-ka" | "feka" => Radiation::FeKa,
        "cr-ka" | "crka" => Radiation::CrKa,
        "ag-ka" | "agka" => Radiation::AgKa,
        other => return Err(crate::invalid_input(format!("unknown radiation {other:?}"))),
    };
    let pattern = calculate_pxrd(
        &structure,
        &PxrdOptions {
            radiation,
            two_theta_range: [5.0, max_two_theta],
            step,
            ..PxrdOptions::default()
        },
    )
    .map_err(crate::map_error)?;
    let dict = PyDict::new(py);
    let reflections = PyList::empty(py);
    for reflection in &pattern.reflections {
        let item = PyDict::new(py);
        item.set_item("hkl", reflection.hkl.to_vec())?;
        item.set_item("d", reflection.d)?;
        item.set_item("two_theta", reflection.two_theta)?;
        item.set_item("multiplicity", reflection.multiplicity)?;
        item.set_item("intensity", reflection.intensity)?;
        reflections.append(item)?;
    }
    dict.set_item("reflections", reflections)?;
    dict.set_item(
        "points",
        pattern
            .points
            .iter()
            .map(|point| (point.two_theta, point.intensity))
            .collect::<Vec<_>>(),
    )?;
    Ok(dict)
}
