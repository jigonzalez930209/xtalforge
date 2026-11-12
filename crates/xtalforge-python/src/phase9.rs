//! Python bindings for comparison and structure edits (Task 9.1-9.2).

use pyo3::prelude::*;
use pyo3::types::{PyDict, PyList};

use xtalforge::core::{RattleOptions, SiteSort, SpeciesReplacement, Structure};

use crate::structures::{load_structure, structure_payload};

fn load(text: &str, hint: Option<&str>) -> PyResult<Structure> {
    load_structure(text, hint)
}

/// Two structures up to basis, origin and permutation.
#[pyfunction]
#[pyo3(signature = (first, second, format_hint = None, ltol = 0.2, stol = 0.3, angle_tol = 5.0, allow_supercell = true))]
#[allow(clippy::too_many_arguments)] // one Python keyword per matcher option
pub(crate) fn match_structures<'py>(
    py: Python<'py>,
    first: &str,
    second: &str,
    format_hint: Option<&str>,
    ltol: f64,
    stol: f64,
    angle_tol: f64,
    allow_supercell: bool,
) -> PyResult<Bound<'py, PyDict>> {
    let options = xtalforge::symmetry::MatcherOptions {
        ltol,
        stol,
        angle_tol,
        allow_supercell,
        ..xtalforge::symmetry::MatcherOptions::default()
    };
    let outcome = xtalforge::symmetry::match_structures(
        &load(first, format_hint)?,
        &load(second, format_hint)?,
        &options,
    )
    .map_err(crate::map_error)?;
    let dict = PyDict::new(py);
    dict.set_item("matched", outcome.matched)?;
    dict.set_item("rms", outcome.rms)?;
    dict.set_item("n_atoms", outcome.n_atoms)?;
    dict.set_item("mapping", outcome.mapping)?;
    dict.set_item("translation", outcome.translation.to_vec())?;
    dict.set_item("lattice_mapping", outcome.lattice_mapping)?;
    dict.set_item("supercell", outcome.supercell)?;
    Ok(dict)
}

/// Periodic RMSD of two matched structures.
#[pyfunction]
#[pyo3(signature = (first, second, mapping, format_hint = None))]
pub(crate) fn periodic_rmsd(
    first: &str,
    second: &str,
    mapping: Vec<usize>,
    format_hint: Option<&str>,
) -> PyResult<f64> {
    xtalforge::core::periodic_rmsd(
        &load(first, format_hint)?,
        &load(second, format_hint)?,
        &mapping,
    )
    .map_err(crate::map_error)
}

/// Kabsch alignment of two point sets.
#[pyfunction]
#[pyo3(signature = (points, reference))]
pub(crate) fn kabsch_alignment<'py>(
    py: Python<'py>,
    points: Vec<Vec<f64>>,
    reference: Vec<Vec<f64>>,
) -> PyResult<Bound<'py, PyDict>> {
    let points = triples(&points)?;
    let reference = triples(&reference)?;
    let alignment =
        xtalforge::core::kabsch_alignment(&points, &reference).map_err(crate::map_error)?;
    let dict = PyDict::new(py);
    dict.set_item("rotation", alignment.rotation)?;
    dict.set_item("translation", alignment.translation.to_vec())?;
    dict.set_item("rmsd", alignment.rmsd)?;
    dict.set_item("aligned", alignment.aligned)?;
    Ok(dict)
}

/// Deterministic equality classes of a list of structures.
#[pyfunction]
#[pyo3(signature = (structures, ltol = 0.2, stol = 0.3))]
pub(crate) fn group_structures<'py>(
    py: Python<'py>,
    structures: Vec<String>,
    ltol: f64,
    stol: f64,
) -> PyResult<Bound<'py, PyDict>> {
    let mut parsed = Vec::with_capacity(structures.len());
    for text in &structures {
        parsed.push(load(text, None)?);
    }
    let options = xtalforge::symmetry::MatcherOptions {
        ltol,
        stol,
        allow_supercell: true,
        ..xtalforge::symmetry::MatcherOptions::default()
    };
    let groups =
        xtalforge::symmetry::group_structures(&parsed, &options).map_err(crate::map_error)?;
    let dict = PyDict::new(py);
    dict.set_item("n_structures", parsed.len())?;
    dict.set_item("n_groups", groups.len())?;
    dict.set_item("groups", groups)?;
    Ok(dict)
}

/// Applies a strain tensor, deformation gradient, volume or density.
#[pyfunction]
#[pyo3(signature = (text, format_hint = None, voigt = None, gradient = None, volume = None, density = None))]
#[allow(clippy::too_many_arguments)] // one Python keyword per deformation mode
pub(crate) fn deform<'py>(
    py: Python<'py>,
    text: &str,
    format_hint: Option<&str>,
    voigt: Option<Vec<f64>>,
    gradient: Option<Vec<Vec<f64>>>,
    volume: Option<f64>,
    density: Option<f64>,
) -> PyResult<Bound<'py, PyDict>> {
    let mut structure = load(text, format_hint)?;
    if let Some(voigt) = voigt {
        let values: [f64; 6] = voigt
            .as_slice()
            .try_into()
            .map_err(|_| crate::invalid_input("voigt needs six components"))?;
        let tensor = xtalforge::core::voigt_to_tensor(values).map_err(crate::map_error)?;
        structure = xtalforge::core::strain(&structure, tensor).map_err(crate::map_error)?;
    }
    if let Some(gradient) = gradient {
        let [a, b, c] = gradient.as_slice() else {
            return Err(crate::invalid_input("gradient needs three rows"));
        };
        let matrix = [triple(a)?, triple(b)?, triple(c)?];
        structure = xtalforge::core::deform(&structure, matrix).map_err(crate::map_error)?;
    }
    if let Some(volume) = volume {
        structure =
            xtalforge::core::scale_to_volume(&structure, volume).map_err(crate::map_error)?;
    }
    if let Some(density) = density {
        structure =
            xtalforge::core::scale_to_density(&structure, density).map_err(crate::map_error)?;
    }
    structure_payload(py, &structure)
}

/// Random displacements with a fixed seed.
#[pyfunction]
#[pyo3(signature = (text, format_hint = None, amplitude = 0.05, seed = 0, bounded = false))]
pub(crate) fn rattle<'py>(
    py: Python<'py>,
    text: &str,
    format_hint: Option<&str>,
    amplitude: f64,
    seed: u64,
    bounded: bool,
) -> PyResult<Bound<'py, PyDict>> {
    let report = xtalforge::core::rattle(
        &load(text, format_hint)?,
        &RattleOptions {
            amplitude,
            seed,
            bounded,
            ..RattleOptions::default()
        },
    )
    .map_err(crate::map_error)?;
    let dict = structure_payload(py, &report.structure)?;
    dict.set_item("max_displacement", report.max_displacement)?;
    dict.set_item("rms_displacement", report.rms_displacement)?;
    Ok(dict)
}

/// Site edits: add, remove, replace, move, species map, sort and wrap.
#[pyfunction]
#[pyo3(signature = (text, format_hint = None, add = None, remove = None, replace = None, move_sites = None, species = None, fraction = 1.0, seed = 0, sort = None, wrap = false))]
#[allow(clippy::too_many_arguments)] // one Python keyword per editing operation
pub(crate) fn edit_structure<'py>(
    py: Python<'py>,
    text: &str,
    format_hint: Option<&str>,
    add: Option<Vec<(String, Vec<f64>)>>,
    remove: Option<Vec<usize>>,
    replace: Option<Vec<(usize, String)>>,
    move_sites: Option<Vec<(usize, Vec<f64>)>>,
    species: Option<Vec<(String, String)>>,
    fraction: f64,
    seed: u64,
    sort: Option<String>,
    wrap: bool,
) -> PyResult<Bound<'py, PyDict>> {
    let mut structure = load(text, format_hint)?;
    if let Some(species) = species {
        let replacement = SpeciesReplacement {
            map: species.into_iter().collect(),
            fraction,
            seed,
        };
        structure =
            xtalforge::core::replace_species(&structure, &replacement).map_err(crate::map_error)?;
    }
    if let Some(add) = add {
        let mut sites = Vec::with_capacity(add.len());
        for (element, coordinates) in add {
            sites.push((element, triple(&coordinates)?));
        }
        structure = xtalforge::core::add_sites(&structure, &sites).map_err(crate::map_error)?;
    }
    if let Some(remove) = remove {
        structure = xtalforge::core::remove_sites(&structure, &remove).map_err(crate::map_error)?;
    }
    if let Some(replace) = replace {
        for (index, element) in replace {
            structure = xtalforge::core::replace_site(&structure, index, &element)
                .map_err(crate::map_error)?;
        }
    }
    if let Some(move_sites) = move_sites {
        for (index, coordinates) in move_sites {
            structure = xtalforge::core::move_site(&structure, index, triple(&coordinates)?)
                .map_err(crate::map_error)?;
        }
    }
    if let Some(sort) = sort {
        let key = match sort.as_str() {
            "element" => SiteSort::Element,
            "electronegativity" | "en" => SiteSort::Electronegativity,
            "mass" => SiteSort::Mass,
            "label" => SiteSort::Label,
            other => return Err(crate::invalid_input(format!("unknown sort key '{other}'"))),
        };
        structure = xtalforge::core::sort_sites(&structure, key).map_err(crate::map_error)?;
    }
    if wrap {
        structure = xtalforge::core::wrap_sites(&structure).map_err(crate::map_error)?;
    }
    structure_payload(py, &structure)
}

/// Linear interpolation between two matched structures.
#[pyfunction]
#[pyo3(signature = (first, second, images = 5, format_hint = None))]
pub(crate) fn interpolate_structures<'py>(
    py: Python<'py>,
    first: &str,
    second: &str,
    images: u32,
    format_hint: Option<&str>,
) -> PyResult<Bound<'py, PyList>> {
    let options = xtalforge::symmetry::MatcherOptions {
        allow_supercell: false,
        ..xtalforge::symmetry::MatcherOptions::default()
    };
    let frames = xtalforge::symmetry::interpolate_matched(
        &load(first, format_hint)?,
        &load(second, format_hint)?,
        images,
        &options,
    )
    .map_err(crate::map_error)?;
    let items: Vec<Bound<'_, PyDict>> = frames
        .iter()
        .map(|frame| structure_payload(py, frame))
        .collect::<PyResult<_>>()?;
    PyList::new(py, items)
}

/// Alias of :func:`interpolate_structures` with the surface-neutral name.
#[pyfunction]
#[pyo3(signature = (first, second, images = 5, format_hint = None))]
pub(crate) fn interpolate<'py>(
    py: Python<'py>,
    first: &str,
    second: &str,
    images: u32,
    format_hint: Option<&str>,
) -> PyResult<Bound<'py, PyList>> {
    interpolate_structures(py, first, second, images, format_hint)
}

/// Supercell closest to a cube for a target cell count.
#[pyfunction]
#[pyo3(signature = (text, cells, format_hint = None))]
pub(crate) fn optimal_supercell<'py>(
    py: Python<'py>,
    text: &str,
    cells: u32,
    format_hint: Option<&str>,
) -> PyResult<Bound<'py, PyDict>> {
    let choice = xtalforge::core::optimal_supercell(&load(text, format_hint)?, cells)
        .map_err(crate::map_error)?;
    let dict = structure_payload(py, &choice.structure)?;
    dict.set_item("matrix", choice.matrix)?;
    dict.set_item("cells", choice.cells)?;
    dict.set_item("score", choice.score)?;
    Ok(dict)
}

pub(crate) fn triple(values: &[f64]) -> PyResult<[f64; 3]> {
    let [x, y, z] = values else {
        return Err(crate::invalid_input("expected three numbers"));
    };
    Ok([*x, *y, *z])
}

fn triples(values: &[Vec<f64>]) -> PyResult<Vec<[f64; 3]>> {
    values.iter().map(|value| triple(value)).collect()
}
