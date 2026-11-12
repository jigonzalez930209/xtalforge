//! Python bindings for local environments, surfaces, defects, interfaces and
//! thermodynamics (Tasks 9.3-9.6).

use pyo3::prelude::*;
use pyo3::types::{PyDict, PyList};

use xtalforge::analysis::{InterfaceOptions, OrderingOptions, SqsOptions};
use xtalforge::core::{EosKind, EosPoint, PhaseEntry, Structure};

use crate::structures::{load_structure, structure_payload};

pub(crate) fn load(text: &str, hint: Option<&str>) -> PyResult<Structure> {
    load_structure(text, hint)
}

fn value_to_py(py: Python<'_>, value: &serde_json::Value) -> PyResult<Py<PyAny>> {
    use serde_json::Value;
    match value {
        Value::Null => Ok(py.None()),
        Value::Bool(flag) => Ok((*flag).into_pyobject(py)?.to_owned().unbind().into_any()),
        Value::Number(number) => {
            if let Some(integer) = number.as_i64() {
                Ok(integer.into_pyobject(py)?.to_owned().unbind().into_any())
            } else {
                Ok(number
                    .as_f64()
                    .unwrap_or(0.0)
                    .into_pyobject(py)?
                    .to_owned()
                    .unbind()
                    .into_any())
            }
        }
        Value::String(text) => Ok(text.into_pyobject(py)?.to_owned().unbind().into_any()),
        Value::Array(items) => {
            let list = PyList::empty(py);
            for item in items {
                list.append(value_to_py(py, item)?)?;
            }
            Ok(list.unbind().into_any())
        }
        Value::Object(map) => {
            let dict = PyDict::new(py);
            for (key, item) in map {
                dict.set_item(key, value_to_py(py, item)?)?;
            }
            Ok(dict.unbind().into_any())
        }
    }
}

pub(crate) fn json_value<T: serde::Serialize>(value: &T) -> PyResult<Py<PyAny>> {
    let value = serde_json::to_value(value)
        .map_err(|error| crate::invalid_input(format!("could not serialize result: {error}")))?;
    Python::attach(|py| value_to_py(py, &value))
}

/// Symmetry-unique Miller indices.
#[pyfunction]
#[pyo3(signature = (text, max_index = 2, symprec = 1e-5, format_hint = None))]
pub(crate) fn unique_miller_indices(
    text: &str,
    max_index: i32,
    symprec: f64,
    format_hint: Option<&str>,
) -> PyResult<Py<PyAny>> {
    let indices =
        xtalforge::analysis::unique_miller_indices(&load(text, format_hint)?, max_index, symprec)
            .map_err(crate::map_error)?;
    json_value(&indices)
}

/// Every termination of a surface.
#[pyfunction]
#[pyo3(signature = (text, h, k, l, layers = 2, vacuum = 10.0, symmetric = false, format_hint = None))]
#[allow(clippy::too_many_arguments)] // one Python keyword per Miller index and option
pub(crate) fn slab_terminations<'py>(
    py: Python<'py>,
    text: &str,
    h: i32,
    k: i32,
    l: i32,
    layers: u32,
    vacuum: f64,
    symmetric: bool,
    format_hint: Option<&str>,
) -> PyResult<Bound<'py, PyList>> {
    let report = xtalforge::analysis::slab_terminations(
        &load(text, format_hint)?,
        [h, k, l],
        layers,
        vacuum,
        symmetric,
    )
    .map_err(crate::map_error)?;
    let items: Vec<Bound<'_, PyDict>> = report
        .terminations
        .iter()
        .map(|termination| {
            let dict = structure_payload(py, &termination.structure)?;
            dict.set_item("index", termination.index)?;
            dict.set_item("shift", termination.shift)?;
            dict.set_item("thickness", termination.thickness)?;
            dict.set_item("dipole", termination.dipole)?;
            dict.set_item("polar", termination.polar)?;
            Ok(dict)
        })
        .collect::<PyResult<_>>()?;
    PyList::new(py, items)
}

/// Wulff construction from surface energies.
#[pyfunction]
#[pyo3(signature = (text, surfaces, format_hint = None))]
pub(crate) fn wulff(
    text: &str,
    surfaces: Vec<Vec<f64>>,
    format_hint: Option<&str>,
) -> PyResult<Py<PyAny>> {
    let structure = load(text, format_hint)?;
    let energies: Vec<xtalforge::core::SurfaceEnergy> = surfaces
        .iter()
        .map(|surface| {
            let [h, k, l, energy] = surface.as_slice() else {
                return Err(crate::invalid_input("each surface needs h,k,l,energy"));
            };
            Ok(xtalforge::core::SurfaceEnergy {
                hkl: [*h as i32, *k as i32, *l as i32],
                energy: *energy,
            })
        })
        .collect::<PyResult<_>>()?;
    let shape = xtalforge::core::wulff(structure.lattice(), &energies).map_err(crate::map_error)?;
    json_value(&shape)
}

/// Top, bridge and hollow adsorption sites.
#[pyfunction]
#[pyo3(signature = (text, tolerance = 0.05, format_hint = None))]
pub(crate) fn adsorption_sites(
    text: &str,
    tolerance: f64,
    format_hint: Option<&str>,
) -> PyResult<Py<PyAny>> {
    let sites = xtalforge::core::adsorption_sites(&load(text, format_hint)?, tolerance)
        .map_err(crate::map_error)?;
    json_value(&sites)
}

/// Places a molecule on one adsorption site.
#[pyfunction]
#[pyo3(signature = (slab, molecule, site_index = 0, offset = 2.0, format_hint = None))]
pub(crate) fn place_molecule<'py>(
    py: Python<'py>,
    slab: &str,
    molecule: &str,
    site_index: usize,
    offset: f64,
    format_hint: Option<&str>,
) -> PyResult<Bound<'py, PyDict>> {
    let slab = load(slab, format_hint)?;
    let molecule = load(molecule, format_hint)?;
    let sites = xtalforge::core::adsorption_sites(&slab, 0.05).map_err(crate::map_error)?;
    let site = sites
        .get(site_index)
        .ok_or_else(|| crate::invalid_input("adsorption site index is out of range"))?;
    let placed = xtalforge::core::place_molecule(&slab, &molecule, site, offset)
        .map_err(crate::map_error)?;
    structure_payload(py, &placed)
}

/// Symmetry-unique vacancies.
#[pyfunction]
#[pyo3(signature = (text, min_image = 10.0, format_hint = None))]
pub(crate) fn vacancies<'py>(
    py: Python<'py>,
    text: &str,
    min_image: f64,
    format_hint: Option<&str>,
) -> PyResult<Bound<'py, PyList>> {
    let defects = xtalforge::analysis::vacancies(&load(text, format_hint)?, min_image)
        .map_err(crate::map_error)?;
    defect_list(py, &defects)
}

/// Symmetry-unique substitutions.
#[pyfunction]
#[pyo3(signature = (text, species, min_image = 10.0, format_hint = None))]
pub(crate) fn substitutions<'py>(
    py: Python<'py>,
    text: &str,
    species: &str,
    min_image: f64,
    format_hint: Option<&str>,
) -> PyResult<Bound<'py, PyList>> {
    let defects = xtalforge::analysis::substitutions(&load(text, format_hint)?, species, min_image)
        .map_err(crate::map_error)?;
    defect_list(py, &defects)
}

/// Interstitial candidates from Voronoi vertices.
#[pyfunction]
#[pyo3(signature = (text, species, min_image = 10.0, cutoff = 5.0, clearance = 1.0, format_hint = None))]
#[allow(clippy::too_many_arguments)] // one Python keyword per defect option
pub(crate) fn interstitials<'py>(
    py: Python<'py>,
    text: &str,
    species: &str,
    min_image: f64,
    cutoff: f64,
    clearance: f64,
    format_hint: Option<&str>,
) -> PyResult<Bound<'py, PyList>> {
    let defects = xtalforge::analysis::interstitials(
        &load(text, format_hint)?,
        species,
        min_image,
        cutoff,
        clearance,
    )
    .map_err(crate::map_error)?;
    defect_list(py, &defects)
}

fn defect_kind_name(kind: xtalforge::analysis::DefectKind) -> &'static str {
    match kind {
        xtalforge::analysis::DefectKind::Vacancy => "vacancy",
        xtalforge::analysis::DefectKind::Substitution => "substitution",
        xtalforge::analysis::DefectKind::Interstitial => "interstitial",
    }
}

fn defect_list<'py>(
    py: Python<'py>,
    defects: &[xtalforge::analysis::PointDefect],
) -> PyResult<Bound<'py, PyList>> {
    let items: Vec<Bound<'_, PyDict>> = defects
        .iter()
        .map(|defect| {
            let dict = structure_payload(py, &defect.structure)?;
            dict.set_item("kind", defect_kind_name(defect.kind))?;
            dict.set_item("site_index", defect.site_index)?;
            dict.set_item("original", defect.original.clone())?;
            dict.set_item("species", defect.species.clone())?;
            dict.set_item("frac", defect.frac.to_vec())?;
            dict.set_item("supercell", defect.supercell)?;
            dict.set_item("image_distance", defect.image_distance)?;
            Ok(dict)
        })
        .collect::<PyResult<_>>()?;
    PyList::new(py, items)
}

/// Orderings of partial occupancies.
#[pyfunction]
#[pyo3(signature = (text, supercell = (1, 1, 1), max_orderings = 64, format_hint = None))]
pub(crate) fn orderings<'py>(
    py: Python<'py>,
    text: &str,
    supercell: (u32, u32, u32),
    max_orderings: usize,
    format_hint: Option<&str>,
) -> PyResult<Bound<'py, PyList>> {
    let results = xtalforge::analysis::orderings(
        &load(text, format_hint)?,
        &OrderingOptions {
            supercell: [supercell.0, supercell.1, supercell.2],
            max_orderings,
            ..OrderingOptions::default()
        },
    )
    .map_err(crate::map_error)?;
    let items: Vec<Bound<'_, PyDict>> = results
        .iter()
        .map(|ordering| {
            let dict = structure_payload(py, &ordering.structure)?;
            dict.set_item("energy", ordering.energy)?;
            dict.set_item("assignment", ordering.assignment.clone())?;
            Ok(dict)
        })
        .collect::<PyResult<_>>()?;
    PyList::new(py, items)
}

/// Special quasirandom structure.
#[pyfunction]
#[pyo3(signature = (text, supercell = (2, 2, 2), seed = 0, cutoff = 6.0, iterations = 5000, format_hint = None))]
#[allow(clippy::too_many_arguments)] // one Python keyword per SQS option
pub(crate) fn special_quasirandom<'py>(
    py: Python<'py>,
    text: &str,
    supercell: (u32, u32, u32),
    seed: u64,
    cutoff: f64,
    iterations: usize,
    format_hint: Option<&str>,
) -> PyResult<Bound<'py, PyDict>> {
    let result = xtalforge::analysis::special_quasirandom(
        &load(text, format_hint)?,
        &SqsOptions {
            supercell: [supercell.0, supercell.1, supercell.2],
            seed,
            cutoff,
            iterations,
            ..SqsOptions::default()
        },
    )
    .map_err(crate::map_error)?;
    let dict = structure_payload(py, &result.structure)?;
    dict.set_item("assignment", result.assignment)?;
    dict.set_item("shells", json_value(&result.shells)?)?;
    dict.set_item("triplet_correlation", result.triplet_correlation)?;
    dict.set_item("objective", result.objective)?;
    dict.set_item("accepted", result.accepted)?;
    dict.set_item("attempted", result.attempted)?;
    Ok(dict)
}

/// Commensurate 2D lattice matches.
#[pyfunction]
#[pyo3(signature = (film, substrate, max_index = 3, max_area = 200.0, tolerance = 0.05, format_hint = None))]
#[allow(clippy::too_many_arguments)] // one Python keyword per interface option
pub(crate) fn lattice_matches(
    film: &str,
    substrate: &str,
    max_index: i32,
    max_area: f64,
    tolerance: f64,
    format_hint: Option<&str>,
) -> PyResult<Py<PyAny>> {
    let matches = xtalforge::analysis::lattice_matches(
        &load(film, format_hint)?,
        &load(substrate, format_hint)?,
        &InterfaceOptions {
            max_index,
            max_area,
            length_tolerance: tolerance,
            ..InterfaceOptions::default()
        },
    )
    .map_err(crate::map_error)?;
    json_value(&matches)
}

/// Stacks a film on a substrate.
#[pyfunction]
#[pyo3(signature = (film, substrate, gap = 2.0, vacuum = 15.0, shift = (0.0, 0.0), format_hint = None))]
#[allow(clippy::too_many_arguments)] // one Python keyword per interface option
pub(crate) fn stack_slabs<'py>(
    py: Python<'py>,
    film: &str,
    substrate: &str,
    gap: f64,
    vacuum: f64,
    shift: (f64, f64),
    format_hint: Option<&str>,
) -> PyResult<Bound<'py, PyDict>> {
    let result = xtalforge::analysis::stack_slabs(
        &load(film, format_hint)?,
        &load(substrate, format_hint)?,
        gap,
        [shift.0, shift.1, 0.0],
        vacuum,
    )
    .map_err(crate::map_error)?;
    let dict = structure_payload(py, &result.structure)?;
    dict.set_item("strain", result.strain)?;
    dict.set_item("gap", result.gap)?;
    dict.set_item("vacuum", result.vacuum)?;
    Ok(dict)
}

/// Elastic tensor from UFF stress finite differences.
#[pyfunction]
#[pyo3(signature = (text, strain = 0.005, format_hint = None))]
pub(crate) fn elastic_tensor(
    text: &str,
    strain: f64,
    format_hint: Option<&str>,
) -> PyResult<Py<PyAny>> {
    let calculator = xtalforge::ff::UffCalculator::new(xtalforge::ff::CalculatorOptions::default())
        .map_err(|error| crate::invalid_input(error.to_string()))?;
    let tensor = xtalforge::ff::elastic_tensor(
        &load(text, format_hint)?,
        &calculator,
        &xtalforge::ff::ElasticOptions { strain },
    )
    .map_err(crate::map_error)?;
    json_value(&tensor)
}

/// Equation-of-state fit.
#[pyfunction]
#[pyo3(signature = (points, model = "birch"))]
pub(crate) fn eos_fit(points: Vec<(f64, f64)>, model: &str) -> PyResult<Py<PyAny>> {
    let kind = match model {
        "birch" | "birch_murnaghan" => EosKind::BirchMurnaghan,
        "vinet" => EosKind::Vinet,
        "murnaghan" => EosKind::Murnaghan,
        other => return Err(crate::invalid_input(format!("unknown EOS model '{other}'"))),
    };
    let points: Vec<EosPoint> = points
        .into_iter()
        .map(|(volume, energy)| EosPoint { volume, energy })
        .collect();
    let fit = xtalforge::core::fit_eos(&points, kind).map_err(crate::map_error)?;
    json_value(&fit)
}

/// Formation energies and convex-hull distances.
#[pyfunction]
#[pyo3(signature = (entries))]
pub(crate) fn convex_hull(
    entries: Vec<(String, std::collections::BTreeMap<String, f64>, f64)>,
) -> PyResult<Py<PyAny>> {
    let entries: Vec<PhaseEntry> = entries
        .into_iter()
        .map(|(name, composition, energy)| PhaseEntry {
            name,
            composition,
            energy,
        })
        .collect();
    let hull = xtalforge::core::convex_hull(&entries, None).map_err(crate::map_error)?;
    json_value(&hull)
}
