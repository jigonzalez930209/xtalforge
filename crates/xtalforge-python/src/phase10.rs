//! Python bindings for the MOF platform (Phase 10).

use std::collections::BTreeMap;

use pyo3::prelude::*;
use pyo3::types::{PyDict, PyList};

use xtalforge::core::Structure;

use crate::phase9_analysis::json_value;
use crate::structures::{load_structure, structure_payload};

fn load(text: &str, hint: Option<&str>) -> PyResult<Structure> {
    load_structure(text, hint)
}

/// MOF quality checker.
#[pyfunction]
#[pyo3(signature = (text, format_hint = None))]
pub(crate) fn mof_check(text: &str, format_hint: Option<&str>) -> PyResult<Py<PyAny>> {
    json_value(&xtalforge::mof::check_mof(&load(text, format_hint)?).map_err(crate::map_error)?)
}

/// Removes free solvent and optionally counter-ions.
#[pyfunction]
#[pyo3(signature = (text, format_hint = None, remove_ions = false))]
pub(crate) fn mof_solvent_remove(
    text: &str,
    format_hint: Option<&str>,
    remove_ions: bool,
) -> PyResult<Py<PyAny>> {
    let report = xtalforge::mof::remove_solvent(
        &load(text, format_hint)?,
        &xtalforge::mof::SolventOptions {
            remove_ions,
            ..xtalforge::mof::SolventOptions::default()
        },
    )
    .map_err(crate::map_error)?;
    json_value(&report)
}

/// Resolves partial occupancies to one ordered configuration.
#[pyfunction]
#[pyo3(signature = (text, format_hint = None))]
pub(crate) fn mof_disorder<'py>(
    py: Python<'py>,
    text: &str,
    format_hint: Option<&str>,
) -> PyResult<Bound<'py, PyDict>> {
    let structure = xtalforge::mof::resolve_disorder(
        &load(text, format_hint)?,
        &xtalforge::mof::DisorderOptions::default(),
    )
    .map_err(crate::map_error)?;
    structure_payload(py, &structure)
}

/// Metal oxidation states and framework charge balance.
#[pyfunction]
#[pyo3(signature = (text, format_hint = None))]
pub(crate) fn mof_charge_balance(text: &str, format_hint: Option<&str>) -> PyResult<Py<PyAny>> {
    json_value(
        &xtalforge::mof::charge_balance(&load(text, format_hint)?).map_err(crate::map_error)?,
    )
}

/// Deconstruction into nodes, linkers, solvent and net.
#[pyfunction]
#[pyo3(signature = (text, format_hint = None))]
pub(crate) fn mof_deconstruct(text: &str, format_hint: Option<&str>) -> PyResult<Py<PyAny>> {
    json_value(&xtalforge::mof::deconstruct(&load(text, format_hint)?).map_err(crate::map_error)?)
}

/// MOFid and MOFkey identity.
#[pyfunction]
#[pyo3(signature = (text, format_hint = None))]
pub(crate) fn mof_mofid(text: &str, format_hint: Option<&str>) -> PyResult<Py<PyAny>> {
    json_value(&xtalforge::mof::mofid(&load(text, format_hint)?).map_err(crate::map_error)?)
}

/// Coordinatively unsaturated metals.
#[pyfunction]
#[pyo3(signature = (text, format_hint = None))]
pub(crate) fn mof_open_metal_sites(text: &str, format_hint: Option<&str>) -> PyResult<Py<PyAny>> {
    json_value(
        &xtalforge::mof::open_metal_sites(
            &load(text, format_hint)?,
            &xtalforge::mof::OpenMetalOptions::default(),
        )
        .map_err(crate::map_error)?,
    )
}

/// Linker, node and framework descriptors.
#[pyfunction]
#[pyo3(signature = (text, format_hint = None))]
pub(crate) fn mof_descriptors(text: &str, format_hint: Option<&str>) -> PyResult<Py<PyAny>> {
    json_value(
        &xtalforge::mof::framework_descriptors(&load(text, format_hint)?)
            .map_err(crate::map_error)?,
    )
}

/// Native pore descriptors.
#[pyfunction]
#[pyo3(signature = (text, spacing = 0.8, probe = 1.2, format_hint = None))]
pub(crate) fn mof_pores(
    text: &str,
    spacing: f64,
    probe: f64,
    format_hint: Option<&str>,
) -> PyResult<Py<PyAny>> {
    json_value(
        &xtalforge::mof::pore_descriptors(
            &load(text, format_hint)?,
            &xtalforge::mof::PoreOptions {
                spacing,
                probe_radius: probe,
                channel_probe: probe,
            },
        )
        .map_err(crate::map_error)?,
    )
}

/// Embedded rigid guests.
#[pyfunction]
pub(crate) fn mof_guests() -> PyResult<Py<PyAny>> {
    json_value(&xtalforge::mof::guests().map_err(crate::map_error)?)
}

/// Host UFF/EQeq force field.
#[pyfunction]
#[pyo3(signature = (text, format_hint = None))]
pub(crate) fn mof_forcefield(text: &str, format_hint: Option<&str>) -> PyResult<Py<PyAny>> {
    json_value(&xtalforge::mof::force_field(&load(text, format_hint)?).map_err(crate::map_error)?)
}

/// Widom insertion and Henry coefficient.
#[pyfunction]
#[pyo3(signature = (text, guest = "CO2", temperature = 298.0, samples = 1000, seed = 0, format_hint = None))]
#[allow(clippy::too_many_arguments)]
pub(crate) fn mof_widom(
    text: &str,
    guest: &str,
    temperature: f64,
    samples: usize,
    seed: u64,
    format_hint: Option<&str>,
) -> PyResult<Py<PyAny>> {
    json_value(
        &xtalforge::mof::widom_insertion(
            &load(text, format_hint)?,
            guest,
            &xtalforge::mof::WidomOptions {
                temperature,
                samples,
                seed,
            },
        )
        .map_err(crate::map_error)?,
    )
}

/// Single-component GCMC isotherm.
#[pyfunction]
#[pyo3(signature = (text, guest = "CO2", temperature = 298.0, pressures = vec![1.0, 5.0], steps = 1000, seed = 0, format_hint = None))]
#[allow(clippy::too_many_arguments)]
pub(crate) fn mof_gcmc(
    text: &str,
    guest: &str,
    temperature: f64,
    pressures: Vec<f64>,
    steps: usize,
    seed: u64,
    format_hint: Option<&str>,
) -> PyResult<Py<PyAny>> {
    json_value(
        &xtalforge::mof::gcmc_isotherm(
            &load(text, format_hint)?,
            guest,
            &xtalforge::mof::GcmcOptions {
                temperature,
                pressures,
                steps,
                seed,
            },
        )
        .map_err(crate::map_error)?,
    )
}

/// Guest energy grid summary.
#[pyfunction]
#[pyo3(signature = (text, guest = "CO2", spacing = 1.0, format_hint = None))]
pub(crate) fn mof_energy_grid(
    text: &str,
    guest: &str,
    spacing: f64,
    format_hint: Option<&str>,
) -> PyResult<Py<PyAny>> {
    let grid = xtalforge::mof::energy_grid(&load(text, format_hint)?, guest, spacing)
        .map_err(crate::map_error)?;
    json_value(&serde_json::json!({
        "guest": grid.guest,
        "shape": grid.grid.shape,
        "minimum": grid.minimum,
        "maximum": grid.maximum,
    }))
}

/// RASPA input set.
#[pyfunction]
#[pyo3(signature = (text, guest = "CO2", temperature = 298.0, pressure = 1.0, format_hint = None))]
pub(crate) fn mof_raspa(
    text: &str,
    guest: &str,
    temperature: f64,
    pressure: f64,
    format_hint: Option<&str>,
) -> PyResult<Py<PyAny>> {
    json_value(
        &xtalforge::mof::raspa_input(&load(text, format_hint)?, guest, temperature, pressure)
            .map_err(crate::map_error)?,
    )
}

/// Post-synthetic functionalisation.
#[pyfunction]
#[pyo3(signature = (text, group = "amino", positions = vec![], fraction = 1.0, seed = 0, format_hint = None))]
#[allow(clippy::too_many_arguments)]
pub(crate) fn mof_functionalize<'py>(
    py: Python<'py>,
    text: &str,
    group: &str,
    positions: Vec<usize>,
    fraction: f64,
    seed: u64,
    format_hint: Option<&str>,
) -> PyResult<Bound<'py, PyDict>> {
    let group = match group {
        "amino" | "nh2" => xtalforge::mof::FunctionalGroup::Amino,
        "hydroxyl" | "oh" => xtalforge::mof::FunctionalGroup::Hydroxyl,
        "methyl" | "ch3" => xtalforge::mof::FunctionalGroup::Methyl,
        "fluoro" | "f" => xtalforge::mof::FunctionalGroup::Fluorine,
        other => return Err(crate::invalid_input(format!("unknown group {other:?}"))),
    };
    let structure = xtalforge::mof::functionalize(
        &load(text, format_hint)?,
        &xtalforge::mof::FunctionalizeOptions {
            group,
            positions,
            seed,
            fraction,
        },
    )
    .map_err(crate::map_error)?;
    structure_payload(py, &structure)
}

/// Missing linker and cluster defects.
#[pyfunction]
#[pyo3(signature = (text, linker_fraction = 0.0, cluster_fraction = 0.0, cap = "formate", seed = 0, format_hint = None))]
#[allow(clippy::too_many_arguments)]
pub(crate) fn mof_defects<'py>(
    py: Python<'py>,
    text: &str,
    linker_fraction: f64,
    cluster_fraction: f64,
    cap: &str,
    seed: u64,
    format_hint: Option<&str>,
) -> PyResult<Bound<'py, PyDict>> {
    let cap = match cap {
        "formate" => xtalforge::mof::CapSpecies::Formate,
        "water" => xtalforge::mof::CapSpecies::Water,
        "hydroxide" | "oh" => xtalforge::mof::CapSpecies::Hydroxide,
        other => return Err(crate::invalid_input(format!("unknown cap {other:?}"))),
    };
    let report = xtalforge::mof::missing_defects(
        &load(text, format_hint)?,
        &xtalforge::mof::DefectOptions {
            linker_fraction,
            cluster_fraction,
            cap,
            seed,
        },
    )
    .map_err(crate::map_error)?;
    let dict = structure_payload(py, &report.structure)?;
    dict.set_item("removed_linkers", report.removed_linkers)?;
    dict.set_item("removed_clusters", report.removed_clusters)?;
    dict.set_item("caps", report.caps)?;
    Ok(dict)
}

/// n-fold interpenetration.
#[pyfunction]
#[pyo3(signature = (text, folds = 2, format_hint = None))]
pub(crate) fn mof_interpenetrate<'py>(
    py: Python<'py>,
    text: &str,
    folds: u32,
    format_hint: Option<&str>,
) -> PyResult<Bound<'py, PyDict>> {
    let structure = xtalforge::mof::interpenetrate(&load(text, format_hint)?, folds)
        .map_err(crate::map_error)?;
    structure_payload(py, &structure)
}

/// AA/AB layer stacking.
#[pyfunction]
#[pyo3(signature = (text, layers = 2, offset = (0.5, 0.5), format_hint = None))]
pub(crate) fn mof_stack_layers<'py>(
    py: Python<'py>,
    text: &str,
    layers: u32,
    offset: (f64, f64),
    format_hint: Option<&str>,
) -> PyResult<Bound<'py, PyDict>> {
    let structure =
        xtalforge::mof::stack_layers(&load(text, format_hint)?, layers, [offset.0, offset.1])
            .map_err(crate::map_error)?;
    structure_payload(py, &structure)
}

/// Mixed-metal or mixed-linker framework.
#[pyfunction]
#[pyo3(signature = (text, replacements, ordered = true, fraction = 1.0, seed = 0, format_hint = None))]
#[allow(clippy::too_many_arguments)]
pub(crate) fn mof_mixed<'py>(
    py: Python<'py>,
    text: &str,
    replacements: Vec<(String, String)>,
    ordered: bool,
    fraction: f64,
    seed: u64,
    format_hint: Option<&str>,
) -> PyResult<Bound<'py, PyDict>> {
    let mut map = BTreeMap::new();
    for (from, to) in replacements {
        map.insert(from, to);
    }
    let structure = xtalforge::mof::mixed_framework(
        &load(text, format_hint)?,
        &xtalforge::mof::MixtureOptions {
            map,
            ordered,
            fraction,
            seed,
        },
    )
    .map_err(crate::map_error)?;
    structure_payload(py, &structure)
}

/// Reads a CoRE MOF/QMOF directory.
#[pyfunction]
#[pyo3(signature = (path, kind = "core"))]
pub(crate) fn mof_read_database<'py>(
    py: Python<'py>,
    path: &str,
    kind: &str,
) -> PyResult<Bound<'py, PyDict>> {
    let (entries, info) = if kind.eq_ignore_ascii_case("qmof") {
        xtalforge::mof::read_qmof(std::path::Path::new(path))
    } else {
        xtalforge::mof::read_core_mof(std::path::Path::new(path))
    }
    .map_err(crate::map_error)?;
    let dict = PyDict::new(py);
    dict.set_item("info", json_value(&info)?)?;
    let items: Vec<Bound<'_, PyDict>> = entries
        .iter()
        .map(|entry| {
            let item = structure_payload(py, &entry.structure)?;
            item.set_item("id", entry.id.clone())?;
            item.set_item("source", json_value(&entry.source)?)?;
            item.set_item("properties", json_value(&entry.properties)?)?;
            Ok(item)
        })
        .collect::<PyResult<_>>()?;
    dict.set_item("entries", PyList::new(py, items)?)?;
    Ok(dict)
}

/// Downloads a database archive into a cache, or reuses the cached copy.
#[pyfunction]
#[pyo3(signature = (url, cache = ".xtalforge-cache", sha256 = None, offline = false, kind = "core"))]
pub(crate) fn mof_fetch_database(
    url: &str,
    cache: &str,
    sha256: Option<&str>,
    offline: bool,
    kind: &str,
) -> PyResult<String> {
    let mut options = xtalforge::mof::FetchOptions::new(cache).offline(offline);
    if let Some(sha256) = sha256 {
        options = options.with_sha256(sha256);
    }
    let fetched = if kind.eq_ignore_ascii_case("qmof") {
        xtalforge::mof::fetch_qmof(url, &options)
    } else {
        xtalforge::mof::fetch_core_mof(url, &options)
    }
    .map_err(crate::map_error)?;
    serde_json::to_string(&fetched)
        .map_err(|error| crate::invalid_input(format!("could not serialize fetch: {error}")))
}

/// Cleans and screens a database directory.
#[pyfunction]
#[pyo3(signature = (path, guest = "CO2", spacing = 0.8, samples = 100, limit = 0))]
pub(crate) fn mof_screen(
    path: &str,
    guest: &str,
    spacing: f64,
    samples: usize,
    limit: usize,
) -> PyResult<Py<PyAny>> {
    let (entries, _) =
        xtalforge::mof::read_core_mof(std::path::Path::new(path)).map_err(crate::map_error)?;
    let records = xtalforge::mof::screen_database(
        &entries,
        &xtalforge::mof::ScreeningOptions {
            spacing,
            guest: guest.to_owned(),
            samples,
            limit,
            ..xtalforge::mof::ScreeningOptions::default()
        },
    )
    .map_err(crate::map_error)?;
    json_value(&records)
}

/// Deduplicates structures by MOFkey and the matcher.
#[pyfunction]
#[pyo3(signature = (structures))]
pub(crate) fn mof_deduplicate(structures: Vec<String>) -> PyResult<Py<PyAny>> {
    let mut parsed = Vec::with_capacity(structures.len());
    for text in &structures {
        parsed.push(load(text, None)?);
    }
    json_value(&xtalforge::mof::deduplicate(&parsed).map_err(crate::map_error)?)
}
