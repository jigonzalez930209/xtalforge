//! Analysis and builder bindings (Tasks 5.1-5.11 and 6.6).

use pyo3::prelude::*;
use pyo3::types::{PyDict, PyList};

use xtalforge::analysis::{
    RcsrIndex, ZeoOptions, default_radii, distance_grid, overlay_net, parse_zeo_channels,
    quotient_graph,
};

use crate::structures::load_structure;

/// Returns the porosity numbers of a structure.
#[pyfunction]
#[pyo3(signature = (text, spacing = 0.35, probe = 1.2, format_hint = None))]
pub(crate) fn porosity<'py>(
    py: Python<'py>,
    text: &str,
    spacing: f64,
    probe: f64,
    format_hint: Option<&str>,
) -> PyResult<Bound<'py, PyDict>> {
    let structure = load_structure(text, format_hint)?;
    let radii = default_radii(&structure).map_err(crate::map_error)?;
    let grid = distance_grid(&structure, &radii, spacing).map_err(crate::map_error)?;
    let dict = PyDict::new(py);
    dict.set_item("cell_volume", grid.cell_volume())?;
    dict.set_item(
        "accessible_volume",
        grid.accessible_volume(probe).map_err(crate::map_error)?,
    )?;
    dict.set_item("pore_volume", grid.pore_volume())?;
    dict.set_item("surface_area", grid.surface_area())?;
    dict.set_item(
        "largest_included_diameter",
        grid.largest_included_diameter(),
    )?;
    Ok(dict)
}

/// Parses Zeo++ `Channels` output.
#[pyfunction]
pub(crate) fn zeo_parse<'py>(py: Python<'py>, text: &str) -> PyResult<Bound<'py, PyDict>> {
    let result = parse_zeo_channels(text).map_err(crate::map_error)?;
    let dict = PyDict::new(py);
    dict.set_item("pore_volume", result.pore_volume)?;
    dict.set_item("accessible_volume", result.accessible_volume)?;
    dict.set_item("surface_area", result.surface_area)?;
    dict.set_item("largest_included_sphere", result.largest_included_sphere)?;
    dict.set_item("largest_free_sphere", result.largest_free_sphere)?;
    dict.set_item("pore_limiting_diameter", result.pore_limiting_diameter)?;
    Ok(dict)
}

/// Returns one dictionary of net invariants per CGD entry.
#[pyfunction]
#[pyo3(signature = (text, tolerance = 0.01))]
pub(crate) fn net_invariants<'py>(
    py: Python<'py>,
    text: &str,
    tolerance: f64,
) -> PyResult<Bound<'py, PyList>> {
    let file = xtalforge::io::read_cgd(text);
    let index = RcsrIndex::from_embedded().map_err(crate::map_error)?;
    let items = PyList::empty(py);
    for entry in &file.entries {
        let graph = quotient_graph(entry, tolerance).map_err(crate::map_error)?;
        let invariants = graph.invariants();
        let dict = PyDict::new(py);
        dict.set_item("name", &entry.name)?;
        dict.set_item("nodes", invariants.nodes)?;
        dict.set_item("edges", invariants.edges)?;
        dict.set_item("cycle_rank", invariants.cycle_rank)?;
        dict.set_item("coordination", invariants.coordination)?;
        dict.set_item(
            "identified",
            identification(py, &index.identify(&graph).map_err(crate::map_error)?)?,
        )?;
        dict.set_item("components", graph.components().len())?;
        items.append(dict)?;
    }
    Ok(items)
}

/// Builds an overlay dictionary from the first CGD entry.
#[pyfunction]
#[pyo3(signature = (text, tolerance = 0.01))]
pub(crate) fn net_overlay<'py>(
    py: Python<'py>,
    text: &str,
    tolerance: f64,
) -> PyResult<Bound<'py, PyDict>> {
    let file = xtalforge::io::read_cgd(text);
    let entry = file
        .entries
        .first()
        .ok_or_else(|| crate::invalid_input("the CGD file has no readable entry"))?;
    let graph = quotient_graph(entry, tolerance).map_err(crate::map_error)?;
    let cell = match entry.cell.as_deref() {
        Some([a, b, c, alpha, beta, gamma, ..]) => {
            let lattice =
                xtalforge::core::Lattice::from_parameters(*a, *b, *c, *alpha, *beta, *gamma)
                    .map_err(crate::map_error)?;
            lattice.matrix()
        }
        _ => [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]],
    };
    let overlay = overlay_net(
        &graph,
        cell,
        xtalforge::analysis::OverlayStyle {
            color: [90, 90, 90],
            radius: 0.2,
        },
    );
    let dict = PyDict::new(py);
    dict.set_item("name", overlay.name)?;
    dict.set_item(
        "cell",
        overlay
            .cell
            .iter()
            .map(|row| row.to_vec())
            .collect::<Vec<_>>(),
    )?;
    dict.set_item(
        "sites",
        overlay
            .sites
            .iter()
            .map(|(point, label)| (point.to_vec(), label.clone()))
            .collect::<Vec<_>>(),
    )?;
    Ok(dict)
}

/// Returns the RCSR names and validates the index.
#[pyfunction]
pub(crate) fn rcsr_names() -> PyResult<Vec<String>> {
    let index = RcsrIndex::from_embedded().map_err(crate::map_error)?;
    index.validate().map_err(crate::map_error)?;
    Ok(index.names().into_iter().map(str::to_owned).collect())
}

/// One identification as a Python dictionary.
fn identification<'py>(
    py: Python<'py>,
    found: &xtalforge::analysis::rcsr::RcsrIdentification,
) -> PyResult<Bound<'py, PyDict>> {
    let dict = PyDict::new(py);
    dict.set_item("verdict", found.verdict())?;
    dict.set_item("names", found.names.clone())?;
    dict.set_item("nearest", found.nearest.clone())?;
    dict.set_item("nearest_by", &found.nearest_by)?;
    dict.set_item("periodicity", found.fingerprint.periodicity)?;
    dict.set_item("sequences", found.fingerprint.sequences.clone())?;
    dict.set_item("symbols", found.fingerprint.symbols.clone())?;
    Ok(dict)
}

/// Runs Zeo++ if the binary is available.
#[pyfunction]
#[pyo3(signature = (text, probe = 1.2, samples = 1000, format_hint = None))]
pub(crate) fn zeo_run<'py>(
    py: Python<'py>,
    text: &str,
    probe: f64,
    samples: u32,
    format_hint: Option<&str>,
) -> PyResult<Bound<'py, PyDict>> {
    let structure = load_structure(text, format_hint)?;
    let result = xtalforge::analysis::run_zeo(
        &structure,
        &ZeoOptions {
            probe_radius: probe,
            samples,
            ..ZeoOptions::default()
        },
    )
    .map_err(crate::map_error)?;
    let dict = PyDict::new(py);
    dict.set_item("pore_volume", result.pore_volume)?;
    dict.set_item("accessible_volume", result.accessible_volume)?;
    dict.set_item("surface_area", result.surface_area)?;
    dict.set_item("largest_included_sphere", result.largest_included_sphere)?;
    dict.set_item("largest_free_sphere", result.largest_free_sphere)?;
    dict.set_item("pore_limiting_diameter", result.pore_limiting_diameter)?;
    Ok(dict)
}

/// Adds hydrogens at the missing valences.
#[pyfunction]
#[pyo3(signature = (text, distance = 1.09, format_hint = None))]
pub(crate) fn add_hydrogens<'py>(
    py: Python<'py>,
    text: &str,
    distance: f64,
    format_hint: Option<&str>,
) -> PyResult<Bound<'py, PyDict>> {
    let structure = load_structure(text, format_hint)?;
    let options = xtalforge::build::hydrogens::HydrogenOptions {
        distance,
        ..Default::default()
    };
    let filled = xtalforge::build::hydrogens::add_hydrogens(&structure, &options)
        .map_err(crate::map_error)?;
    crate::structures::structure_payload(py, &filled)
}
