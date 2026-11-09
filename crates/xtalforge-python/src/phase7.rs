//! Powder refinement bindings (Tasks 7.1-7.11).

use pyo3::prelude::*;
use pyo3::types::{PyDict, PyList};

use xtalforge::analysis::{
    AutomaticOptions, IndexingOptions, PeakOptions, RefineOptions, RietveldOptions,
    extinction_symbol, find_peaks, free_cell_parameters, index_pattern, refine_phase,
    rietx_inventory, run_automatic,
};

use crate::structures::load_structure;

fn pattern(text: &str) -> PyResult<Vec<xtalforge::io::PatternPoint>> {
    xtalforge::io::read_xy(text).map_err(crate::map_error)
}

/// Runs the automatic peak/index/Le Bail/Rietveld pipeline.
#[pyfunction]
#[pyo3(signature = (text, pattern_text, format_hint = None))]
pub(crate) fn powder_automatic<'py>(
    py: Python<'py>,
    text: &str,
    pattern_text: &str,
    format_hint: Option<&str>,
) -> PyResult<Bound<'py, PyDict>> {
    let structure = load_structure(text, format_hint)?;
    let observed = pattern(pattern_text)?;
    let result = run_automatic(&structure, &observed, &AutomaticOptions::default())
        .map_err(crate::map_error)?;
    let dict = PyDict::new(py);
    dict.set_item("peaks", result.peaks.len())?;
    dict.set_item("steps", result.steps)?;
    dict.set_item("index", result.index.map(|index| index.a))?;
    dict.set_item("lebail", result.lebail.map(|lebail| lebail.r_factor))?;
    dict.set_item("rwp", result.phase.map(|phase| phase.rwp))?;
    Ok(dict)
}

/// Refines positions with a whole-pattern least-squares step.
#[pyfunction]
#[pyo3(signature = (text, pattern_text, refine_positions = true, format_hint = None))]
pub(crate) fn powder_refine<'py>(
    py: Python<'py>,
    text: &str,
    pattern_text: &str,
    refine_positions: bool,
    format_hint: Option<&str>,
) -> PyResult<Bound<'py, PyDict>> {
    let structure = load_structure(text, format_hint)?;
    let observed = pattern(pattern_text)?;
    let result = refine_phase(
        &structure,
        &observed,
        &RietveldOptions::default(),
        &RefineOptions {
            refine_positions,
            ..RefineOptions::default()
        },
    )
    .map_err(crate::map_error)?;
    let dict = PyDict::new(py);
    dict.set_item("rwp", result.rwp)?;
    dict.set_item("converged", result.converged)?;
    dict.set_item("iterations", result.iterations)?;
    dict.set_item("parameters", result.parameters)?;
    Ok(dict)
}

/// Returns the free cell parameters and the extinction symbol.
#[pyfunction]
pub(crate) fn powder_cells<'py>(py: Python<'py>, number: u32) -> PyResult<Bound<'py, PyDict>> {
    let dict = PyDict::new(py);
    dict.set_item(
        "free_parameters",
        free_cell_parameters(number).map_err(crate::map_error)?,
    )?;
    dict.set_item(
        "extinction",
        extinction_symbol(number).map_err(crate::map_error)?,
    )?;
    Ok(dict)
}

/// Returns the rietx bridge inventory.
#[pyfunction]
pub(crate) fn rietx<'py>(py: Python<'py>) -> PyResult<Bound<'py, PyDict>> {
    let bridge = rietx_inventory();
    let dict = PyDict::new(py);
    dict.set_item("name", bridge.name)?;
    dict.set_item("capabilities", bridge.capabilities)?;
    dict.set_item("license", bridge.license)?;
    dict.set_item("source", bridge.source)?;
    Ok(dict)
}

/// Returns the RCSR names (kept here for the powder page's cross-links).
#[pyfunction]
pub(crate) fn powder_supported_groups<'py>(py: Python<'py>) -> PyResult<Bound<'py, PyList>> {
    let list = PyList::empty(py);
    for number in [1u32, 2, 3, 4, 14, 15, 16, 75, 143, 168, 195, 225] {
        list.append(number)?;
    }
    Ok(list)
}

/// Finds peaks and indexes the metric families of a measured pattern.
///
/// Returns ``peaks`` and the ranked ``candidates`` with the M20 and F_N
/// figures of merit.
#[pyfunction]
#[pyo3(signature = (pattern_text, wavelength = 1.54184, max_index = 3, tolerance = 0.05, background = false, refine = false))]
pub(crate) fn powder_index<'py>(
    py: Python<'py>,
    pattern_text: &str,
    wavelength: f64,
    max_index: i32,
    tolerance: f64,
    background: bool,
    refine: bool,
) -> PyResult<Bound<'py, PyDict>> {
    let observed = pattern(pattern_text)?;
    let peaks = find_peaks(
        &observed,
        &PeakOptions {
            background,
            refine,
            ..PeakOptions::default()
        },
    )
    .map_err(crate::map_error)?;
    let result = index_pattern(
        &peaks,
        wavelength,
        &IndexingOptions {
            max_index,
            tolerance,
            ..IndexingOptions::default()
        },
    )
    .map_err(crate::map_error)?;
    let dict = PyDict::new(py);
    let peak_list = PyList::empty(py);
    for peak in &peaks {
        let item = PyDict::new(py);
        item.set_item("two_theta", peak.two_theta)?;
        item.set_item("intensity", peak.intensity)?;
        item.set_item("fwhm", peak.fwhm)?;
        peak_list.append(item)?;
    }
    dict.set_item("peaks", peak_list)?;
    let candidates = PyList::empty(py);
    for candidate in &result.candidates {
        let item = PyDict::new(py);
        item.set_item("system", candidate.system.name())?;
        item.set_item("a", candidate.cell.a)?;
        item.set_item("b", candidate.cell.b)?;
        item.set_item("c", candidate.cell.c)?;
        item.set_item("alpha", candidate.cell.alpha)?;
        item.set_item("beta", candidate.cell.beta)?;
        item.set_item("gamma", candidate.cell.gamma)?;
        item.set_item("m20", candidate.m20)?;
        item.set_item("f_n", candidate.f_n)?;
        item.set_item("coverage", candidate.coverage)?;
        item.set_item("completeness", candidate.completeness)?;
        candidates.append(item)?;
    }
    dict.set_item("candidates", candidates)?;
    Ok(dict)
}
