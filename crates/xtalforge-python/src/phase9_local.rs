//! Python bindings for local environments (Task 9.3).

use pyo3::prelude::*;

use xtalforge::analysis::BvsOptions;
use xtalforge::core::{CoordinationStrategy, RdfOptions};

use crate::phase9_analysis::{json_value, load};

/// Periodic Voronoi cells.
#[pyfunction]
#[pyo3(signature = (text, cutoff = 6.0, format_hint = None))]
pub(crate) fn voronoi(text: &str, cutoff: f64, format_hint: Option<&str>) -> PyResult<Py<PyAny>> {
    let report =
        xtalforge::core::voronoi(&load(text, format_hint)?, cutoff).map_err(crate::map_error)?;
    json_value(&report)
}

/// Coordination shell and shape measures.
#[pyfunction]
#[pyo3(signature = (text, index, strategy = "voronoi", cutoff = 5.0, solid_angle = 0.1, format_hint = None))]
#[allow(clippy::too_many_arguments)] // one Python keyword per strategy option
pub(crate) fn coordination(
    text: &str,
    index: usize,
    strategy: &str,
    cutoff: f64,
    solid_angle: f64,
    format_hint: Option<&str>,
) -> PyResult<Py<PyAny>> {
    let strategy = match strategy {
        "cutoff" => CoordinationStrategy::Cutoff(cutoff),
        "voronoi" => CoordinationStrategy::Voronoi {
            cutoff,
            solid_angle_min: solid_angle,
        },
        "crystalnn" | "crystal_nn" => CoordinationStrategy::CrystalNn { cutoff },
        other => {
            return Err(crate::invalid_input(format!(
                "unknown coordination strategy '{other}'"
            )));
        }
    };
    let shell = xtalforge::core::coordination(&load(text, format_hint)?, index, strategy)
        .map_err(crate::map_error)?;
    json_value(&shell)
}

/// Bond-valence sums.
#[pyfunction]
#[pyo3(signature = (text, cutoff = 5.0, format_hint = None))]
pub(crate) fn bond_valence_sums(
    text: &str,
    cutoff: f64,
    format_hint: Option<&str>,
) -> PyResult<Py<PyAny>> {
    let report = xtalforge::analysis::bond_valence_sums(
        &load(text, format_hint)?,
        &BvsOptions {
            cutoff,
            ..BvsOptions::default()
        },
    )
    .map_err(crate::map_error)?;
    json_value(&report)
}

/// Oxidation states from bond valence.
#[pyfunction]
#[pyo3(signature = (text, cutoff = 5.0, format_hint = None))]
pub(crate) fn oxidation_states(
    text: &str,
    cutoff: f64,
    format_hint: Option<&str>,
) -> PyResult<Py<PyAny>> {
    let report = xtalforge::analysis::oxidation_states(
        &load(text, format_hint)?,
        &BvsOptions {
            cutoff,
            ..BvsOptions::default()
        },
    )
    .map_err(crate::map_error)?;
    json_value(&report)
}

/// Total and partial radial distribution functions.
#[pyfunction]
#[pyo3(signature = (text, rmax = 10.0, bin = 0.05, smoothing = 0, format_hint = None))]
pub(crate) fn radial_distribution(
    text: &str,
    rmax: f64,
    bin: f64,
    smoothing: usize,
    format_hint: Option<&str>,
) -> PyResult<Py<PyAny>> {
    let report = xtalforge::core::radial_distribution(
        &load(text, format_hint)?,
        &RdfOptions {
            r_max: rmax,
            bin_width: bin,
            smoothing,
        },
    )
    .map_err(crate::map_error)?;
    json_value(&report)
}

/// Composition, cell and RDF fingerprint.
#[pyfunction]
#[pyo3(signature = (text, rmax = 10.0, bin = 0.1, format_hint = None))]
pub(crate) fn structure_fingerprint(
    text: &str,
    rmax: f64,
    bin: f64,
    format_hint: Option<&str>,
) -> PyResult<Py<PyAny>> {
    let fingerprint = xtalforge::core::fingerprint(
        &load(text, format_hint)?,
        &RdfOptions {
            r_max: rmax,
            bin_width: bin,
            smoothing: 1,
        },
    )
    .map_err(crate::map_error)?;
    json_value(&fingerprint)
}

/// Similarity of two structure fingerprints.
#[pyfunction]
#[pyo3(signature = (first, second, rmax = 10.0, bin = 0.1))]
pub(crate) fn fingerprint_similarity(
    first: &str,
    second: &str,
    rmax: f64,
    bin: f64,
) -> PyResult<f64> {
    let options = RdfOptions {
        r_max: rmax,
        bin_width: bin,
        smoothing: 1,
    };
    let first =
        xtalforge::core::fingerprint(&load(first, None)?, &options).map_err(crate::map_error)?;
    let second =
        xtalforge::core::fingerprint(&load(second, None)?, &options).map_err(crate::map_error)?;
    Ok(first.similarity(&second))
}
