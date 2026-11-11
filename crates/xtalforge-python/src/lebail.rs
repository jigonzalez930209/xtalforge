//! Le Bail intensity extraction for Python (functional audit 2026-11-11: the
//! operation existed only in the CLI).

use pyo3::prelude::*;
use pyo3::types::PyDict;

use crate::structures::load_structure;

/// Extracts reflection intensities with the Le Bail method.
///
/// Returns ``positions``, ``intensities``, ``intensity_esds``, ``r_factor``,
/// ``rwp`` and ``chi_squared``, the fields of ``xtalforge --json powder lebail``.
#[pyfunction]
#[pyo3(signature = (text, pattern_text, format_hint = None))]
pub(crate) fn powder_lebail<'py>(
    py: Python<'py>,
    text: &str,
    pattern_text: &str,
    format_hint: Option<&str>,
) -> PyResult<Bound<'py, PyDict>> {
    let structure = load_structure(text, format_hint)?;
    let observed = xtalforge::io::read_xy(pattern_text).map_err(crate::map_error)?;
    let result = xtalforge::analysis::le_bail(
        &structure,
        &observed,
        &xtalforge::analysis::LebailOptions::default(),
    )
    .map_err(crate::map_error)?;
    let dict = PyDict::new(py);
    dict.set_item("positions", result.positions)?;
    dict.set_item("intensities", result.intensities)?;
    dict.set_item("intensity_esds", result.intensity_esds)?;
    dict.set_item("r_factor", result.r_factor)?;
    dict.set_item("rwp", result.rwp)?;
    dict.set_item("chi_squared", result.chi_squared)?;
    Ok(dict)
}
