//! Symmetry detection bindings (Task 2.1.4).

use pyo3::prelude::*;
use pyo3::types::{PyDict, PyList};

use xtalforge::symmetry::{SettingChoice, SymmetryOptions, detect_symmetry as detect};

use crate::structures::load_structure;

/// Detects the space group, operations, orbits and Wyckoff letters.
///
/// Parameters
/// ----------
/// text:
///     Structure text.
/// format_hint:
///     Format name or extension; ``None`` auto-detects it.
/// symprec:
///     Position tolerance in Angstrom.
/// angle_tolerance:
///     Angle tolerance in degrees; ``None`` uses the symprec-driven default.
/// hall_number:
///     Preferred Hall setting (1 to 530); ``None`` uses the spglib default.
#[pyfunction]
#[pyo3(signature = (text, format_hint = None, symprec = 1e-3, angle_tolerance = None, hall_number = None))]
pub(crate) fn detect_symmetry<'py>(
    py: Python<'py>,
    text: &str,
    format_hint: Option<&str>,
    symprec: f64,
    angle_tolerance: Option<f64>,
    hall_number: Option<u32>,
) -> PyResult<Bound<'py, PyDict>> {
    let structure = load_structure(text, format_hint)?;
    let options = SymmetryOptions {
        symprec,
        angle_tolerance: angle_tolerance.map(f64::to_radians),
        setting: hall_number.map_or(SettingChoice::Spglib, SettingChoice::HallNumber),
    };
    let dataset = detect(&structure, options).map_err(crate::map_error)?;
    let dict = PyDict::new(py);
    dict.set_item("number", dataset.number)?;
    dict.set_item("hall_number", dataset.hall_number)?;
    dict.set_item("hm_symbol", &dataset.hm_symbol)?;
    dict.set_item("pearson_symbol", &dataset.pearson_symbol)?;
    let operations = PyList::empty(py);
    for operation in &dataset.operations {
        let item = PyDict::new(py);
        item.set_item(
            "rotation",
            operation
                .rotation()
                .iter()
                .map(|row| row.to_vec())
                .collect::<Vec<_>>(),
        )?;
        item.set_item("translation", operation.translation().to_vec())?;
        item.set_item("triplet", operation.triplet())?;
        operations.append(item)?;
    }
    dict.set_item("operations", operations)?;
    dict.set_item("orbits", &dataset.orbits)?;
    dict.set_item(
        "wyckoffs",
        dataset
            .wyckoffs
            .iter()
            .map(char::to_string)
            .collect::<Vec<_>>(),
    )?;
    dict.set_item("site_symmetry_symbols", &dataset.site_symmetry_symbols)?;
    dict.set_item(
        "std_lattice",
        dataset
            .std_lattice
            .iter()
            .map(|row| row.to_vec())
            .collect::<Vec<_>>(),
    )?;
    dict.set_item(
        "std_positions",
        dataset
            .std_positions
            .iter()
            .map(|point| point.to_vec())
            .collect::<Vec<_>>(),
    )?;
    dict.set_item("std_origin_shift", dataset.std_origin_shift.to_vec())?;
    dict.set_item(
        "std_rotation_matrix",
        dataset
            .std_rotation_matrix
            .iter()
            .map(|row| row.to_vec())
            .collect::<Vec<_>>(),
    )?;
    Ok(dict)
}
