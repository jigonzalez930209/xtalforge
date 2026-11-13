//! Symmetry detection bindings (Task 2.1.4).

use pyo3::prelude::*;
use pyo3::types::{PyDict, PyList};

use xtalforge::symmetry::{SettingChoice, SymmetryOptions, detect_symmetry as detect};

use crate::structures::load_structure;

/// Detects the space group, operations, orbits and Wyckoff letters.
///
/// Also reports the standardized elements, the primitive standardized
/// lattice and positions, the input-to-standard linear map, the Hall symbol,
/// the setting description and the point group, so the adapter can build a
/// spglib-shaped dataset without a second search.
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
    dict.set_item(
        "std_linear",
        dataset
            .std_linear
            .iter()
            .map(|row| row.to_vec())
            .collect::<Vec<_>>(),
    )?;
    dict.set_item("std_elements", &dataset.std_elements)?;
    dict.set_item(
        "prim_std_lattice",
        dataset
            .prim_std_lattice
            .iter()
            .map(|row| row.to_vec())
            .collect::<Vec<_>>(),
    )?;
    dict.set_item("mapping_to_primitive", &dataset.mapping_std_prim)?;
    let group = xtalforge::symmetry::SpaceGroup::from_hall_number(dataset.hall_number)
        .map_err(crate::map_error)?;
    dict.set_item("hall_symbol", &group.hall_symbol)?;
    dict.set_item("setting", &group.setting)?;
    dict.set_item("point_group", &group.point_group)?;
    Ok(dict)
}

/// Reduces a structure to its detected asymmetric unit.
///
/// The returned dictionary has the same shape as :func:`read_structure`.
#[pyfunction]
#[pyo3(signature = (text, format_hint = None, symprec = 1e-3))]
pub(crate) fn asymmetrize<'py>(
    py: Python<'py>,
    text: &str,
    format_hint: Option<&str>,
    symprec: f64,
) -> PyResult<Bound<'py, PyDict>> {
    let structure = load_structure(text, format_hint)?;
    let reduced = xtalforge::symmetry::asymmetrize(
        &structure,
        SymmetryOptions {
            symprec,
            ..SymmetryOptions::default()
        },
    )
    .map_err(crate::map_error)?;
    crate::structures::structure_payload(py, &reduced)
}

/// Expands a structure to P1.
///
/// The returned dictionary has the same shape as :func:`read_structure`.
#[pyfunction]
#[pyo3(signature = (text, format_hint = None))]
pub(crate) fn to_p1<'py>(
    py: Python<'py>,
    text: &str,
    format_hint: Option<&str>,
) -> PyResult<Bound<'py, PyDict>> {
    let structure = load_structure(text, format_hint)?;
    let expanded = xtalforge::symmetry::to_p1(&structure).map_err(crate::map_error)?;
    crate::structures::structure_payload(py, &expanded)
}

/// Applies a space group in ``reinterpret`` or ``impose`` mode.
///
/// The dictionary also carries a ``set_report`` entry with the input/output
/// site counts and the largest imposed displacement.
#[pyfunction]
#[pyo3(signature = (text, symbol, mode = "reinterpret".to_string(), format_hint = None))]
pub(crate) fn set_space_group<'py>(
    py: Python<'py>,
    text: &str,
    symbol: &str,
    mode: String,
    format_hint: Option<&str>,
) -> PyResult<Bound<'py, PyDict>> {
    let structure = load_structure(text, format_hint)?;
    let group = xtalforge::symmetry::SpaceGroup::parse(symbol).map_err(crate::map_error)?;
    let mode = match mode.as_str() {
        "reinterpret" => xtalforge::symmetry::SetMode::Reinterpret,
        "impose" => xtalforge::symmetry::SetMode::Impose,
        other => {
            return Err(crate::invalid_input(format!(
                "unknown set mode {other:?}; use reinterpret or impose"
            )));
        }
    };
    let report =
        xtalforge::symmetry::set_space_group(&structure, &group, mode).map_err(crate::map_error)?;
    let dict = crate::structures::structure_payload(py, &report.structure)?;
    let summary = report.summary(&group);
    let report_dict = PyDict::new(py);
    report_dict.set_item("n_input", summary.n_input)?;
    report_dict.set_item("n_output", summary.n_output)?;
    report_dict.set_item("max_displacement", summary.max_displacement)?;
    report_dict.set_item("hall", summary.hall)?;
    report_dict.set_item("number", summary.number)?;
    dict.set_item("set_report", report_dict)?;
    Ok(dict)
}

/// Lists the translationengleiche subgroups of a Hall setting.
#[pyfunction]
#[pyo3(signature = (hall_number, epsilon = 1e-5))]
pub(crate) fn subgroups<'py>(
    py: Python<'py>,
    hall_number: u32,
    epsilon: f64,
) -> PyResult<Bound<'py, PyList>> {
    let groups = xtalforge::symmetry::translationengleiche_subgroups(hall_number, epsilon)
        .map_err(crate::map_error)?;
    let items = PyList::empty(py);
    for group in groups {
        let dict = PyDict::new(py);
        dict.set_item("kind", format!("{:?}", group.kind))?;
        dict.set_item("index", group.index)?;
        dict.set_item("depth", group.depth)?;
        dict.set_item("conjugates", group.conjugates)?;
        dict.set_item("subgroup_number", group.subgroup_number)?;
        dict.set_item("subgroup_hall_number", group.subgroup_hall_number)?;
        dict.set_item(
            "operations",
            group
                .operations
                .iter()
                .map(|operation| operation.triplet())
                .collect::<Vec<_>>(),
        )?;
        items.append(dict)?;
    }
    Ok(items)
}

/// Lowers a structure to a subgroup number.
#[pyfunction]
#[pyo3(signature = (text, subgroup, tolerance = 0.5, format_hint = None))]
pub(crate) fn lower_to_subgroup<'py>(
    py: Python<'py>,
    text: &str,
    subgroup: u32,
    tolerance: f64,
    format_hint: Option<&str>,
) -> PyResult<Bound<'py, PyDict>> {
    let structure = load_structure(text, format_hint)?;
    let lowered = xtalforge::symmetry::lower_to_subgroup(&structure, subgroup, tolerance)
        .map_err(crate::map_error)?;
    crate::structures::structure_payload(py, &lowered)
}
