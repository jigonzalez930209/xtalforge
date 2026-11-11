//! Symmetry exports (Task 2.1-2.7).

use wasm_bindgen::prelude::*;

use crate::common::{array_i32, array3, parse_any, set, summary_object, to_js_error};

/// Detects the space group and returns the dataset as a JS object.
///
/// The object carries `number`, `hall_number`, `hm_symbol`, `pearson_symbol`,
/// `operations` (rotation, translation, triplet), `orbits`, `wyckoffs`,
/// `site_symmetry_symbols`, `std_lattice`, `std_positions`,
/// `std_origin_shift` and `std_rotation_matrix`.
///
/// # Errors
///
/// Returns a JavaScript `Error` for invalid text, a bad tolerance or a failed
/// symmetry search.
#[wasm_bindgen]
pub fn detect_symmetry(text: &str, symprec: f64) -> Result<JsValue, JsValue> {
    let structure = parse_any(text).map_err(to_js_error)?;
    let dataset = xtalforge::symmetry::detect_symmetry(
        &structure,
        xtalforge::symmetry::SymmetryOptions {
            symprec,
            ..Default::default()
        },
    )
    .map_err(to_js_error)?;
    let object = js_sys::Object::new();
    set(
        &object,
        "number",
        JsValue::from_f64(f64::from(dataset.number)),
    );
    set(
        &object,
        "hall_number",
        JsValue::from_f64(f64::from(dataset.hall_number)),
    );
    set(&object, "hm_symbol", JsValue::from_str(&dataset.hm_symbol));
    set(
        &object,
        "pearson_symbol",
        JsValue::from_str(&dataset.pearson_symbol),
    );
    let operations = js_sys::Array::new();
    for operation in &dataset.operations {
        let item = js_sys::Object::new();
        let rotation = js_sys::Array::new();
        for row in operation.rotation() {
            rotation.push(&array_i32(row).into());
        }
        set(&item, "rotation", rotation.into());
        set(&item, "translation", array3(operation.translation()).into());
        set(&item, "triplet", JsValue::from_str(&operation.triplet()));
        operations.push(&item.into());
    }
    set(&object, "operations", operations.into());
    let orbits = js_sys::Array::new();
    for orbit in &dataset.orbits {
        orbits.push(&JsValue::from_f64(*orbit as f64));
    }
    set(&object, "orbits", orbits.into());
    let wyckoffs = js_sys::Array::new();
    for wyckoff in &dataset.wyckoffs {
        wyckoffs.push(&JsValue::from_str(&wyckoff.to_string()));
    }
    set(&object, "wyckoffs", wyckoffs.into());
    let symbols = js_sys::Array::new();
    for symbol in &dataset.site_symmetry_symbols {
        symbols.push(&JsValue::from_str(symbol));
    }
    set(&object, "site_symmetry_symbols", symbols.into());
    let lattice = js_sys::Array::new();
    for row in dataset.std_lattice {
        lattice.push(&array3(row).into());
    }
    set(&object, "std_lattice", lattice.into());
    let positions = js_sys::Array::new();
    for point in &dataset.std_positions {
        positions.push(&array3(*point).into());
    }
    set(&object, "std_positions", positions.into());
    set(
        &object,
        "std_origin_shift",
        array3(dataset.std_origin_shift).into(),
    );
    let rotation = js_sys::Array::new();
    for row in dataset.std_rotation_matrix {
        rotation.push(&array3(row).into());
    }
    set(&object, "std_rotation_matrix", rotation.into());
    Ok(object.into())
}

/// Reduces a structure to its detected asymmetric unit.
///
/// # Errors
///
/// Returns a JavaScript `Error` for invalid text or a failed search.
#[wasm_bindgen]
pub fn asymmetrize_info(text: &str, symprec: f64) -> Result<JsValue, JsValue> {
    let structure = parse_any(text).map_err(to_js_error)?;
    let reduced = xtalforge::symmetry::asymmetrize(
        &structure,
        xtalforge::symmetry::SymmetryOptions {
            symprec,
            ..Default::default()
        },
    )
    .map_err(to_js_error)?;
    summary_object(&reduced)
}

/// Expands a structure to P1.
///
/// # Errors
///
/// Returns a JavaScript `Error` for invalid text or an empty structure.
#[wasm_bindgen]
pub fn to_p1_info(text: &str) -> Result<JsValue, JsValue> {
    let structure = parse_any(text).map_err(to_js_error)?;
    let expanded = xtalforge::symmetry::to_p1(&structure).map_err(to_js_error)?;
    summary_object(&expanded)
}

/// Applies a space group in `reinterpret` or `impose` mode.
///
/// The object carries the structure summary plus `set_report` with the
/// input/output site counts and the largest imposed displacement.
///
/// # Errors
///
/// Returns a JavaScript `Error` for invalid text, an unknown symbol or mode.
#[wasm_bindgen]
pub fn set_space_group_info(text: &str, symbol: &str, mode: &str) -> Result<JsValue, JsValue> {
    let structure = parse_any(text).map_err(to_js_error)?;
    let group = xtalforge::symmetry::SpaceGroup::parse(symbol).map_err(to_js_error)?;
    let mode = match mode {
        "reinterpret" => xtalforge::symmetry::SetMode::Reinterpret,
        "impose" => xtalforge::symmetry::SetMode::Impose,
        other => {
            return Err(crate::common::invalid_input(&format!(
                "unknown set mode {other:?}; use reinterpret or impose"
            )));
        }
    };
    let report =
        xtalforge::symmetry::set_space_group(&structure, &group, mode).map_err(to_js_error)?;
    let object = summary_object(&report.structure)?;
    let summary = report.summary(&group);
    let set_report = js_sys::Object::new();
    set(
        &set_report,
        "n_input",
        JsValue::from_f64(summary.n_input as f64),
    );
    set(
        &set_report,
        "n_output",
        JsValue::from_f64(summary.n_output as f64),
    );
    set(
        &set_report,
        "max_displacement",
        JsValue::from_f64(summary.max_displacement),
    );
    set(&set_report, "hall", JsValue::from_str(&summary.hall));
    set(
        &set_report,
        "number",
        JsValue::from_f64(f64::from(summary.number)),
    );
    let object: js_sys::Object = object.unchecked_into();
    set(&object, "set_report", set_report.into());
    Ok(object.into())
}

/// Minimum-image distance between two P1 atoms, in Angstrom.
///
/// # Errors
///
/// Returns a JavaScript `Error` for invalid text or an out-of-range index.
#[wasm_bindgen]
pub fn distance(text: &str, i: u32, j: u32) -> Result<f64, JsValue> {
    let structure = parse_any(text).map_err(to_js_error)?;
    xtalforge::core::distance(&structure, i as usize, j as usize).map_err(to_js_error)
}

/// Angle `i-j-k` in degrees, using minimum images.
///
/// # Errors
///
/// Returns a JavaScript `Error` for invalid text, an out-of-range index or a
/// degenerate angle.
#[wasm_bindgen]
pub fn angle(text: &str, i: u32, j: u32, k: u32) -> Result<f64, JsValue> {
    let structure = parse_any(text).map_err(to_js_error)?;
    xtalforge::core::angle(&structure, i as usize, j as usize, k as usize).map_err(to_js_error)
}

/// Torsion `i-j-k-l` in degrees, using minimum images.
///
/// # Errors
///
/// Returns a JavaScript `Error` for invalid text, an out-of-range index or a
/// collinear segment.
#[wasm_bindgen]
pub fn torsion(text: &str, i: u32, j: u32, k: u32, l: u32) -> Result<f64, JsValue> {
    let structure = parse_any(text).map_err(to_js_error)?;
    xtalforge::core::torsion(&structure, i as usize, j as usize, k as usize, l as usize)
        .map_err(to_js_error)
}

/// Lists the translationengleiche subgroups of a Hall setting.
///
/// # Errors
///
/// Returns a JavaScript `Error` when enumeration fails.
#[wasm_bindgen]
pub fn subgroups_info(hall_number: u32, epsilon: f64) -> Result<js_sys::Array, JsValue> {
    let groups = xtalforge::symmetry::translationengleiche_subgroups(hall_number, epsilon)
        .map_err(to_js_error)?;
    let array = js_sys::Array::new();
    for group in groups {
        let object = js_sys::Object::new();
        set(
            &object,
            "kind",
            JsValue::from_str(&format!("{:?}", group.kind)),
        );
        set(&object, "index", JsValue::from_f64(group.index as f64));
        set(&object, "depth", JsValue::from_f64(group.depth as f64));
        set(
            &object,
            "conjugates",
            JsValue::from_f64(group.conjugates as f64),
        );
        set(
            &object,
            "subgroup_number",
            group
                .subgroup_number
                .map_or(JsValue::NULL, |value| JsValue::from_f64(f64::from(value))),
        );
        set(
            &object,
            "operations",
            JsValue::from_f64(group.operations.len() as f64),
        );
        array.push(&object.into());
    }
    Ok(array)
}

/// Lowers a structure to a subgroup number.
///
/// # Errors
///
/// Returns a JavaScript `Error` for invalid text or a lowering failure.
#[wasm_bindgen]
pub fn lower_to_subgroup_info(
    text: &str,
    subgroup: u32,
    tolerance: f64,
) -> Result<JsValue, JsValue> {
    let structure = parse_any(text).map_err(to_js_error)?;
    let lowered = xtalforge::symmetry::lower_to_subgroup(&structure, subgroup, tolerance)
        .map_err(to_js_error)?;
    summary_object(&lowered)
}
