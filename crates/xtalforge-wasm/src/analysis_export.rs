//! Analysis and builder exports (Tasks 5.11.2 and 6.6).

use wasm_bindgen::prelude::*;

use crate::common::{parse_any, set, summary_object, to_js_error};

/// Returns the porosity numbers of a structure.
///
/// # Errors
///
/// Returns a JavaScript `Error` for invalid text or a bad spacing.
#[wasm_bindgen]
pub fn porosity_info(text: &str, spacing: f64, probe: f64) -> Result<JsValue, JsValue> {
    let structure = parse_any(text).map_err(to_js_error)?;
    let radii = xtalforge::analysis::default_radii(&structure).map_err(to_js_error)?;
    let grid =
        xtalforge::analysis::distance_grid(&structure, &radii, spacing).map_err(to_js_error)?;
    let object = js_sys::Object::new();
    set(
        &object,
        "cell_volume",
        JsValue::from_f64(grid.cell_volume()),
    );
    set(
        &object,
        "accessible_volume",
        JsValue::from_f64(grid.accessible_volume(probe).map_err(to_js_error)?),
    );
    set(
        &object,
        "pore_volume",
        JsValue::from_f64(grid.pore_volume()),
    );
    set(
        &object,
        "surface_area",
        JsValue::from_f64(grid.surface_area()),
    );
    set(
        &object,
        "largest_included_diameter",
        JsValue::from_f64(grid.largest_included_diameter()),
    );
    Ok(object.into())
}

/// Parses Zeo++ `Channels` output.
///
/// # Errors
///
/// Returns a JavaScript `Error` for missing lines.
#[wasm_bindgen]
pub fn zeo_parse_info(text: &str) -> Result<JsValue, JsValue> {
    let result = xtalforge::analysis::parse_zeo_channels(text).map_err(to_js_error)?;
    let object = js_sys::Object::new();
    set(
        &object,
        "pore_volume",
        JsValue::from_f64(result.pore_volume),
    );
    set(
        &object,
        "accessible_volume",
        JsValue::from_f64(result.accessible_volume),
    );
    set(
        &object,
        "surface_area",
        JsValue::from_f64(result.surface_area),
    );
    set(
        &object,
        "largest_included_sphere",
        JsValue::from_f64(result.largest_included_sphere),
    );
    set(
        &object,
        "largest_free_sphere",
        JsValue::from_f64(result.largest_free_sphere),
    );
    Ok(object.into())
}

/// Returns the invariants of every CGD entry.
///
/// # Errors
///
/// Returns a JavaScript `Error` for an invalid net.
#[wasm_bindgen]
pub fn net_invariants_info(text: &str, tolerance: f64) -> Result<js_sys::Array, JsValue> {
    let file = xtalforge::io::read_cgd(text);
    let index = xtalforge::analysis::RcsrIndex::from_embedded().map_err(to_js_error)?;
    let array = js_sys::Array::new();
    for entry in &file.entries {
        let graph = xtalforge::analysis::quotient_graph(entry, tolerance).map_err(to_js_error)?;
        let invariants = graph.invariants();
        let object = js_sys::Object::new();
        set(&object, "name", JsValue::from_str(&entry.name));
        set(&object, "nodes", JsValue::from_f64(invariants.nodes as f64));
        set(&object, "edges", JsValue::from_f64(invariants.edges as f64));
        set(
            &object,
            "cycle_rank",
            JsValue::from_f64(invariants.cycle_rank as f64),
        );
        let coordination = js_sys::Array::new();
        for value in &invariants.coordination {
            coordination.push(&JsValue::from_f64(*value as f64));
        }
        set(&object, "coordination", coordination.into());
        set(
            &object,
            "components",
            JsValue::from_f64(graph.components().len() as f64),
        );
        let found = index.identify(&graph).map_err(to_js_error)?;
        let identified = js_sys::Object::new();
        set(&identified, "verdict", JsValue::from_str(found.verdict()));
        let names = js_sys::Array::new();
        for name in &found.names {
            names.push(&JsValue::from_str(name));
        }
        set(&identified, "names", names.into());
        let nearest = js_sys::Array::new();
        for name in &found.nearest {
            nearest.push(&JsValue::from_str(name));
        }
        set(&identified, "nearest", nearest.into());
        set(
            &identified,
            "nearest_by",
            JsValue::from_str(&found.nearest_by),
        );
        set(
            &identified,
            "periodicity",
            JsValue::from_f64(found.fingerprint.periodicity as f64),
        );
        set(&object, "identified", identified.into());
        array.push(&object.into());
    }
    Ok(array)
}

/// Adds hydrogens at the missing valences.
///
/// # Errors
///
/// Returns a JavaScript `Error` for invalid text or a bad distance.
#[wasm_bindgen]
pub fn hydrogens_info(text: &str, distance: f64) -> Result<JsValue, JsValue> {
    let structure = parse_any(text).map_err(to_js_error)?;
    let options = xtalforge::build::hydrogens::HydrogenOptions {
        distance,
        ..Default::default()
    };
    let filled =
        xtalforge::build::hydrogens::add_hydrogens(&structure, &options).map_err(to_js_error)?;
    summary_object(&filled)
}
