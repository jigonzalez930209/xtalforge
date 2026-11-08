//! Powder refinement exports (Tasks 7.1-7.11).

use wasm_bindgen::prelude::*;

use crate::common::{parse_any, set, to_js_error};

fn pattern(text: &str) -> Result<Vec<xtalforge::io::PatternPoint>, JsValue> {
    xtalforge::io::read_xy(text).map_err(to_js_error)
}

/// Runs the automatic peak/index/Le Bail/Rietveld pipeline.
///
/// # Errors
///
/// Returns a JavaScript `Error` for invalid text, pattern or options.
#[wasm_bindgen]
pub fn powder_automatic_info(text: &str, pattern_text: &str) -> Result<JsValue, JsValue> {
    let structure = parse_any(text).map_err(to_js_error)?;
    let observed = pattern(pattern_text)?;
    let result = xtalforge::analysis::run_automatic(
        &structure,
        &observed,
        &xtalforge::analysis::AutomaticOptions::default(),
    )
    .map_err(to_js_error)?;
    let steps = js_sys::Array::new();
    for step in &result.steps {
        steps.push(&JsValue::from_str(step));
    }
    let object = js_sys::Object::new();
    set(
        &object,
        "peaks",
        JsValue::from_f64(result.peaks.len() as f64),
    );
    set(&object, "steps", steps.into());
    set(
        &object,
        "cell",
        result
            .index
            .map_or(JsValue::NULL, |index| JsValue::from_f64(index.a)),
    );
    set(
        &object,
        "rwp",
        result
            .phase
            .map_or(JsValue::NULL, |phase| JsValue::from_f64(phase.rwp)),
    );
    Ok(object.into())
}

/// Refines positions with a whole-pattern least-squares step.
///
/// # Errors
///
/// Returns a JavaScript `Error` for invalid text or a failed fit.
#[wasm_bindgen]
pub fn powder_refine_info(
    text: &str,
    pattern_text: &str,
    refine_positions: bool,
) -> Result<JsValue, JsValue> {
    let structure = parse_any(text).map_err(to_js_error)?;
    let observed = pattern(pattern_text)?;
    let result = xtalforge::analysis::refine_phase(
        &structure,
        &observed,
        &xtalforge::analysis::RietveldOptions::default(),
        &xtalforge::analysis::RefineOptions {
            refine_positions,
            ..Default::default()
        },
    )
    .map_err(to_js_error)?;
    let object = js_sys::Object::new();
    set(&object, "rwp", JsValue::from_f64(result.rwp));
    set(&object, "converged", JsValue::from_bool(result.converged));
    set(
        &object,
        "iterations",
        JsValue::from_f64(result.iterations as f64),
    );
    Ok(object.into())
}

/// Returns the free cell parameters and extinction symbol.
///
/// # Errors
///
/// Returns a JavaScript `Error` for an invalid space-group number.
#[wasm_bindgen]
pub fn powder_cells_info(number: u32) -> Result<JsValue, JsValue> {
    let parameters = xtalforge::analysis::free_cell_parameters(number).map_err(to_js_error)?;
    let extinction = xtalforge::analysis::extinction_symbol(number).map_err(to_js_error)?;
    let array = js_sys::Array::new();
    for parameter in parameters {
        array.push(&JsValue::from_str(parameter));
    }
    let object = js_sys::Object::new();
    set(&object, "free_parameters", array.into());
    set(&object, "extinction", JsValue::from_str(extinction));
    Ok(object.into())
}

/// Returns the rietx bridge inventory.
#[wasm_bindgen]
pub fn rietx_inventory_info() -> Result<JsValue, JsValue> {
    let bridge = xtalforge::analysis::rietx_inventory();
    let object = js_sys::Object::new();
    set(&object, "name", JsValue::from_str(bridge.name));
    set(&object, "license", JsValue::from_str(bridge.license));
    set(&object, "source", JsValue::from_str(bridge.source));
    let capabilities = js_sys::Array::new();
    for capability in bridge.capabilities {
        capabilities.push(&JsValue::from_str(capability));
    }
    set(&object, "capabilities", capabilities.into());
    Ok(object.into())
}
