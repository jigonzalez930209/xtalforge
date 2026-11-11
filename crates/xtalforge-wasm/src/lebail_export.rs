//! Le Bail intensity extraction for the WASM surface (functional audit
//! 2026-11-11: the operation existed only in the CLI).

use wasm_bindgen::prelude::*;

use crate::common::{parse_any, set, to_js_error};

fn numbers(values: &[f64]) -> JsValue {
    let array = js_sys::Array::new();
    for value in values {
        array.push(&JsValue::from_f64(*value));
    }
    array.into()
}

/// Le Bail extraction with the default options.
///
/// Returns `{ positions, intensities, intensity_esds, r_factor, rwp,
/// chi_squared }`, the same fields as `xtalforge --json powder lebail`.
///
/// # Errors
///
/// Returns a JavaScript `Error` for invalid structure or pattern text.
#[wasm_bindgen]
pub fn powder_lebail_info(text: &str, pattern_text: &str) -> Result<JsValue, JsValue> {
    let structure = parse_any(text).map_err(to_js_error)?;
    let observed = xtalforge::io::read_xy(pattern_text).map_err(to_js_error)?;
    let result = xtalforge::analysis::le_bail(
        &structure,
        &observed,
        &xtalforge::analysis::LebailOptions::default(),
    )
    .map_err(to_js_error)?;
    let object = js_sys::Object::new();
    set(&object, "positions", numbers(&result.positions));
    set(&object, "intensities", numbers(&result.intensities));
    set(&object, "intensity_esds", numbers(&result.intensity_esds));
    set(&object, "r_factor", JsValue::from_f64(result.r_factor));
    set(&object, "rwp", JsValue::from_f64(result.rwp));
    set(
        &object,
        "chi_squared",
        JsValue::from_f64(result.chi_squared),
    );
    Ok(object.into())
}
