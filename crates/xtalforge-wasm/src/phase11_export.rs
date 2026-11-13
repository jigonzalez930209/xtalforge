//! WebAssembly exports for the Phase 11 usability helpers.
//!
//! The builder half (SMILES embedding, fragments, substitution and guest
//! placement) lives in `phase11_build_export` and belongs to the `build`
//! package of the split WASM build (Task 15.3.1).

use wasm_bindgen::prelude::*;

/// Reads a structure from text; alias of `read_structure`.
#[wasm_bindgen]
pub fn read(text: &str, format_hint: &str) -> Result<JsValue, JsValue> {
    crate::read_structure(text, format_hint)
}

/// Writes one `console.warn` deprecation message.
#[wasm_bindgen]
pub fn deprecation_warning(name: &str, replacement: &str) {
    let message = format!(
        "warning[deprecated_{name}]: use {replacement} instead; see docs/guide/deprecation.md"
    );
    let Ok(console) = js_sys::Reflect::get(&js_sys::global(), &JsValue::from_str("console")) else {
        return;
    };
    if console.is_undefined() {
        return;
    }
    if let Ok(warn) = js_sys::Reflect::get(&console, &JsValue::from_str("warn"))
        && let Ok(warn) = warn.dyn_into::<js_sys::Function>()
    {
        let _ = warn.call1(&console, &JsValue::from_str(&message));
    }
}

/// Machine-readable category for a stable error code.
#[wasm_bindgen]
pub fn error_category(code: &str) -> String {
    match code {
        "invalid_input" | "invalid_argument" | "unknown_verb" | "unknown_document"
        | "no_document" | "io" | "typing" | "unknown_element" | "unknown_type" | "topology"
        | "unknown_calculator" | "unknown_topology" | "unknown_block" => "input".to_owned(),
        "limit_exceeded" => "limit".to_owned(),
        "quantum" => "engine".to_owned(),
        "unsupported" => "unsupported".to_owned(),
        _ => "internal".to_owned(),
    }
}
