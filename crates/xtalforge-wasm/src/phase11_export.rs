//! WebAssembly exports for the Phase 11 usability helpers.

use serde::Deserialize;
use wasm_bindgen::prelude::*;

use crate::common::to_js_error;
use crate::phase9_export::{document, json, load, options};

/// Reads a structure from text; alias of `read_structure`.
#[wasm_bindgen]
pub fn read(text: &str, format_hint: &str) -> Result<JsValue, JsValue> {
    crate::read_structure(text, format_hint)
}

#[derive(Debug, Default, Deserialize)]
#[serde(default)]
struct EmbedOptions {
    seed: u64,
    iterations: usize,
}

/// Parses and embeds a SMILES string; returns structure JSON.
#[wasm_bindgen]
pub fn molecule_from_smiles(smiles: &str, options_json: &str) -> Result<String, JsValue> {
    let options: EmbedOptions = options(options_json)?;
    let graph = xtalforge::build::parse_smiles(smiles).map_err(to_js_error)?;
    let structure = xtalforge::build::embed_smiles(
        &graph,
        &xtalforge::build::EmbedOptions {
            seed: options.seed,
            iterations: if options.iterations > 0 {
                options.iterations
            } else {
                1500
            },
            ..xtalforge::build::EmbedOptions::default()
        },
    )
    .map_err(to_js_error)?;
    document(&structure)
}

/// Round-trips a SMILES through the sketch model.
#[wasm_bindgen]
pub fn sketch_roundtrip(smiles: &str) -> Result<String, JsValue> {
    let sketch = xtalforge::build::from_smiles(smiles).map_err(to_js_error)?;
    sketch.to_smiles().map_err(to_js_error)
}

/// The embedded fragment library as JSON.
#[wasm_bindgen]
pub fn fragment_library_info() -> Result<String, JsValue> {
    json(&xtalforge::build::fragments().map_err(to_js_error)?)
}

/// Embeds one fragment by name; returns structure JSON.
#[wasm_bindgen]
pub fn fragment_info(name: &str, options_json: &str) -> Result<String, JsValue> {
    let options: EmbedOptions = options(options_json)?;
    let fragment = xtalforge::build::find_fragment(name)
        .map_err(to_js_error)?
        .ok_or_else(|| JsValue::from_str("unknown fragment"))?;
    let graph = xtalforge::build::parse_smiles(&fragment.smiles).map_err(to_js_error)?;
    let structure = xtalforge::build::embed_smiles(
        &graph,
        &xtalforge::build::EmbedOptions {
            seed: options.seed,
            iterations: if options.iterations > 0 {
                options.iterations
            } else {
                1500
            },
            ..xtalforge::build::EmbedOptions::default()
        },
    )
    .map_err(to_js_error)?;
    document(&structure)
}

/// Substitutes one site after a clearance check; returns structure JSON.
#[wasm_bindgen]
pub fn substitute_info(
    text: &str,
    format_hint: &str,
    index: usize,
    element: &str,
) -> Result<String, JsValue> {
    let structure = xtalforge::build::substitute(&load(text, format_hint)?, index, element)
        .map_err(to_js_error)?;
    document(&structure)
}

/// Places a guest after a clearance check; returns structure JSON.
#[wasm_bindgen]
pub fn place_guest_info(
    host: &str,
    format_hint: &str,
    guest_json: &str,
    offset: f64,
    clearance: f64,
) -> Result<String, JsValue> {
    let guest: Vec<(String, [f64; 3])> = serde_json::from_str(guest_json)
        .map_err(|error| JsValue::from_str(&format!("invalid guest JSON: {error}")))?;
    let structure = xtalforge::build::place_guest(
        &load(host, format_hint)?,
        &guest,
        [offset, offset, offset],
        if clearance > 0.0 { clearance } else { 1.0 },
    )
    .map_err(to_js_error)?;
    document(&structure)
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
