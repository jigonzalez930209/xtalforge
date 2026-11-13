//! Builder exports for the Phase 11 usability helpers.
//!
//! SMILES embedding, the fragment library, substitution and guest placement
//! belong to the `build` package of the split WASM build (Task 15.3.1).

use serde::Deserialize;
use wasm_bindgen::prelude::*;

use crate::common::to_js_error;
use crate::phase9_export::{document, json, load, options};

#[derive(Debug, Default, Deserialize)]
#[serde(default)]
struct EmbedOptions {
    seed: u64,
    iterations: usize,
}

fn embed(smiles: &str, options_json: &str) -> Result<String, JsValue> {
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

/// Parses and embeds a SMILES string; returns structure JSON.
#[wasm_bindgen]
pub fn molecule_from_smiles(smiles: &str, options_json: &str) -> Result<String, JsValue> {
    embed(smiles, options_json)
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
    let fragment = xtalforge::build::find_fragment(name)
        .map_err(to_js_error)?
        .ok_or_else(|| JsValue::from_str("unknown fragment"))?;
    embed(&fragment.smiles, options_json)
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
