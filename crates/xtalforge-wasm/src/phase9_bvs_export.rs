//! Bond-valence and oxidation-state exports (Task 9.3.3).
//!
//! These live apart from the other local-environment exports because they
//! use the analysis crate and therefore belong to the `analysis` package in
//! the split WASM build (Task 15.3.1).

use wasm_bindgen::prelude::*;

use crate::common::to_js_error;
use crate::phase9_export::{json, load};

/// Bond-valence sums as JSON.
///
/// # Errors
///
/// Returns a JavaScript `Error` for invalid text or cutoff.
#[wasm_bindgen]
pub fn bond_valence_sums(text: &str, format_hint: &str, cutoff: f64) -> Result<String, JsValue> {
    let report = xtalforge::analysis::bond_valence_sums(
        &load(text, format_hint)?,
        &xtalforge::analysis::BvsOptions {
            cutoff,
            ..xtalforge::analysis::BvsOptions::default()
        },
    )
    .map_err(to_js_error)?;
    json(&report)
}

/// Oxidation states as JSON.
///
/// # Errors
///
/// Returns a JavaScript `Error` for invalid text or cutoff.
#[wasm_bindgen]
pub fn oxidation_states(text: &str, format_hint: &str, cutoff: f64) -> Result<String, JsValue> {
    let report = xtalforge::analysis::oxidation_states(
        &load(text, format_hint)?,
        &xtalforge::analysis::BvsOptions {
            cutoff,
            ..xtalforge::analysis::BvsOptions::default()
        },
    )
    .map_err(to_js_error)?;
    json(&report)
}
