//! WebAssembly exports for elasticity and thermodynamics (Task 9.6).

use serde::Deserialize;
use wasm_bindgen::prelude::*;

use crate::phase9_export::{json, load, options};

/// Elastic tensor from UFF stress finite differences as JSON.
///
/// # Errors
///
/// Returns a JavaScript `Error` for invalid text or strain.
#[wasm_bindgen]
pub fn elastic_tensor(text: &str, format_hint: &str, strain: f64) -> Result<String, JsValue> {
    let calculator = xtalforge::ff::UffCalculator::new(xtalforge::ff::CalculatorOptions::default())
        .map_err(|error| JsValue::from_str(&error.to_string()))?;
    let tensor = xtalforge::ff::elastic_tensor(
        &load(text, format_hint)?,
        &calculator,
        &xtalforge::ff::ElasticOptions {
            strain: if strain > 0.0 { strain } else { 0.005 },
        },
    )
    .map_err(crate::common::to_js_error)?;
    json(&tensor)
}

#[derive(Debug, Default, Deserialize)]
#[serde(default)]
pub(crate) struct EosOptions {
    pub points: Vec<(f64, f64)>,
    pub model: String,
}

/// Equation-of-state fit as JSON.
///
/// # Errors
///
/// Returns a JavaScript `Error` for an unknown model or too few points.
#[wasm_bindgen]
pub fn eos_fit(options_json: &str) -> Result<String, JsValue> {
    let options: EosOptions = options(options_json)?;
    let kind = match options.model.as_str() {
        "vinet" => xtalforge::core::EosKind::Vinet,
        "murnaghan" => xtalforge::core::EosKind::Murnaghan,
        _ => xtalforge::core::EosKind::BirchMurnaghan,
    };
    let points: Vec<xtalforge::core::EosPoint> = options
        .points
        .into_iter()
        .map(|(volume, energy)| xtalforge::core::EosPoint { volume, energy })
        .collect();
    let fit = xtalforge::core::fit_eos(&points, kind).map_err(crate::common::to_js_error)?;
    json(&fit)
}

/// Convex hull as JSON from a JSON list of phases.
///
/// # Errors
///
/// Returns a JavaScript `Error` for an invalid phase list.
#[wasm_bindgen]
pub fn convex_hull(entries_json: &str) -> Result<String, JsValue> {
    let entries: Vec<xtalforge::core::PhaseEntry> = serde_json::from_str(entries_json)
        .map_err(|error| JsValue::from_str(&format!("invalid phase entries: {error}")))?;
    let hull = xtalforge::core::convex_hull(&entries, None).map_err(crate::common::to_js_error)?;
    json(&hull)
}
