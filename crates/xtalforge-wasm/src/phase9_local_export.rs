//! WebAssembly exports for local environments (Task 9.3).

use wasm_bindgen::prelude::*;

use crate::common::to_js_error;
use crate::phase9_export::{json, load};

/// Periodic Voronoi cells as JSON.
///
/// # Errors
///
/// Returns a JavaScript `Error` for invalid text or cutoff.
#[wasm_bindgen]
pub fn voronoi(text: &str, format_hint: &str, cutoff: f64) -> Result<String, JsValue> {
    let report =
        xtalforge::core::voronoi(&load(text, format_hint)?, cutoff).map_err(to_js_error)?;
    json(&report)
}

/// Coordination shell and shape measures as JSON.
///
/// # Errors
///
/// Returns a JavaScript `Error` for invalid text or an out-of-range index.
#[wasm_bindgen]
pub fn coordination(
    text: &str,
    format_hint: &str,
    index: usize,
    strategy: &str,
    cutoff: f64,
    solid_angle: f64,
) -> Result<String, JsValue> {
    let strategy = match strategy {
        "cutoff" => xtalforge::core::CoordinationStrategy::Cutoff(cutoff),
        "voronoi" => xtalforge::core::CoordinationStrategy::Voronoi {
            cutoff,
            solid_angle_min: solid_angle,
        },
        "crystalnn" | "crystal_nn" => xtalforge::core::CoordinationStrategy::CrystalNn { cutoff },
        other => {
            return Err(JsValue::from_str(&format!(
                "unknown coordination strategy '{other}'"
            )));
        }
    };
    let shell = xtalforge::core::coordination(&load(text, format_hint)?, index, strategy)
        .map_err(to_js_error)?;
    json(&shell)
}

/// Total and partial RDF as JSON.
///
/// # Errors
///
/// Returns a JavaScript `Error` for invalid text or range.
#[wasm_bindgen]
pub fn radial_distribution(
    text: &str,
    format_hint: &str,
    r_max: f64,
    bin_width: f64,
    smoothing: usize,
) -> Result<String, JsValue> {
    let distribution = xtalforge::core::radial_distribution(
        &load(text, format_hint)?,
        &xtalforge::core::RdfOptions {
            r_max,
            bin_width,
            smoothing,
        },
    )
    .map_err(to_js_error)?;
    json(&distribution)
}

/// Structure fingerprint as JSON.
///
/// # Errors
///
/// Returns a JavaScript `Error` for invalid text or range.
#[wasm_bindgen]
pub fn structure_fingerprint(
    text: &str,
    format_hint: &str,
    r_max: f64,
    bin_width: f64,
) -> Result<String, JsValue> {
    let fingerprint = xtalforge::core::fingerprint(
        &load(text, format_hint)?,
        &xtalforge::core::RdfOptions {
            r_max,
            bin_width,
            smoothing: 1,
        },
    )
    .map_err(to_js_error)?;
    json(&fingerprint)
}

/// Similarity of two fingerprints.
///
/// # Errors
///
/// Returns a JavaScript `Error` for invalid text.
#[wasm_bindgen]
pub fn fingerprint_similarity(
    first: &str,
    second: &str,
    format_hint: &str,
    r_max: f64,
    bin_width: f64,
) -> Result<f64, JsValue> {
    let options = xtalforge::core::RdfOptions {
        r_max,
        bin_width,
        smoothing: 1,
    };
    let first =
        xtalforge::core::fingerprint(&load(first, format_hint)?, &options).map_err(to_js_error)?;
    let second =
        xtalforge::core::fingerprint(&load(second, format_hint)?, &options).map_err(to_js_error)?;
    Ok(first.similarity(&second))
}
