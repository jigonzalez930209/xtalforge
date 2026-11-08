//! Extra carbon and polymer builder exports (Tasks 6.5, 6.7-6.9).

use wasm_bindgen::prelude::*;

use xtalforge::build::{
    PolymerOptions, ZtcOptions, build_polymer, build_schwarzite, build_ztc, commensurate_twist,
    face_spiral, nanoribbon,
};

use crate::common::{set, summary_object, to_js_error};

fn ztc_options(samples: usize, cell: f64, isovalue: f64) -> ZtcOptions {
    ZtcOptions {
        samples,
        cell,
        isovalue,
    }
}

/// Builds a C20 dodecahedron.
///
/// # Errors
///
/// Returns a JavaScript `Error` for a builder failure.
#[wasm_bindgen]
pub fn build_c20_info() -> Result<JsValue, JsValue> {
    let structure = face_spiral(20).map_err(to_js_error)?;
    summary_object(&structure)
}

/// Builds two twisted graphene layers.
///
/// # Errors
///
/// Returns a JavaScript `Error` for invalid repeats or angle.
#[wasm_bindgen]
pub fn twist_info(n: u32, m: u32, theta: f64, layers: u32) -> Result<JsValue, JsValue> {
    let structure = commensurate_twist(n, m, theta, layers).map_err(to_js_error)?;
    summary_object(&structure)
}

/// Builds a honeycomb nanoribbon with vacuum.
///
/// # Errors
///
/// Returns a JavaScript `Error` for invalid dimensions.
#[wasm_bindgen]
pub fn ribbon_info(width: u32, length: u32, vacuum: f64) -> Result<JsValue, JsValue> {
    let structure = nanoribbon(width, length, vacuum).map_err(to_js_error)?;
    summary_object(&structure)
}

/// Builds a gyroid ZTC carbon model.
///
/// # Errors
///
/// Returns a JavaScript `Error` for invalid options.
#[wasm_bindgen]
pub fn ztc_info(samples: usize, cell: f64, isovalue: f64) -> Result<JsValue, JsValue> {
    let structure = build_ztc(&ztc_options(samples, cell, isovalue)).map_err(to_js_error)?;
    summary_object(&structure)
}

/// Builds a Schwarz P carbon model.
///
/// # Errors
///
/// Returns a JavaScript `Error` for invalid options.
#[wasm_bindgen]
pub fn schwarzite_info(samples: usize, cell: f64, isovalue: f64) -> Result<JsValue, JsValue> {
    let structure = build_schwarzite(&ztc_options(samples, cell, isovalue)).map_err(to_js_error)?;
    summary_object(&structure)
}

/// Builds a polymer chain and returns its steps.
///
/// # Errors
///
/// Returns a JavaScript `Error` for an unsupported template.
#[wasm_bindgen]
pub fn polymer_info(smiles: &str, monomers: u32, density: f64) -> Result<JsValue, JsValue> {
    let build =
        build_polymer(smiles, monomers, &PolymerOptions { density }).map_err(to_js_error)?;
    let object = summary_object(&build.structure)?;
    let object: js_sys::Object = object.unchecked_into();
    let steps = js_sys::Array::new();
    for step in &build.steps {
        steps.push(&JsValue::from_str(step));
    }
    set(&object, "steps", steps.into());
    Ok(object.into())
}
