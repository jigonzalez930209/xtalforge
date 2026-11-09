//! Extra carbon and polymer builder exports (Tasks 6.5, 6.7-6.9).

use wasm_bindgen::prelude::*;

use xtalforge::build::{
    CarbonRecipe, PolymerOptions, ZtcOptions, build_carbon, build_polymer, build_schwarzite,
    build_ztc, commensurate_twist, face_spiral, nanoribbon,
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

/// Builds a disordered carbon from a recipe.
///
/// # Errors
///
/// Returns a JavaScript `Error` for an unsupported net, a bad repeat or a
/// surface that cannot be meshed.
#[wasm_bindgen]
#[allow(clippy::too_many_arguments)] // one argument per recipe field
pub fn build_carbon_recipe_info(
    net: &str,
    repeat_x: u32,
    repeat_y: u32,
    repeat_z: u32,
    density: f64,
    radius_ratio: f64,
    coverage: f64,
    seed: u32,
    hydrogen: f64,
    fluorine: f64,
    oxygen: f64,
) -> Result<JsValue, JsValue> {
    let recipe = CarbonRecipe {
        net: net.to_owned(),
        repeat: [repeat_x as usize, repeat_y as usize, repeat_z as usize],
        density,
        radius_ratio,
        coverage,
        seed: u64::from(seed),
        hydrogen,
        fluorine,
        oxygen,
        ..CarbonRecipe::default()
    };
    let build = build_carbon(&recipe).map_err(to_js_error)?;
    let object = summary_object(&build.structure)?;
    let object: js_sys::Object = object.unchecked_into();
    set(&object, "density", JsValue::from_f64(build.density));
    set(&object, "mean_bond", JsValue::from_f64(build.mean_bond));
    set(
        &object,
        "edge_fraction",
        JsValue::from_f64(build.edge_fraction),
    );
    set(&object, "pieces", JsValue::from_f64(build.pieces as f64));
    set(
        &object,
        "periodicity",
        JsValue::from_f64(build.periodicity as f64),
    );
    let census = js_sys::Array::new();
    for (size, count) in &build.census {
        let row = js_sys::Array::new();
        row.push(&JsValue::from_f64(*size as f64));
        row.push(&JsValue::from_f64(*count as f64));
        census.push(&row);
    }
    set(&object, "census", census.into());
    let notes = js_sys::Array::new();
    for note in &build.notes {
        notes.push(&JsValue::from_str(note));
    }
    set(&object, "notes", notes.into());
    Ok(object.into())
}
