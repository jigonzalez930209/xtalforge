//! WebAssembly bindings for `xtalforge`.
//!
//! The module exports the same hello-world surface as the Rust, CLI and
//! Python packages. The wasm-bindgen attributes are compiled only for
//! `wasm32-unknown-unknown`; a host build of this crate is intentionally
//! empty so `cargo test --workspace` stays green.

#[cfg(target_arch = "wasm32")]
use wasm_bindgen::prelude::*;

/// Installs the panic hook in debug builds only.
#[cfg(all(target_arch = "wasm32", debug_assertions))]
#[wasm_bindgen(start)]
pub fn install_panic_hook() {
    console_error_panic_hook::set_once();
}

/// Version of the workspace.
#[cfg(target_arch = "wasm32")]
#[wasm_bindgen]
pub fn version() -> String {
    xtalforge::core::version().to_owned()
}

/// Volume of a cell in cubic Angstrom.
///
/// # Errors
///
/// Returns a JavaScript `Error` whose `code` property is the stable error
/// code from the Rust core.
#[cfg(target_arch = "wasm32")]
#[wasm_bindgen]
pub fn cell_volume(
    a: f64,
    b: f64,
    c: f64,
    alpha: f64,
    beta: f64,
    gamma: f64,
) -> Result<f64, JsValue> {
    xtalforge::core::cell_volume(a, b, c, alpha, beta, gamma).map_err(to_js_error)
}

/// Unit-cell information from the six parameters.
///
/// # Errors
///
/// Returns a JavaScript `Error` with a stable `code` property for an invalid
/// cell.
#[cfg(target_arch = "wasm32")]
#[wasm_bindgen]
pub fn lattice_info(
    a: f64,
    b: f64,
    c: f64,
    alpha: f64,
    beta: f64,
    gamma: f64,
) -> Result<JsValue, JsValue> {
    let lattice =
        xtalforge::core::lattice_info(a, b, c, alpha, beta, gamma).map_err(to_js_error)?;
    let object = js_sys::Object::new();
    set(&object, "parameters", array6(lattice.parameters()).into());
    set(&object, "volume", JsValue::from_f64(lattice.volume()));
    let matrix = js_sys::Array::new();
    for row in lattice.matrix() {
        matrix.push(&array3(row).into());
    }
    set(&object, "matrix", matrix.into());
    set(
        &object,
        "right_handed",
        JsValue::from_bool(lattice.is_right_handed()),
    );
    Ok(object.into())
}

/// Minimum-image fractional vector and distance.
///
/// # Errors
///
/// Returns a JavaScript `Error` with a stable `code` property for an invalid
/// cell or a delta that is not three numbers.
#[cfg(target_arch = "wasm32")]
#[wasm_bindgen]
pub fn lattice_minimum_image(
    a: f64,
    b: f64,
    c: f64,
    alpha: f64,
    beta: f64,
    gamma: f64,
    delta: Vec<f64>,
) -> Result<JsValue, JsValue> {
    let [dx, dy, dz] = delta.as_slice() else {
        return Err(JsValue::from_str(
            "invalid_input: delta must contain three fractional coordinates",
        ));
    };
    let (vector, distance) =
        xtalforge::core::lattice_minimum_image(a, b, c, alpha, beta, gamma, [*dx, *dy, *dz])
            .map_err(to_js_error)?;
    let object = js_sys::Object::new();
    set(&object, "delta", array3(vector).into());
    set(&object, "distance", JsValue::from_f64(distance));
    Ok(object.into())
}

#[cfg(target_arch = "wasm32")]
fn array3(values: [f64; 3]) -> js_sys::Array {
    let array = js_sys::Array::new();
    for value in values {
        array.push(&JsValue::from_f64(value));
    }
    array
}

#[cfg(target_arch = "wasm32")]
fn array6(values: [f64; 6]) -> js_sys::Array {
    let array = js_sys::Array::new();
    for value in values {
        array.push(&JsValue::from_f64(value));
    }
    array
}

/// Soft and hard bounds for a profile and budget.
///
/// # Errors
///
/// Returns a JavaScript `Error` for an unknown budget.
#[cfg(target_arch = "wasm32")]
#[wasm_bindgen]
pub fn limits_bounds(profile: &str, budget: &str) -> Result<js_sys::Array, JsValue> {
    let (soft, hard) = xtalforge::core::limits_bounds(profile, budget).map_err(to_js_error)?;
    let array = js_sys::Array::new();
    array.push(&JsValue::from_f64(soft as f64));
    array.push(&hard.map_or(JsValue::NULL, |value| JsValue::from_f64(value as f64)));
    Ok(array)
}

/// Generic limit check for a profile, budget and estimate.
///
/// # Errors
///
/// Returns a JavaScript `Error` for an unknown budget.
#[cfg(target_arch = "wasm32")]
#[wasm_bindgen]
pub fn limits_check(profile: &str, budget: &str, estimate: f64) -> Result<JsValue, JsValue> {
    let estimate = if estimate.is_finite() && estimate > 0.0 {
        estimate.round() as u64
    } else {
        0
    };
    let verdict = xtalforge::core::limits_check(profile, budget, estimate).map_err(to_js_error)?;
    let object = js_sys::Object::new();
    set(&object, "level", JsValue::from_str(verdict.level.name()));
    set(
        &object,
        "estimate",
        JsValue::from_f64(verdict.estimate as f64),
    );
    set(&object, "sentence", JsValue::from_str(&verdict.sentence));
    set(&object, "refused", JsValue::from_bool(verdict.refused()));
    set(&object, "warned", JsValue::from_bool(verdict.warned()));
    Ok(object.into())
}

/// Summary of a structure JSON document.
///
/// # Errors
///
/// Returns a JavaScript `Error` with a stable `code` property for invalid
/// input.
#[cfg(target_arch = "wasm32")]
#[wasm_bindgen]
pub fn structure_info(text: &str) -> Result<JsValue, JsValue> {
    let summary = xtalforge::core::structure_info(text).map_err(to_js_error)?;
    let object = js_sys::Object::new();
    set(&object, "name", JsValue::from_str(&summary.name));
    set(&object, "formula", JsValue::from_str(&summary.formula));
    set(&object, "hall", JsValue::from_str(&summary.hall));
    set(
        &object,
        "n_sites",
        JsValue::from_f64(summary.n_sites as f64),
    );
    set(&object, "volume", JsValue::from_f64(summary.volume));
    set(&object, "density", JsValue::from_f64(summary.density));
    set(
        &object,
        "total_charge",
        JsValue::from_f64(summary.total_charge),
    );
    set(&object, "is_p1", JsValue::from_bool(summary.is_p1));
    Ok(object.into())
}

/// Reference data for an element by symbol or atomic number.
///
/// # Errors
///
/// Returns a JavaScript `Error` with a stable `code` property when the input
/// is not a known element.
#[cfg(target_arch = "wasm32")]
#[wasm_bindgen]
pub fn element_info(input: &str) -> Result<JsValue, JsValue> {
    let element = xtalforge::core::element_info(input).map_err(to_js_error)?;
    let object = js_sys::Object::new();
    set(&object, "z", JsValue::from_f64(f64::from(element.z)));
    set(&object, "symbol", JsValue::from_str(&element.symbol));
    set(&object, "name", JsValue::from_str(&element.name));
    set(&object, "mass", JsValue::from_f64(element.mass));
    set(
        &object,
        "covalent_radius",
        JsValue::from_f64(element.covalent_radius),
    );
    set(&object, "vdw_radius", JsValue::from_f64(element.vdw_radius));
    let color = js_sys::Array::new();
    for channel in element.color {
        color.push(&JsValue::from_f64(f64::from(channel)));
    }
    set(&object, "color", color.into());
    set(&object, "is_metal", JsValue::from_bool(element.is_metal));
    set(
        &object,
        "valence",
        element
            .valence
            .map_or(JsValue::NULL, |value| JsValue::from_f64(f64::from(value))),
    );
    set(
        &object,
        "electronegativity",
        element
            .electronegativity
            .map_or(JsValue::NULL, JsValue::from_f64),
    );
    set(
        &object,
        "pyykko_single",
        element
            .pyykko_single
            .map_or(JsValue::NULL, JsValue::from_f64),
    );
    set(
        &object,
        "pyykko_double",
        element
            .pyykko_double
            .map_or(JsValue::NULL, JsValue::from_f64),
    );
    set(
        &object,
        "pyykko_aromatic",
        element
            .pyykko_aromatic
            .map_or(JsValue::NULL, JsValue::from_f64),
    );
    Ok(object.into())
}

#[cfg(target_arch = "wasm32")]
fn set(object: &js_sys::Object, key: &str, value: JsValue) {
    let _ = js_sys::Reflect::set(object, &JsValue::from_str(key), &value);
}

#[cfg(target_arch = "wasm32")]
fn to_js_error(error: xtalforge::core::Error) -> JsValue {
    let js_error = js_sys::Error::new(&error.message());
    let _ = js_sys::Reflect::set(
        &js_error,
        &JsValue::from_str("code"),
        &JsValue::from_str(error.code()),
    );
    js_error.into()
}
