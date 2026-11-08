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
