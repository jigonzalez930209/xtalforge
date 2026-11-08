//! Python bindings for `xtalforge`.
//!
//! The module layout is `xtalforge._xtalforge`, re-exported by the Python
//! package in `python/xtalforge`. The public surface is kept at parity with
//! the Rust, CLI and WASM APIs.

#![allow(unsafe_code)] // pyo3 macros expand to FFI code; this crate writes none.

use pyo3::create_exception;
use pyo3::exceptions::PyException;
use pyo3::prelude::*;

create_exception!(
    _xtalforge,
    XtalforgeError,
    PyException,
    "Base class for every xtalforge error."
);
create_exception!(
    _xtalforge,
    InputError,
    XtalforgeError,
    "A geometry or value does not satisfy the API preconditions."
);

/// Version of the workspace.
#[pyfunction]
fn version() -> &'static str {
    xtalforge::core::version()
}

/// Volume of a cell in cubic Angstrom.
#[pyfunction]
#[pyo3(signature = (a, b, c, alpha, beta, gamma))]
fn cell_volume(a: f64, b: f64, c: f64, alpha: f64, beta: f64, gamma: f64) -> PyResult<f64> {
    xtalforge::core::cell_volume(a, b, c, alpha, beta, gamma).map_err(|error| {
        let args = (error.code(), error.message());
        match error {
            xtalforge::core::Error::InvalidInput(_) => InputError::new_err(args),
            _ => XtalforgeError::new_err(args),
        }
    })
}

#[pymodule]
fn _xtalforge(module: &Bound<'_, PyModule>) -> PyResult<()> {
    module.add("__version__", xtalforge::core::version())?;
    module.add("XtalforgeError", module.py().get_type::<XtalforgeError>())?;
    module.add("InputError", module.py().get_type::<InputError>())?;
    module.add_function(wrap_pyfunction!(version, module)?)?;
    module.add_function(wrap_pyfunction!(cell_volume, module)?)?;
    Ok(())
}
