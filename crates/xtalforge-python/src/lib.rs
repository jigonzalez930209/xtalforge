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

fn map_error(error: xtalforge::core::Error) -> PyErr {
    let args = (error.code(), error.message());
    match error {
        xtalforge::core::Error::InvalidInput(_) => InputError::new_err(args),
        _ => XtalforgeError::new_err(args),
    }
}

/// Volume of a cell in cubic Angstrom.
#[pyfunction]
#[pyo3(signature = (a, b, c, alpha, beta, gamma))]
fn cell_volume(a: f64, b: f64, c: f64, alpha: f64, beta: f64, gamma: f64) -> PyResult<f64> {
    xtalforge::core::cell_volume(a, b, c, alpha, beta, gamma).map_err(map_error)
}

/// Reference data for an element by symbol or atomic number.
#[pyfunction]
fn element_info<'py>(
    py: Python<'py>,
    input: &Bound<'_, PyAny>,
) -> PyResult<Bound<'py, pyo3::types::PyDict>> {
    let text = if let Ok(z) = input.extract::<u8>() {
        z.to_string()
    } else {
        input.extract::<String>()?
    };
    let element = xtalforge::core::element_info(&text).map_err(map_error)?;
    let dict = pyo3::types::PyDict::new(py);
    dict.set_item("z", element.z)?;
    dict.set_item("symbol", element.symbol)?;
    dict.set_item("name", element.name)?;
    dict.set_item("mass", element.mass)?;
    dict.set_item("covalent_radius", element.covalent_radius)?;
    dict.set_item("vdw_radius", element.vdw_radius)?;
    dict.set_item("color", element.color.to_vec())?;
    dict.set_item("is_metal", element.is_metal)?;
    dict.set_item("valence", element.valence)?;
    dict.set_item("electronegativity", element.electronegativity)?;
    dict.set_item("pyykko_single", element.pyykko_single)?;
    dict.set_item("pyykko_double", element.pyykko_double)?;
    dict.set_item("pyykko_aromatic", element.pyykko_aromatic)?;
    Ok(dict)
}

#[pymodule]
fn _xtalforge(module: &Bound<'_, PyModule>) -> PyResult<()> {
    module.add("__version__", xtalforge::core::version())?;
    module.add("XtalforgeError", module.py().get_type::<XtalforgeError>())?;
    module.add("InputError", module.py().get_type::<InputError>())?;
    module.add_function(wrap_pyfunction!(version, module)?)?;
    module.add_function(wrap_pyfunction!(cell_volume, module)?)?;
    module.add_function(wrap_pyfunction!(element_info, module)?)?;
    Ok(())
}
