//! Integration tests for the facade hello-world surface.
#![allow(
    clippy::unwrap_used,
    clippy::expect_used,
    clippy::panic,
    clippy::indexing_slicing
)]

use xtalforge::prelude::{Error, cell_volume, element_info, version};

#[test]
fn facade_reexports_the_hello_world_api() {
    assert_eq!(version(), env!("CARGO_PKG_VERSION"));
    assert!((cell_volume(1.0, 1.0, 1.0, 90.0, 90.0, 90.0).unwrap() - 1.0).abs() < 1e-12);
    let error = Error::InvalidInput("bad cell".to_owned());
    assert_eq!(error.code(), "invalid_input");
}

#[test]
fn facade_reexports_element_info() {
    assert_eq!(element_info("Fe").unwrap().symbol, "Fe");
    assert_eq!(element_info("26").unwrap().z, 26);
}

#[test]
fn facade_modules_expose_the_core_api() {
    let _: fn() -> &'static str = xtalforge::core::version;
    let _: fn(f64, f64, f64, f64, f64, f64) -> Result<f64, xtalforge::core::Error> =
        xtalforge::core::cell_volume;
}
