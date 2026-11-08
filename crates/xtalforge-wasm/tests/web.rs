//! WASM integration tests, run with `wasm-pack test --node`.

#![cfg(target_arch = "wasm32")]
#![allow(
    clippy::unwrap_used,
    clippy::expect_used,
    clippy::panic,
    clippy::indexing_slicing
)]

use serde::Deserialize;
use wasm_bindgen_test::*;
use xtalforge_wasm::{cell_volume, version};

#[derive(Debug, Deserialize)]
struct Vectors {
    cell_volume: Vec<Case>,
}

#[derive(Debug, Deserialize)]
struct Case {
    a: f64,
    b: f64,
    c: f64,
    alpha: f64,
    beta: f64,
    gamma: f64,
    expected: f64,
    tolerance: f64,
}

#[wasm_bindgen_test]
fn version_is_not_empty() {
    let value = version();
    assert!(!value.is_empty());
    assert_eq!(value, xtalforge::core::version());
}

#[wasm_bindgen_test]
fn cell_volume_shared_vectors() {
    let vectors: Vectors =
        serde_json::from_str(include_str!("../../../tests/vectors/hello_world.json"))
            .expect("shared hello-world vectors must be valid JSON");
    for case in vectors.cell_volume {
        let got = cell_volume(case.a, case.b, case.c, case.alpha, case.beta, case.gamma)
            .expect("valid vector case");
        assert!((got - case.expected).abs() <= case.tolerance);
    }
}

#[wasm_bindgen_test]
fn invalid_cell_returns_an_error() {
    assert!(cell_volume(0.0, 1.0, 1.0, 90.0, 90.0, 90.0).is_err());
}
