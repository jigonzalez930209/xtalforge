//! WASM integration tests, run with `wasm-pack test --node`.

#![cfg(target_arch = "wasm32")]
#![allow(
    clippy::unwrap_used,
    clippy::expect_used,
    clippy::panic,
    clippy::indexing_slicing
)]

use serde::Deserialize;
use wasm_bindgen::JsValue;
use wasm_bindgen_test::*;
use xtalforge_wasm::{cell_volume, element_info, lattice_info, lattice_minimum_image, version};

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

#[derive(Debug, Deserialize)]
struct ElementVectors {
    elements: Vec<ElementCase>,
}

#[derive(Debug, Deserialize)]
struct ElementCase {
    input: String,
    symbol: String,
    z: u8,
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

fn property(value: &JsValue, name: &str) -> JsValue {
    js_sys::Reflect::get(value, &JsValue::from_str(name)).expect("property must exist")
}

#[wasm_bindgen_test]
fn element_info_shared_vectors() {
    let vectors: ElementVectors =
        serde_json::from_str(include_str!("../../../tests/vectors/elements.json"))
            .expect("element vectors must be valid JSON");
    for case in vectors.elements {
        let value = element_info(&case.input).expect("valid element input");
        assert_eq!(
            property(&value, "symbol").as_string().as_deref(),
            Some(case.symbol.as_str())
        );
        assert_eq!(property(&value, "z").as_f64(), Some(f64::from(case.z)));
    }
}

#[wasm_bindgen_test]
fn invalid_element_returns_an_error() {
    assert!(element_info("Zz").is_err());
}

#[wasm_bindgen_test]
fn lattice_shared_vectors() {
    let vectors: serde_json::Value =
        serde_json::from_str(include_str!("../../../tests/vectors/lattice.json"))
            .expect("lattice vectors must be valid JSON");
    for case in vectors["cases"].as_array().expect("cases array") {
        let number = |name: &str| case[name].as_f64().expect("numeric field");
        let info = lattice_info(
            number("a"),
            number("b"),
            number("c"),
            number("alpha"),
            number("beta"),
            number("gamma"),
        )
        .expect("valid cell");
        let volume = property(&info, "volume").as_f64().expect("volume");
        assert!((volume - number("volume")).abs() <= 1e-9);

        let delta: Vec<f64> = case["minimum_image"]["delta"]
            .as_array()
            .expect("delta")
            .iter()
            .map(|value| value.as_f64().expect("number"))
            .collect();
        let result = lattice_minimum_image(
            number("a"),
            number("b"),
            number("c"),
            number("alpha"),
            number("beta"),
            number("gamma"),
            delta,
        )
        .expect("valid cell");
        let distance = property(&result, "distance").as_f64().expect("distance");
        assert!(
            (distance
                - case["minimum_image"]["distance"]
                    .as_f64()
                    .expect("expected"))
            .abs()
                <= 1e-9
        );
    }
}

#[wasm_bindgen_test]
fn invalid_lattice_returns_an_error() {
    assert!(lattice_info(0.0, 1.0, 1.0, 90.0, 90.0, 90.0).is_err());
}
