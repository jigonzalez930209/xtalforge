//! WASM tests for the extra carbon and polymer builders.

#![cfg(target_arch = "wasm32")]
#![allow(
    clippy::unwrap_used,
    clippy::expect_used,
    clippy::panic,
    clippy::indexing_slicing
)]

use wasm_bindgen::JsValue;
use wasm_bindgen_test::*;
use xtalforge_wasm::{
    build_c20_info, polymer_info, ribbon_info, schwarzite_info, twist_info, ztc_info,
};

fn sites(value: &JsValue) -> f64 {
    js_sys::Reflect::get(value, &JsValue::from_str("n_sites"))
        .expect("n_sites")
        .as_f64()
        .expect("number")
}

#[wasm_bindgen_test]
fn extra_builders_export() {
    assert_eq!(sites(&build_c20_info().expect("c20")), 20.0);
    assert_eq!(sites(&twist_info(2, 2, 0.0, 2).expect("twist")), 16.0);
    assert!(sites(&ribbon_info(2, 3, 12.0).expect("ribbon")) > 0.0);
    assert!(sites(&ztc_info(8, 12.0, 0.0).expect("ztc")) > 0.0);
    assert!(sites(&schwarzite_info(8, 12.0, 0.0).expect("schwarzite")) > 0.0);
    let polymer = polymer_info("C=C", 3, 1.0).expect("polymer");
    assert!(sites(&polymer) > 0.0);
    assert!(ztc_info(1, 12.0, 0.0).is_err());
}
