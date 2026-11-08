//! WASM MOF tests (Tasks 6.11-6.17).

#![cfg(target_arch = "wasm32")]
#![allow(
    clippy::unwrap_used,
    clippy::expect_used,
    clippy::panic,
    clippy::indexing_slicing
)]

use wasm_bindgen::JsValue;
use wasm_bindgen_test::*;
use xtalforge_wasm::{build_mof_info, mof_blocks_info, mof_topologies_info};

fn property(value: &JsValue, name: &str) -> JsValue {
    js_sys::Reflect::get(value, &JsValue::from_str(name)).expect("property")
}

#[wasm_bindgen_test]
fn mof_exports() {
    let topologies = mof_topologies_info().expect("topologies");
    assert!(topologies.length() > 0);
    let blocks = mof_blocks_info().expect("blocks");
    assert!(blocks.length() > 0);
    let built = build_mof_info("pcu", "Zn4O", "benzene").expect("build");
    assert!(property(&built, "n_sites").as_f64().unwrap_or(0.0) > 0.0);
    assert!(build_mof_info("nope", "Zn4O", "benzene").is_err());
}
