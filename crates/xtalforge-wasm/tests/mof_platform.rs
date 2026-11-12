//! WASM integration tests for the MOF platform exports.

#![cfg(target_arch = "wasm32")]
#![allow(
    clippy::unwrap_used,
    clippy::expect_used,
    clippy::panic,
    clippy::indexing_slicing
)]

use serde_json::Value;
use wasm_bindgen_test::*;
use xtalforge_wasm::{mof_check, mof_deconstruct, mof_guests, mof_pores};

const SI_CIF: &str = include_str!("../../../tests/corpus/simple/Si.cif");

#[wasm_bindgen_test]
fn mof_check_and_pores_return_json() {
    let check = mof_check(SI_CIF, "cif").expect("check");
    let check: Value = serde_json::from_str(&check).expect("json");
    assert!(check["findings"].is_array());
    let pores = mof_pores(SI_CIF, "cif", "{\"spacing\":1.5}").expect("pores");
    let pores: Value = serde_json::from_str(&pores).expect("json");
    assert!(pores["descriptors"]["lcd"].is_number());
}

#[wasm_bindgen_test]
fn mof_guests_and_deconstruction_exports() {
    let guests = mof_guests().expect("guests");
    let guests: Value = serde_json::from_str(&guests).expect("json");
    assert!(guests.as_array().is_some_and(|list| list.len() >= 5));
    let units = mof_deconstruct(SI_CIF, "cif").expect("deconstruct");
    let units: Value = serde_json::from_str(&units).expect("json");
    assert!(units["units"].is_array());
}
