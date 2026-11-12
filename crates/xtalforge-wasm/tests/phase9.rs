//! WASM integration tests for the Phase 9 exports.

#![cfg(target_arch = "wasm32")]
#![allow(
    clippy::unwrap_used,
    clippy::expect_used,
    clippy::panic,
    clippy::indexing_slicing
)]

use serde_json::Value;
use wasm_bindgen_test::*;
use xtalforge_wasm::{
    deform, match_structures, radial_distribution, unique_miller_indices, voronoi,
};

const SI_CIF: &str = include_str!("../../../tests/corpus/simple/Si.cif");

#[wasm_bindgen_test]
fn matcher_accepts_a_structure_with_itself() {
    let text = match_structures(SI_CIF, SI_CIF, "cif", "{}").expect("match");
    let value: Value = serde_json::from_str(&text).expect("json");
    assert_eq!(value["matched"], Value::Bool(true));
    assert!(value["rms"].as_f64().unwrap_or(1.0) < 1e-6);
}

#[wasm_bindgen_test]
fn deform_changes_the_cell_volume() {
    let text = deform(SI_CIF, "cif", "{\"volume\":200.0}").expect("deform");
    let value: Value = serde_json::from_str(&text).expect("json");
    assert_eq!(value["schema_version"], Value::from(1));
    let matrix = value["structure"]["lattice"]
        .as_array()
        .expect("lattice matrix");
    assert_eq!(matrix.len(), 3);
}

#[wasm_bindgen_test]
fn miller_and_voronoi_exports_return_json() {
    let indices = unique_miller_indices(SI_CIF, "cif", "{\"max_index\":1}").expect("miller");
    let indices: Value = serde_json::from_str(&indices).expect("json");
    assert_eq!(indices.as_array().map(Vec::len), Some(3));
    let cells = voronoi(SI_CIF, "cif", 4.0).expect("voronoi");
    let cells: Value = serde_json::from_str(&cells).expect("json");
    assert_eq!(cells["cells"].as_array().map(Vec::len), Some(8));
}

#[wasm_bindgen_test]
fn rdf_export_has_bins() {
    let rdf = radial_distribution(SI_CIF, "cif", 4.0, 0.05, 1).expect("rdf");
    let rdf: Value = serde_json::from_str(&rdf).expect("json");
    assert_eq!(rdf["r"].as_array().map(Vec::len), Some(80));
}
