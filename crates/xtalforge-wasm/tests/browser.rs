//! Browser tests, run with `wasm-pack test --headless --chrome`.
//!
//! The other suites in this directory run in Node; this one is configured for
//! a real browser so the release surface is exercised there too.

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
    bonds, build_graphene_info, build_nanotube_info, cell_volume, read_cif, rings, scene_data,
    structure_info, version,
};

wasm_bindgen_test_configure!(run_in_browser);

const NACL_CIF: &str = include_str!("../../../tests/corpus/simple/NaCl.cif");

fn property(value: &JsValue, name: &str) -> JsValue {
    js_sys::Reflect::get(value, &JsValue::from_str(name)).expect("property must exist")
}

#[wasm_bindgen_test]
fn the_browser_reports_the_workspace_version() {
    assert_eq!(version(), xtalforge::core::version());
}

#[wasm_bindgen_test]
fn cell_volume_matches_the_documented_example() {
    let volume = cell_volume(3.0, 4.0, 5.0, 60.0, 70.0, 80.0).unwrap();
    assert!((volume - 48.827_604_320_084_97).abs() < 1e-9);
    assert!(cell_volume(-1.0, 1.0, 1.0, 90.0, 90.0, 90.0).is_err());
}

#[wasm_bindgen_test]
fn cif_text_reads_bonds_and_rings_in_the_browser() {
    let summary = read_cif(NACL_CIF).unwrap();
    assert_eq!(
        property(&summary, "formula").as_string().as_deref(),
        Some("Na4 Cl4")
    );
    let list = bonds(NACL_CIF).unwrap();
    assert_eq!(js_sys::Array::from(&list).length(), 24);
    assert!(js_sys::Array::from(&rings(NACL_CIF, 8).unwrap()).length() > 0);
}

#[wasm_bindgen_test]
fn scene_data_returns_the_full_geometry() {
    let scene = scene_data(NACL_CIF, "ball-stick", 1.0).unwrap();
    let atoms = js_sys::Array::from(&property(&scene, "atoms"));
    let bonds = js_sys::Array::from(&property(&scene, "bonds"));
    let edges = js_sys::Array::from(&property(&scene, "cell_edges"));
    assert_eq!(atoms.length(), 8);
    assert_eq!(bonds.length(), 24);
    assert_eq!(edges.length(), 12);
    let atom = atoms.get(0);
    assert!(property(&atom, "element").as_string().is_some());
    assert_eq!(
        js_sys::Array::from(&property(&atom, "position")).length(),
        3
    );
    assert!(property(&atom, "radius").as_f64().unwrap_or(0.0) > 0.0);
}

#[wasm_bindgen_test]
fn structure_info_accepts_a_versioned_envelope() {
    let structure = xtalforge::io::read_cif(NACL_CIF).unwrap();
    let wrapped = structure.to_document().to_json().unwrap();
    let summary = structure_info(&wrapped).unwrap();
    assert_eq!(property(&summary, "n_sites").as_f64(), Some(2.0));
}

#[wasm_bindgen_test]
fn carbon_builders_have_the_expected_atom_counts() {
    let graphene = build_graphene_info(3, 2, 1, 3.35, 15.0).unwrap();
    assert_eq!(property(&graphene, "n_sites").as_f64(), Some(12.0));
    let tube = build_nanotube_info(5, 5, 1, 10.0).unwrap();
    assert_eq!(property(&tube, "n_sites").as_f64(), Some(20.0));
}
