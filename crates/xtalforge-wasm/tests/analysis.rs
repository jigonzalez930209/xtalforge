//! WASM analysis and hydrogens tests (Tasks 5.11.2 and 6.6).

#![cfg(target_arch = "wasm32")]
#![allow(
    clippy::unwrap_used,
    clippy::expect_used,
    clippy::panic,
    clippy::indexing_slicing
)]

use wasm_bindgen::JsValue;
use wasm_bindgen_test::*;
use xtalforge_wasm::{hydrogens_info, net_invariants_info, porosity_info, zeo_parse_info};

const SI_CIF: &str = include_str!("../../../tests/corpus/simple/Si.cif");

fn property(value: &JsValue, name: &str) -> JsValue {
    js_sys::Reflect::get(value, &JsValue::from_str(name)).expect("property")
}

#[wasm_bindgen_test]
fn analysis_exports() {
    let metrics = porosity_info(SI_CIF, 0.5, 1.2).expect("porosity");
    assert!(property(&metrics, "cell_volume").as_f64().unwrap_or(0.0) > 0.0);
    let cgd =
        "CRYSTAL\nNAME pcu\nGROUP P1\nNODE A 6 0 0 0\nNODE B 6 0.5 0 0\nEDGE 0 0 0 0.5 0 0\nEND\n";
    let nets = net_invariants_info(cgd, 0.01).expect("nets");
    assert_eq!(nets.length(), 1);
    let identified = property(&nets.get(0), "identified");
    assert!(property(&identified, "names").is_object());
    let methane = "1\nmethane\nC 0 0 0\n";
    let filled = hydrogens_info(methane, 1.09).expect("hydrogens");
    assert_eq!(property(&filled, "n_sites").as_f64(), Some(5.0));
    assert!(porosity_info("not a structure", 0.5, 1.2).is_err());
    assert!(zeo_parse_info("nothing").is_err());
}
