//! WASM Pawley test, run with `wasm-pack test --node`.

#![cfg(target_arch = "wasm32")]
#![allow(
    clippy::unwrap_used,
    clippy::expect_used,
    clippy::panic,
    clippy::indexing_slicing
)]

use wasm_bindgen::JsCast;
use wasm_bindgen::JsValue;
use wasm_bindgen_test::*;
use xtalforge_wasm::powder_pawley_info;

fn property(value: &JsValue, name: &str) -> JsValue {
    js_sys::Reflect::get(value, &JsValue::from_str(name)).expect("property must exist")
}

fn as_array(value: &JsValue) -> js_sys::Array {
    value
        .clone()
        .dyn_into::<js_sys::Array>()
        .expect("the export must return a JavaScript array")
}

#[wasm_bindgen_test]
fn powder_pawley_recovers_the_golden_pattern() {
    let golden: serde_json::Value = serde_json::from_str(include_str!(
        "../../../tests/golden/phase-7/pawley/one-atom-cubic.json"
    ))
    .unwrap();
    let cell = golden["cell"].as_array().unwrap();
    let value = |index: usize| cell[index].as_f64().unwrap();
    let cif = format!(
        "data_one\n_cell_length_a {}\n_cell_length_b {}\n_cell_length_c {}\n\
         _cell_angle_alpha {}\n_cell_angle_beta {}\n_cell_angle_gamma {}\n\
         loop_\n_atom_site_label\n_atom_site_type_symbol\n_atom_site_fract_x\n\
         _atom_site_fract_y\n_atom_site_fract_z\n_atom_site_occupancy\nC1 C 0 0 0 1\n",
        value(0),
        value(1),
        value(2),
        value(3),
        value(4),
        value(5)
    );
    let mut pattern = String::new();
    for peak in golden["peaks"].as_array().unwrap() {
        pattern.push_str(&format!(
            "{} {}\n",
            peak[0].as_f64().unwrap(),
            peak[1].as_f64().unwrap()
        ));
    }
    let result = powder_pawley_info(&cif, &pattern, "gaussian").expect("Pawley fit");
    let a = as_array(&property(&result, "cell"))
        .get(0)
        .as_f64()
        .unwrap();
    assert!((a - value(0)).abs() < 5e-3);
    let zero = property(&result, "zero_shift").as_f64().unwrap();
    assert!((zero - golden["zero_shift"].as_f64().unwrap()).abs() < 5e-3);
    assert!(property(&result, "rwp").as_f64().unwrap() < 0.02);
    assert_eq!(as_array(&property(&result, "reflections")).length(), 80);
}

#[wasm_bindgen_test]
fn powder_refine_reports_the_mccusker_stages() {
    let cif = include_str!("../../../tests/corpus/simple/TiO2.cif");
    let pattern = include_str!("../../../tests/golden/phase-7/rietveld/rutile.xy");
    let result = xtalforge_wasm::powder_refine_info(
        cif, pattern, true, false, false, true, true, true, false, true,
    )
    .expect("Rietveld refinement");
    let stages = as_array(&property(&result, "stages"));
    assert_eq!(stages.length(), 4);
    assert_eq!(stages.get(0).as_string().unwrap(), "background_scale");
    let esds = property(&result, "parameter_esds");
    assert!(!property(&esds, "zero").is_undefined());
    let rwp = property(&result, "rwp").as_f64().unwrap();
    assert!(rwp < 0.3, "Rwp {rwp}");
}

#[wasm_bindgen_test]
fn powder_automatic_returns_the_leading_cell_table() {
    let cif = include_str!("../../../tests/corpus/simple/TiO2.cif");
    let pattern = include_str!("../../../tests/golden/phase-7/rietveld/rutile.xy");
    let result = xtalforge_wasm::powder_automatic_info(cif, pattern).expect("automatic run");
    let rows = as_array(&property(&result, "leading_cells"));
    assert!(rows.length() > 0);
    let first = rows.get(0);
    assert!(
        !property(&first, "class_symbol")
            .as_string()
            .unwrap()
            .is_empty()
    );
    assert!(property(&first, "class_rank").as_f64().unwrap() >= 1.0);
    let groups = as_array(&property(&first, "space_groups"));
    assert!(groups.length() > 0);
    let group = groups.get(0);
    assert!(property(&group, "number").as_f64().unwrap() > 0.0);
    assert!(!property(&group, "symbol").as_string().unwrap().is_empty());
}
