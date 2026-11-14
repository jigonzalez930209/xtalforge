//! WASM ML, Parquet and ONNX tests (Phase 21).

#![cfg(all(target_arch = "wasm32", feature = "full"))]
#![allow(
    clippy::unwrap_used,
    clippy::expect_used,
    clippy::panic,
    clippy::indexing_slicing
)]

use wasm_bindgen::JsValue;
use wasm_bindgen_test::*;
use xtalforge_wasm::{
    crystal_graph_info, ml_parquet_read_info, ml_parquet_write_info, onnx_energy_info,
    onnx_predict_info, soap_descriptor_info,
};

const MODEL: &[u8] = include_bytes!("../../xtalforge-ml/data/onnx_linear.onnx");
const SI_CIF: &str = include_str!("../../../tests/corpus/simple/Si.cif");

fn property(value: &JsValue, name: &str) -> JsValue {
    js_sys::Reflect::get(value, &JsValue::from_str(name)).expect("property")
}

#[wasm_bindgen_test]
fn ml_exports_run() {
    let graph = crystal_graph_info(SI_CIF, 4.0, Some(4)).expect("graph");
    assert!(property(&graph, "edge_index").is_object());
    let soap = soap_descriptor_info(SI_CIF, 4.0, 2, 1, 0.5, true).expect("soap");
    assert_eq!(property(&soap, "n_features").as_f64(), Some(6.0));

    let dataset = xtalforge::analysis::dataset::Dataset::from_records(
        "wasm",
        vec![
            xtalforge::analysis::dataset::DatasetRecord::new(
                "a",
                xtalforge::io::read_structure(SI_CIF, "cif").expect("structure"),
            )
            .expect("record"),
        ],
    )
    .expect("dataset");
    let jsonl = dataset.to_json_lines().expect("jsonl");
    let bytes = ml_parquet_write_info(&jsonl).expect("parquet write");
    let restored = ml_parquet_read_info(&bytes).expect("parquet read");
    let text = property(&restored, "dataset").as_string().expect("dataset");
    assert!(text.contains("\"id\":\"a\""));
}

#[wasm_bindgen_test]
fn onnx_exports_run() {
    let values = onnx_predict_info(MODEL, "[0,0,0,2,0,0]", 6).expect("predict");
    let array = js_sys::Array::from(&values);
    let value = array.get(0).as_f64().expect("value");
    assert!((value - 0.2).abs() < 1e-5);
    let poscar = "C2\n1.0\n20 0 0\n0 20 0\n0 0 20\nC\n2\nDirect\n0 0 0\n0.1 0 0\n";
    let energy = onnx_energy_info(poscar, MODEL, 6, 0.001).expect("energy");
    assert!((energy - 0.2).abs() < 1e-5);
}
