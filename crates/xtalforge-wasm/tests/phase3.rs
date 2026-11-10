//! WASM tests for projects, parameters, documents and workspaces (Phase 3).

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
    WasmDocument, param_schema, project_info, project_pack, workspace_archive,
    workspace_read_archive,
};

const SI_CIF: &str = include_str!("../../../tests/corpus/simple/Si.cif");

fn property(value: &JsValue, name: &str) -> JsValue {
    js_sys::Reflect::get(value, &JsValue::from_str(name)).expect("property must exist")
}

#[wasm_bindgen_test]
fn params_and_projects_export() {
    let schema = param_schema("symmetry").expect("schema");
    assert!(property(&schema, "properties").is_object());
    assert!(param_schema("unknown").is_err());

    let bytes = project_pack(SI_CIF).expect("pack");
    assert_eq!(&bytes[..2], b"PK");
    let info = project_info(&bytes).expect("info");
    assert_eq!(
        property(&info, "formula").as_string().as_deref(),
        Some("Si8")
    );
    assert_eq!(property(&info, "n_sites").as_f64(), Some(1.0));
    assert!(project_info(b"not a zip").is_err());
}

#[wasm_bindgen_test]
fn document_supports_undo_and_redo() {
    let mut document = WasmDocument::new(SI_CIF).expect("document");
    assert!(!document.can_undo());
    document.move_sites(vec![0], 0.01, 0.0, 0.0).expect("move");
    assert!(document.can_undo());
    assert!(document.dirty());
    assert_eq!(document.history().length(), 1);
    assert!(document.undo());
    assert!(document.can_redo());
    assert!(document.redo());
    assert!(!document.can_redo());
    let summary = document.summary().expect("summary");
    // The site moved off the 8a special position, so the whole cell holds its
    // 48-atom general orbit; the formula counts the cell, not the asymmetric unit.
    assert_eq!(
        property(&summary, "formula").as_string().as_deref(),
        Some("Si48")
    );
    assert!(document.to_cif().expect("cif").starts_with("data_"));
    assert!(document.set_element(99, "N").is_err());
}

#[wasm_bindgen_test]
fn workspace_archive_round_trip() {
    let records = r#"[{"id":"run-1","module":"demo","input":"Si.cif","status":"Done","outputs":[],"data":{},"created":1}]"#;
    let bytes = workspace_archive(records).expect("archive");
    assert_eq!(&bytes[..2], b"PK");
    let ids = workspace_read_archive(&bytes).expect("read");
    assert_eq!(ids.length(), 1);
    assert_eq!(ids.get(0).as_string().as_deref(), Some("run-1"));
    assert!(workspace_archive("{}").is_err());
}
