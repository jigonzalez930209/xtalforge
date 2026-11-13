//! WASM integration tests for the Phase 11 usability exports.

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
    Structure, deprecation_warning, fragment_library_info, mesh_buffers, molecule_from_smiles,
    pxrd_buffers, read, scene_buffers, sketch_roundtrip, substitute_info,
};

const SI_CIF: &str = include_str!("../../../tests/corpus/simple/Si.cif");

#[wasm_bindgen_test]
fn read_alias_and_molecule_export() {
    let summary = read(SI_CIF, "cif").expect("read");
    assert!(!summary.is_undefined());
    let molecule = molecule_from_smiles("c1ccccc1", "{}").expect("molecule");
    let molecule: Value = serde_json::from_str(&molecule).expect("json");
    assert_eq!(molecule["schema_version"], 1);
}

#[wasm_bindgen_test]
fn sketch_and_fragments_and_substitution() {
    assert!(!sketch_roundtrip("c1ccccc1").unwrap_or_default().is_empty());
    let fragments = fragment_library_info().expect("fragments");
    let fragments: Value = serde_json::from_str(&fragments).expect("json");
    assert!(fragments.as_array().is_some_and(|list| !list.is_empty()));
    let replaced = substitute_info(SI_CIF, "cif", 0, "Ge").expect("substitute");
    let replaced: Value = serde_json::from_str(&replaced).expect("json");
    assert!(
        replaced["structure"]["sites"]
            .as_array()
            .is_some_and(|list| list.iter().any(|site| site["element"] == "Ge"))
    );
    deprecation_warning("old_name", "new_name");
}

#[wasm_bindgen_test]
fn structure_handle_and_typed_buffers() {
    let handle = Structure::from_text(SI_CIF, "cif").expect("handle");
    assert_eq!(handle.n_atoms().unwrap_or(0), 8);
    assert_eq!(
        handle.positions().map(|array| array.length()).unwrap_or(0),
        24
    );
    assert_eq!(handle.cell().length(), 9);
    let supercell = handle.supercell(2, 1, 1).expect("supercell");
    assert_eq!(supercell.n_atoms().unwrap_or(0), 16);
    let buffers = scene_buffers(SI_CIF, "cif").expect("scene buffers");
    assert!(!buffers.is_undefined());
    let pattern = pxrd_buffers(SI_CIF, "cif").expect("pattern buffers");
    assert!(!pattern.is_undefined());
    let mesh = mesh_buffers(SI_CIF, "cif").expect("mesh buffers");
    assert!(!mesh.is_undefined());
}
