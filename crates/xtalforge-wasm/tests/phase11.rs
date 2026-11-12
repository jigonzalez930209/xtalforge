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
    deprecation_warning, fragment_library_info, molecule_from_smiles, read, sketch_roundtrip,
    substitute_info,
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
