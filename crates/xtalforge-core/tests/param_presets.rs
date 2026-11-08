//! Integration tests for the built-in parameter sets (Task 3.8.2).

#![allow(
    clippy::unwrap_used,
    clippy::expect_used,
    clippy::panic,
    clippy::indexing_slicing
)]

use serde_json::json;
use xtalforge_core::{PRESET_NAMES, preset, presets};

#[test]
fn every_preset_exists_and_validates_its_defaults() {
    assert_eq!(presets().unwrap().len(), PRESET_NAMES.len());
    for name in PRESET_NAMES {
        let set = preset(name).unwrap();
        assert!(set.check_required().is_ok());
        let schema = set.schema();
        assert_eq!(schema["type"], "object");
        assert!(schema["properties"].is_object());
    }
    assert!(preset("unknown").is_err());
}

#[test]
fn presets_validate_user_values() {
    let mut symmetry = preset("symmetry").unwrap();
    assert!(symmetry.set("symprec", json!(0.01)).is_ok());
    assert!(symmetry.set("symprec", json!(-1.0)).is_err());
    assert!(symmetry.set("hall_number", json!(531)).is_err());
    let mut pxrd = preset("pxrd").unwrap();
    assert!(pxrd.set("radiation", json!("mo-ka")).is_ok());
    assert!(pxrd.set("radiation", json!("gamma")).is_err());
    let mut uff = preset("uff").unwrap();
    assert!(uff.set("vdw", json!(false)).is_ok());
    assert!(uff.set("vdw", json!("no")).is_err());
}
