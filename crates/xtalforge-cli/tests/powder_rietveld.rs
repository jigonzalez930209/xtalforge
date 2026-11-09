//! CLI test for `powder rietveld --refine` (Task 7.8).

#![allow(
    clippy::unwrap_used,
    clippy::expect_used,
    clippy::panic,
    clippy::indexing_slicing
)]

use std::path::PathBuf;
use std::process::Command;

fn repository() -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR"))
        .parent()
        .and_then(|path| path.parent())
        .expect("the crate lives under crates/")
        .to_path_buf()
}

#[test]
fn powder_rietveld_reports_stages_and_esds() {
    let root = repository();
    let structure = root.join("tests/corpus/simple/TiO2.cif");
    let pattern = root.join("tests/golden/phase-7/rietveld/rutile.xy");
    let output = Command::new(env!("CARGO_BIN_EXE_xtalforge"))
        .args([
            "--json",
            "powder",
            "rietveld",
            &structure.to_string_lossy(),
            &pattern.to_string_lossy(),
            "--refine",
        ])
        .output()
        .expect("the CLI must run");
    assert!(
        output.status.success(),
        "{}",
        String::from_utf8_lossy(&output.stderr)
    );
    let value: serde_json::Value = serde_json::from_slice(&output.stdout).expect("JSON");
    assert_eq!(
        value["stages"],
        serde_json::json!(["background_scale", "zero_cell", "width", "structure"])
    );
    assert!(value["rwp"].as_f64().unwrap() < 0.3);
    let parameters = value["parameters"].as_array().unwrap();
    let esds = value["parameter_esds"].as_array().unwrap();
    assert!(!parameters.is_empty());
    assert!(esds.iter().any(|entry| entry[0] == "zero"));
}
