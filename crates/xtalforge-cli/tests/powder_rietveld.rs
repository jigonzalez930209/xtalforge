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

fn rietveld_rwp(extra: &[&str]) -> f64 {
    let root = repository();
    let structure = root.join("tests/corpus/simple/TiO2.cif");
    let pattern = root.join("tests/golden/phase-7/rietveld/rutile.xy");
    let mut command = Command::new(env!("CARGO_BIN_EXE_xtalforge"));
    command.args([
        "--json",
        "powder",
        "rietveld",
        &structure.to_string_lossy(),
        &pattern.to_string_lossy(),
        "--refine",
    ]);
    command.args(extra);
    let output = command.output().expect("the CLI must run");
    assert!(output.status.success());
    let value: serde_json::Value = serde_json::from_slice(&output.stdout).expect("JSON");
    value["rwp"].as_f64().unwrap()
}

#[test]
fn no_doublet_models_one_wavelength() {
    // The golden pattern was generated with the Cu K alpha doublet, so dropping
    // it from the model must fit worse than the default.
    let with_doublet = rietveld_rwp(&[]);
    let single = rietveld_rwp(&["--no-doublet"]);
    assert!(single > with_doublet, "{single} vs {with_doublet}");
}

#[test]
fn no_doublet_requires_refine() {
    let root = repository();
    let output = Command::new(env!("CARGO_BIN_EXE_xtalforge"))
        .args([
            "powder",
            "rietveld",
            &root.join("tests/corpus/simple/TiO2.cif").to_string_lossy(),
            &root
                .join("tests/golden/phase-7/rietveld/rutile.xy")
                .to_string_lossy(),
            "--no-doublet",
        ])
        .output()
        .expect("the CLI must run");
    assert!(!output.status.success());
}
