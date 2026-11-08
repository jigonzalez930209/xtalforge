//! Integration tests for the Phase 4 CLI commands.

#![allow(
    clippy::unwrap_used,
    clippy::expect_used,
    clippy::panic,
    clippy::indexing_slicing
)]

use std::path::Path;
use std::process::{Command, Output};

fn run(args: &[&str]) -> Output {
    Command::new(env!("CARGO_BIN_EXE_xtalforge"))
        .args(args)
        .output()
        .expect("the CLI binary must run")
}

fn corpus(name: &str) -> String {
    Path::new(env!("CARGO_MANIFEST_DIR"))
        .join("../../tests/corpus")
        .join(name)
        .to_string_lossy()
        .into_owned()
}

#[test]
fn ff_energy_and_types() {
    let path = corpus("simple/Si.cif");
    let output = run(&["--json", "ff", "energy", &path]);
    assert!(output.status.success());
    let value: serde_json::Value = serde_json::from_slice(&output.stdout).expect("JSON");
    assert_eq!(value["calculator"], "uff");
    assert!(value["energy"].as_f64().expect("energy").is_finite());

    let output = run(&["--json", "ff", "types", &path]);
    assert!(output.status.success());
    let value: serde_json::Value = serde_json::from_slice(&output.stdout).expect("JSON");
    assert_eq!(value[0]["type"], "Si3");
}

#[test]
fn ff_engines_lists_uff() {
    let output = run(&["--json", "ff", "engines"]);
    assert!(output.status.success());
    let value: serde_json::Value = serde_json::from_slice(&output.stdout).expect("JSON");
    assert!(
        value
            .as_array()
            .expect("array")
            .iter()
            .any(|engine| engine["name"] == "uff")
    );
}

#[test]
fn ff_charges_are_neutral() {
    let path = corpus("simple/NaCl.cif");
    let output = run(&["--json", "ff", "charges", &path, "--method", "qeq"]);
    assert!(output.status.success());
    let value: serde_json::Value = serde_json::from_slice(&output.stdout).expect("JSON");
    assert!(value["total"].as_f64().expect("total").abs() < 1e-6);
    assert!(value["per_atom"].as_array().expect("charges").len() >= 2);
}

#[test]
fn ff_optimize_writes_a_structure() {
    let path = corpus("simple/Si.cif");
    let output = run(&["--json", "ff", "optimize", &path, "--max-steps", "2"]);
    assert!(output.status.success());
    let value: serde_json::Value = serde_json::from_slice(&output.stdout).expect("JSON");
    assert!(value["energy"].as_f64().expect("energy").is_finite());
    // A structure already at a stationary point legitimately takes zero steps.
    assert!(value["steps"].is_u64());
}

#[test]
fn ff_vibrations_and_thermo() {
    let path = corpus("simple/Si.cif");
    let output = run(&["--json", "ff", "vibrations", &path]);
    assert!(output.status.success());
    let value: serde_json::Value = serde_json::from_slice(&output.stdout).expect("JSON");
    assert_eq!(
        value["frequencies"].as_array().expect("frequencies").len(),
        24
    );
    let output = run(&["--json", "ff", "thermo", &path, "--temperature", "300"]);
    assert!(output.status.success());
    let value: serde_json::Value = serde_json::from_slice(&output.stdout).expect("JSON");
    assert!(value["gibbs"].as_f64().expect("gibbs").is_finite());
}
