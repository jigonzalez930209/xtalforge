//! Integration tests for the Phase 10 MOF commands.

#![allow(
    clippy::unwrap_used,
    clippy::expect_used,
    clippy::panic,
    clippy::indexing_slicing
)]

use std::process::{Command, Output};

const SILICON: &str = concat!(
    env!("CARGO_MANIFEST_DIR"),
    "/../../tests/corpus/simple/Si.cif"
);

fn run(args: &[&str]) -> Output {
    Command::new(env!("CARGO_BIN_EXE_xtalforge"))
        .args(args)
        .output()
        .expect("the CLI binary must run")
}

fn json(output: Output) -> serde_json::Value {
    assert!(
        output.status.success(),
        "stderr: {}",
        String::from_utf8_lossy(&output.stderr)
    );
    serde_json::from_slice(&output.stdout).expect("stdout must be JSON")
}

#[test]
fn mof_check_reports_findings() {
    let value = json(run(&["--json", "mof", "check", SILICON]));
    assert!(value["findings"].is_array());
    assert!(value["counts"].is_object());
}

#[test]
fn mof_pores_returns_native_descriptors() {
    let value = json(run(&[
        "--json",
        "mof",
        "pores",
        SILICON,
        "--spacing",
        "1.5",
    ]));
    assert!(value["descriptors"]["pld"].is_number());
    assert!(value["descriptors"]["lcd"].is_number());
    assert!(value["descriptors"]["channel_dimension"].is_number());
}

#[test]
fn mof_deconstruct_returns_a_net() {
    let value = json(run(&["--json", "mof", "deconstruct", SILICON]));
    assert!(value["units"].is_array());
    assert!(value["net"]["nodes"].is_array());
}

#[test]
fn mof_adsorption_guests_are_listed() {
    let value = json(run(&["--json", "mof", "adsorption", "guests"]));
    let names: Vec<&str> = value
        .as_array()
        .expect("list")
        .iter()
        .filter_map(|entry| entry["name"].as_str())
        .collect();
    assert!(names.contains(&"CO2"));
    assert!(names.contains(&"H2O"));
}

#[test]
fn mof_widom_is_finite() {
    let value = json(run(&[
        "--json",
        "mof",
        "adsorption",
        "widom",
        SILICON,
        "--samples",
        "20",
        "--guest",
        "CO2",
    ]));
    assert!(value["henry"].as_f64().is_some_and(f64::is_finite));
}
