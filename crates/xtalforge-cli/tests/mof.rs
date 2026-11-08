//! Integration tests for the MOF CLI (Tasks 6.11-6.17).

#![allow(
    clippy::unwrap_used,
    clippy::expect_used,
    clippy::panic,
    clippy::indexing_slicing
)]

use std::process::{Command, Output};

fn run(args: &[&str]) -> Output {
    Command::new(env!("CARGO_BIN_EXE_xtalforge"))
        .args(args)
        .output()
        .expect("the CLI binary must run")
}

#[test]
fn mof_catalogues_and_build() {
    let output = run(&["--json", "mof", "topologies"]);
    assert!(output.status.success());
    let value: serde_json::Value = serde_json::from_slice(&output.stdout).expect("JSON");
    assert!(
        value
            .as_array()
            .expect("array")
            .iter()
            .any(|item| item["name"] == "pcu")
    );

    let output = run(&["--json", "mof", "blocks"]);
    assert!(output.status.success());
    let value: serde_json::Value = serde_json::from_slice(&output.stdout).expect("JSON");
    assert!(
        value
            .as_array()
            .expect("array")
            .iter()
            .any(|item| item["name"] == "Zn4O")
    );

    let directory = std::env::temp_dir().join(format!("xtalforge-mof-{}", std::process::id()));
    std::fs::create_dir_all(&directory).unwrap();
    let out = directory.join("pcu.cif");
    let output = run(&[
        "--json",
        "mof",
        "build",
        "pcu",
        "Zn4O",
        "benzene",
        "--out",
        &out.to_string_lossy(),
    ]);
    assert!(output.status.success());
    assert!(out.exists());
    std::fs::remove_dir_all(directory).unwrap();
}

#[test]
fn unknown_names_error() {
    let output = run(&["mof", "build", "nope", "Zn4O", "benzene"]);
    assert!(!output.status.success());
    let text = String::from_utf8(output.stderr).expect("stderr");
    assert!(text.contains("unknown topology"));
}
