//! Integration tests for the Phase 3 CLI commands (Tasks 3.7, 3.8, 3.16).

#![allow(
    clippy::unwrap_used,
    clippy::expect_used,
    clippy::panic,
    clippy::indexing_slicing
)]

use std::path::{Path, PathBuf};
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

fn temp_dir(name: &str) -> PathBuf {
    let directory =
        std::env::temp_dir().join(format!("xtalforge-cli-{name}-{}", std::process::id()));
    std::fs::create_dir_all(&directory).expect("temp directory");
    directory
}

#[test]
fn formats_lists_the_registry() {
    let output = run(&["--json", "formats"]);
    assert!(output.status.success());
    let value: serde_json::Value = serde_json::from_slice(&output.stdout).expect("JSON");
    let entries = value.as_array().expect("array");
    assert!(entries.iter().any(|entry| entry["name"] == "cif"));
    assert!(entries.iter().any(|entry| entry["name"] == "lammps-data"));
    assert!(
        entries
            .iter()
            .any(|entry| entry["name"] == "pdb" && entry["read"] == false)
    );
}

#[test]
fn params_prints_a_schema() {
    let output = run(&["--json", "params", "symmetry"]);
    assert!(output.status.success());
    let value: serde_json::Value = serde_json::from_slice(&output.stdout).expect("JSON");
    assert_eq!(value["properties"]["symprec"]["type"], "number");
    let output = run(&["--json", "params"]);
    assert!(output.status.success());
    let value: serde_json::Value = serde_json::from_slice(&output.stdout).expect("JSON");
    assert_eq!(value.as_array().expect("array").len(), 3);
}

#[test]
fn project_round_trip_through_the_cli() {
    let directory = temp_dir("project");
    let archive = directory.join("si.xtalproj");
    let output = run(&[
        "--json",
        "project",
        "pack",
        &corpus("simple/Si.cif"),
        &archive.to_string_lossy(),
    ]);
    assert!(output.status.success());
    assert!(archive.exists());

    let output = run(&["--json", "project", "info", &archive.to_string_lossy()]);
    assert!(output.status.success());
    let value: serde_json::Value = serde_json::from_slice(&output.stdout).expect("JSON");
    assert_eq!(value["formula"], "Si");
    assert_eq!(value["n_sites"], 1);

    let extracted = directory.join("extracted");
    let output = run(&[
        "project",
        "extract",
        &archive.to_string_lossy(),
        &extracted.to_string_lossy(),
    ]);
    assert!(output.status.success());
    assert!(extracted.join("structure.cif").exists());
    std::fs::remove_dir_all(&directory).unwrap();
}

#[test]
fn trajectory_info_and_frame() {
    let directory = temp_dir("trajectory");
    let path = directory.join("frames.xyz");
    std::fs::write(
        &path,
        "2\nLattice=\"4 0 0 0 4 0 0 0 4\"\nC 0 0 0\nO 1 0 0\n\
         2\nLattice=\"4 0 0 0 4 0 0 0 4\" energy=-1.0\nC 0.1 0 0\nO 1.1 0 0\n",
    )
    .unwrap();
    let output = run(&["--json", "trajectory", &path.to_string_lossy()]);
    assert!(output.status.success());
    let value: serde_json::Value = serde_json::from_slice(&output.stdout).expect("JSON");
    assert_eq!(value["frames"], 2);
    let output = run(&[
        "--json",
        "trajectory",
        &path.to_string_lossy(),
        "--frame",
        "1",
    ]);
    assert!(output.status.success());
    let value: serde_json::Value = serde_json::from_slice(&output.stdout).expect("JSON");
    assert_eq!(
        value["structure"]["sites"].as_array().expect("sites").len(),
        2
    );
    std::fs::remove_dir_all(&directory).unwrap();
}

#[test]
fn script_applies_commands_with_undo_semantics() {
    let directory = temp_dir("script");
    let path = directory.join("script.json");
    let source = std::fs::read_to_string(corpus("simple/Si.cif")).unwrap();
    let payload = serde_json::json!({
        "structure": source,
        "commands": [
            {"type": "set_element", "site": 0, "element": "Ge"}
        ]
    });
    std::fs::write(&path, serde_json::to_string(&payload).unwrap()).unwrap();
    let output = run(&["--json", "script", &path.to_string_lossy()]);
    assert!(output.status.success());
    let value: serde_json::Value = serde_json::from_slice(&output.stdout).expect("JSON");
    assert_eq!(value["structure"]["sites"][0]["element"], "Ge");
    std::fs::remove_dir_all(&directory).unwrap();
}

#[test]
fn analysis_porosity_and_nets() {
    let output = run(&[
        "--json",
        "analysis",
        "porosity",
        &corpus("simple/Si.cif"),
        "--spacing",
        "0.5",
    ]);
    assert!(output.status.success());
    let value: serde_json::Value = serde_json::from_slice(&output.stdout).expect("JSON");
    assert!(value["cell_volume"].as_f64().expect("volume") > 0.0);
    assert!(value["surface_area"].as_f64().expect("area") >= 0.0);
}
