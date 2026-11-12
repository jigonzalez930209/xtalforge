//! Integration tests for the Phase 11 usability fixes.

#![allow(
    clippy::unwrap_used,
    clippy::expect_used,
    clippy::panic,
    clippy::indexing_slicing
)]

use std::io::Write;
use std::process::{Command, Output, Stdio};

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

#[test]
fn help_lists_the_registered_formats_and_builders() {
    let help = String::from_utf8(run(&["--help"]).stdout).expect("UTF-8");
    for format in [
        "cif",
        "poscar",
        "xyz",
        "pmg-json",
        "cssr",
        "gen",
        "lammps-data",
        "shelx",
        "pdb",
    ] {
        assert!(help.contains(format), "missing format {format}");
    }
    for builder in [
        "graphene",
        "nanotube",
        "C20",
        "C60",
        "twist",
        "ribbon",
        "ZTC",
        "schwarzite",
        "carbon",
        "polymer",
    ] {
        assert!(
            help.to_ascii_lowercase()
                .contains(&builder.to_ascii_lowercase()),
            "missing builder {builder}"
        );
    }
}

#[test]
fn doctor_reports_capabilities_and_diagnostics() {
    let output = run(&["--json", "doctor", SILICON]);
    assert!(output.status.success());
    let value: serde_json::Value = serde_json::from_slice(&output.stdout).expect("JSON");
    assert!(value["version"].is_string());
    assert!(value["diagnostics"]["diagnostics"].is_array());
    assert!(
        value["suggestions"]
            .as_array()
            .is_some_and(|list| !list.is_empty())
    );
}

#[test]
fn a_pipeline_reads_stdin_and_writes_stdout() {
    let first = run(&["supercell", SILICON, "2", "2", "2", "--out", "-"]);
    assert!(first.status.success());
    let mut child = Command::new(env!("CARGO_BIN_EXE_xtalforge"))
        .args([
            "slab", "-", "1", "1", "1", "--layers", "1", "--vacuum", "10", "--out", "-",
        ])
        .stdin(Stdio::piped())
        .stdout(Stdio::piped())
        .stderr(Stdio::piped())
        .spawn()
        .expect("spawn");
    child
        .stdin
        .as_mut()
        .expect("stdin")
        .write_all(&first.stdout)
        .expect("write stdin");
    let second = child.wait_with_output().expect("wait");
    assert!(
        second.status.success(),
        "stderr: {}",
        String::from_utf8_lossy(&second.stderr)
    );
    let value: serde_json::Value = serde_json::from_slice(&second.stdout).expect("JSON");
    assert_eq!(value["schema_version"], 1);
}

#[test]
fn json_errors_carry_a_category() {
    let output = run(&["--json", "info", "/nonexistent/file.cif"]);
    assert!(!output.status.success());
    let value: serde_json::Value = serde_json::from_slice(&output.stderr).expect("JSON");
    assert_eq!(value["error"]["code"], "invalid_input");
    assert_eq!(value["error"]["category"], "input");
}

#[test]
fn molecule_sketch_fragment_and_guest_commands_run() {
    let molecule = run(&["--json", "molecule", "--smiles", "c1ccccc1", "--out", "-"]);
    assert!(molecule.status.success());
    let value: serde_json::Value = serde_json::from_slice(&molecule.stdout).expect("JSON");
    assert_eq!(value["schema_version"], 1);

    let sketch = run(&["--json", "sketch", "c1ccccc1"]);
    assert!(sketch.status.success());
    let value: serde_json::Value = serde_json::from_slice(&sketch.stdout).expect("JSON");
    assert!(
        value["roundtrip"]
            .as_str()
            .is_some_and(|text| !text.is_empty())
    );

    let fragments = run(&["--json", "fragment"]);
    assert!(fragments.status.success());
    let value: serde_json::Value = serde_json::from_slice(&fragments.stdout).expect("JSON");
    assert!(value.as_array().is_some_and(|list| !list.is_empty()));

    let guest = run(&[
        "guest",
        SILICON,
        "--smiles",
        "O",
        "--offset",
        "2.7,2.7,2.7",
        "--clearance",
        "0.3",
        "--out",
        "-",
    ]);
    assert!(guest.status.success());
}
