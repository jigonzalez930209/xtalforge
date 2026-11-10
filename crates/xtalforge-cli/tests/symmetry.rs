//! Integration tests for the symmetry command surface (Task 2.1-2.3).

#![allow(
    clippy::unwrap_used,
    clippy::expect_used,
    clippy::panic,
    clippy::indexing_slicing
)]

use std::io::Write;
use std::path::Path;
use std::process::{Command, Output, Stdio};

fn run(args: &[&str]) -> Output {
    Command::new(env!("CARGO_BIN_EXE_xtalforge"))
        .args(args)
        .output()
        .expect("the CLI binary must run")
}

fn run_with_stdin(args: &[&str], input: &str) -> Output {
    let mut child = Command::new(env!("CARGO_BIN_EXE_xtalforge"))
        .args(args)
        .stdin(Stdio::piped())
        .stdout(Stdio::piped())
        .stderr(Stdio::piped())
        .spawn()
        .expect("the CLI binary must start");
    child
        .stdin
        .take()
        .expect("stdin must be piped")
        .write_all(input.as_bytes())
        .expect("stdin must accept the input");
    child.wait_with_output().expect("the CLI binary must run")
}

#[test]
fn symmetry_detect_reports_the_space_group() {
    let path = Path::new(env!("CARGO_MANIFEST_DIR")).join("../../tests/corpus/simple/Si.cif");
    let output = run(&[
        "--json",
        "symmetry",
        "detect",
        &path.to_string_lossy(),
        "--symprec",
        "0.001",
    ]);
    assert!(output.status.success());
    let value: serde_json::Value =
        serde_json::from_slice(&output.stdout).expect("stdout must be JSON");
    assert_eq!(value["number"], 227);
    assert_eq!(value["hall_number"], 525);
    assert_eq!(
        value["operations"].as_array().expect("operations").len(),
        192
    );
}

#[test]
fn spacegroup_show_and_list() {
    let output = run(&["--json", "spacegroup", "show", "225"]);
    assert!(output.status.success());
    let value: serde_json::Value =
        serde_json::from_slice(&output.stdout).expect("stdout must be JSON");
    assert_eq!(value["number"], 225);
    assert_eq!(value["point_group"], "m-3m");
    assert_eq!(value["bravais"], "cF");
    assert_eq!(
        value["operations"].as_array().expect("operations").len(),
        192
    );

    let output = run(&[
        "--json",
        "spacegroup",
        "list",
        "--system",
        "cubic",
        "--centering",
        "F",
    ]);
    assert!(output.status.success());
    let value: serde_json::Value =
        serde_json::from_slice(&output.stdout).expect("stdout must be JSON");
    let items = value.as_array().expect("array");
    assert!(!items.is_empty());
    assert!(items.iter().all(|item| item["centering"] == "F"));

    let output = run(&["--json", "spacegroup", "search", "Fd"]);
    assert!(output.status.success());
    let value: serde_json::Value =
        serde_json::from_slice(&output.stdout).expect("stdout must be JSON");
    assert!(
        value
            .as_array()
            .expect("array")
            .iter()
            .any(|item| item["number"] == 227)
    );
}

#[test]
fn symmetry_p1_and_asymmetrize() {
    let path = Path::new(env!("CARGO_MANIFEST_DIR")).join("../../tests/corpus/simple/Si.cif");
    let output = run(&["--json", "symmetry", "p1", &path.to_string_lossy()]);
    assert!(output.status.success());
    let value: serde_json::Value =
        serde_json::from_slice(&output.stdout).expect("stdout must be JSON");
    assert_eq!(
        value["structure"]["sites"].as_array().expect("sites").len(),
        8
    );
    assert_eq!(value["structure"]["hall"], "P 1");

    let output = run(&[
        "--json",
        "symmetry",
        "asymmetrize",
        &path.to_string_lossy(),
        "--symprec",
        "0.001",
    ]);
    assert!(output.status.success());
    let value: serde_json::Value =
        serde_json::from_slice(&output.stdout).expect("stdout must be JSON");
    assert_eq!(
        value["structure"]["sites"].as_array().expect("sites").len(),
        1
    );
    assert_eq!(value["structure"]["hall"], "F 4d 2 3 -1d");
}

#[test]
fn symmetry_set_impose() {
    let path = Path::new(env!("CARGO_MANIFEST_DIR")).join("../../tests/corpus/simple/Si.cif");
    let output = run(&[
        "--json",
        "symmetry",
        "asymmetrize",
        &path.to_string_lossy(),
        "--symprec",
        "0.001",
    ]);
    assert!(output.status.success());
    let asymmetric = String::from_utf8(output.stdout).expect("UTF-8");
    let document: serde_json::Value = serde_json::from_str(&asymmetric).expect("JSON");
    let structure = serde_json::to_string(&document["structure"]).expect("structure");
    let output = run_with_stdin(
        &[
            "--json",
            "symmetry",
            "set",
            "-",
            "227",
            "--mode",
            "reinterpret",
        ],
        &structure,
    );
    assert!(output.status.success());
    let value: serde_json::Value =
        serde_json::from_slice(&output.stdout).expect("stdout must be JSON");
    assert_eq!(value["sites"].as_array().expect("sites").len(), 1);
}

#[test]
fn symmetry_subgroups_and_lowering() {
    let output = run(&["--json", "symmetry", "subgroups", "525"]);
    assert!(output.status.success());
    let value: serde_json::Value =
        serde_json::from_slice(&output.stdout).expect("stdout must be JSON");
    let groups = value.as_array().expect("array");
    assert!(!groups.is_empty());
    assert!(groups.iter().any(|group| group["subgroup_number"] == 216));

    let path = Path::new(env!("CARGO_MANIFEST_DIR")).join("../../tests/corpus/simple/Si.cif");
    let output = run(&[
        "--json",
        "symmetry",
        "lower",
        &path.to_string_lossy(),
        "216",
    ]);
    assert!(output.status.success());
    let value: serde_json::Value =
        serde_json::from_slice(&output.stdout).expect("stdout must be JSON");
    assert_eq!(value["structure"]["hall"], "F -4 2 3");
}

#[test]
fn symmetry_analysis_commands() {
    let path = Path::new(env!("CARGO_MANIFEST_DIR")).join("../../tests/corpus/simple/Si.cif");
    let path = path.to_string_lossy();

    let output = run(&[
        "--json",
        "symmetry",
        "standardize",
        &path,
        "--symprec",
        "0.001",
    ]);
    assert!(output.status.success());
    let value: serde_json::Value = serde_json::from_slice(&output.stdout).expect("JSON");
    assert_eq!(
        value["structure"]["sites"].as_array().expect("sites").len(),
        1
    );

    let output = run(&["--json", "symmetry", "primitive", &path]);
    assert!(output.status.success());

    let output = run(&[
        "--json",
        "symmetry",
        "reduce",
        &path,
        "--algorithm",
        "niggli",
    ]);
    assert!(output.status.success());
    let value: serde_json::Value = serde_json::from_slice(&output.stdout).expect("JSON");
    assert!(value["structure"]["lattice"].is_array());

    let output = run(&["--json", "symmetry", "wyckoff", &path, "--symprec", "0.001"]);
    assert!(output.status.success());
    let value: serde_json::Value = serde_json::from_slice(&output.stdout).expect("JSON");
    assert!(
        value
            .as_array()
            .expect("array")
            .iter()
            .all(|item| item["letter"] == "b")
    );

    let output = run(&["--json", "symmetry", "brillouin", &path]);
    assert!(output.status.success());
    let value: serde_json::Value = serde_json::from_slice(&output.stdout).expect("JSON");
    assert_eq!(value["faces"].as_array().expect("faces").len(), 14);
    assert_eq!(value["vertices"].as_array().expect("vertices").len(), 24);

    let output = run(&["--json", "symmetry", "scan", &path]);
    assert!(output.status.success());

    let output = run(&["--json", "symmetry", "diagnose", &path, "227"]);
    assert!(output.status.success());
    let value: serde_json::Value = serde_json::from_slice(&output.stdout).expect("JSON");
    assert!(value["max_displacement"].as_f64().expect("max") < 1e-6);
}
