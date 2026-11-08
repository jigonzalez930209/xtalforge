//! Integration tests for the Phase 1 command surface (Task 1.17.1).

#![allow(
    clippy::unwrap_used,
    clippy::expect_used,
    clippy::panic,
    clippy::indexing_slicing
)]

use std::io::Write;
use std::path::{Path, PathBuf};
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

fn fixture(name: &str) -> String {
    Path::new(env!("CARGO_MANIFEST_DIR"))
        .join("tests/fixtures")
        .join(name)
        .to_string_lossy()
        .into_owned()
}

fn corpus(relative: &str) -> String {
    Path::new(env!("CARGO_MANIFEST_DIR"))
        .join("../../tests/corpus")
        .join(relative)
        .to_string_lossy()
        .into_owned()
}

fn temp_path(name: &str) -> PathBuf {
    std::env::temp_dir().join(format!("xtalforge-cli-{}-{name}", std::process::id()))
}

fn stdout_json(output: &Output) -> serde_json::Value {
    assert!(
        output.status.success(),
        "stderr: {}",
        String::from_utf8_lossy(&output.stderr)
    );
    serde_json::from_slice(&output.stdout).expect("stdout must be a single JSON value")
}

fn assert_close(value: Option<&serde_json::Value>, expected: f64, tolerance: f64) {
    let value = value.and_then(serde_json::Value::as_f64).expect("number");
    assert!(
        (value - expected).abs() <= tolerance,
        "got {value}, expected {expected} +/- {tolerance}"
    );
}

#[test]
fn convert_cif_to_poscar_writes_a_file() {
    let out = temp_path("convert.poscar");
    let out_text = out.to_string_lossy().into_owned();
    let output = run(&[
        "convert",
        &corpus("simple/Si.cif"),
        "--format",
        "poscar",
        "--out",
        &out_text,
    ]);
    assert!(output.status.success());
    let text = std::fs::read_to_string(&out).expect("the converted file must exist");
    assert!(text.contains("Direct"));
    assert!(text.contains("Si"));
    let human = String::from_utf8(output.stdout).expect("stdout is UTF-8");
    assert!(human.contains("wrote poscar"));
    std::fs::remove_file(out).expect("the temporary file must be removable");
}

#[test]
fn convert_json_to_stdout_carries_the_content() {
    let output = run(&[
        "--json",
        "convert",
        &fixture("silicon.poscar"),
        "--format",
        "xyz",
    ]);
    let value = stdout_json(&output);
    assert_eq!(value["format"], "xyz");
    assert!(
        value["content"]
            .as_str()
            .expect("content")
            .contains("Lattice=")
    );
    assert!(value["bytes"].as_u64().expect("bytes") > 0);
}

#[test]
fn convert_unknown_input_exits_with_code_three() {
    let output = run(&[
        "--json",
        "convert",
        "/definitely/not/here.unknown",
        "--format",
        "cif",
    ]);
    assert_eq!(output.status.code(), Some(3));
    let value: serde_json::Value =
        serde_json::from_slice(&output.stderr).expect("stderr must be JSON");
    assert_eq!(value["error"]["code"], "invalid_input");
}

#[test]
fn supercell_json_is_a_schema_document() {
    let output = run(&[
        "--json",
        "supercell",
        &corpus("simple/Si.cif"),
        "2",
        "1",
        "1",
    ]);
    let value = stdout_json(&output);
    assert_eq!(value["schema_version"], 1);
    assert_eq!(
        value["structure"]["sites"]
            .as_array()
            .expect("sites array")
            .len(),
        2
    );
}

#[test]
fn supercell_out_infers_the_format_from_the_extension() {
    let out = temp_path("supercell.xyz");
    let out_text = out.to_string_lossy().into_owned();
    let output = run(&[
        "--quiet",
        "supercell",
        &corpus("simple/Si.cif"),
        "2",
        "1",
        "1",
        "--out",
        &out_text,
    ]);
    assert!(output.status.success());
    assert!(output.stdout.is_empty(), "quiet writes no confirmation");
    let text = std::fs::read_to_string(&out).expect("the supercell file must exist");
    assert!(text.starts_with("2\n"), "got: {text}");
    std::fs::remove_file(out).expect("the temporary file must be removable");
}

#[test]
fn slab_json_contains_the_cut_structure() {
    let output = run(&[
        "--json",
        "slab",
        &corpus("simple/Si.cif"),
        "1",
        "0",
        "0",
        "--layers",
        "2",
        "--vacuum",
        "10",
    ]);
    let value = stdout_json(&output);
    assert_eq!(value["schema_version"], 1);
    let sites = value["structure"]["sites"]
        .as_array()
        .expect("sites array")
        .len();
    assert!(sites > 0);
    assert_eq!(value["structure"]["hall"], "P 1");
}

#[test]
fn bonds_reports_the_water_bonds() {
    let output = run(&["--json", "bonds", &fixture("water.xyz")]);
    let value = stdout_json(&output);
    assert_eq!(value["n_atoms"], 3);
    assert_eq!(value["n_bonds"], 2);
    for bond in value["bonds"].as_array().expect("bonds array") {
        assert!(bond["distance"].as_f64().expect("distance") < 1.2);
    }
}

#[test]
fn rings_finds_the_six_membered_ring() {
    let output = run(&[
        "--json",
        "rings",
        &fixture("hexagon.xyz"),
        "--max-size",
        "6",
    ]);
    let value = stdout_json(&output);
    assert_eq!(value["max_size"], 6);
    assert!(value["n_rings"].as_u64().expect("n_rings") >= 1);
    assert_eq!(value["rings"][0]["size"], 6);
    assert_eq!(
        value["rings"][0]["atoms"]
            .as_array()
            .expect("atoms array")
            .len(),
        6
    );
}

#[test]
fn groups_finds_the_methanol_hydroxyl() {
    let output = run(&["--json", "groups", &fixture("methanol.xyz")]);
    let value = stdout_json(&output);
    let kinds: Vec<&str> = value["groups"]
        .as_array()
        .expect("groups array")
        .iter()
        .filter_map(|group| group["kind"].as_str())
        .collect();
    assert!(kinds.contains(&"hydroxyl"), "got {kinds:?}");
}

#[test]
fn neighbors_respects_the_cutoff() {
    let output = run(&[
        "--json",
        "neighbors",
        &fixture("water.xyz"),
        "--cutoff",
        "2.0",
    ]);
    let value = stdout_json(&output);
    assert!(value["n_neighbors"].as_u64().expect("n_neighbors") > 0);
    for neighbor in value["neighbors"].as_array().expect("neighbors array") {
        assert!(neighbor["distance"].as_f64().expect("distance") <= 2.0 + 1e-9);
    }
}

#[test]
fn distance_reports_the_water_oh_length() {
    let output = run(&["--json", "distance", &fixture("water.xyz"), "0", "1"]);
    let value = stdout_json(&output);
    assert_eq!(value["i"], 0);
    assert_eq!(value["j"], 1);
    assert_close(value.get("distance"), 0.96, 1e-9);
}

#[test]
fn angle_reports_the_water_bond_angle() {
    let output = run(&["--json", "angle", &fixture("water.xyz"), "1", "0", "2"]);
    let value = stdout_json(&output);
    assert_close(value.get("angle"), 104.470_294_100_065_9, 1e-9);
}

#[test]
fn torsion_reports_the_peroxide_dihedral() {
    let output = run(&[
        "--json",
        "torsion",
        &fixture("hydrogen_peroxide.xyz"),
        "2",
        "0",
        "1",
        "3",
    ]);
    let value = stdout_json(&output);
    assert_close(value.get("torsion"), -146.148_989_696_441_8, 1e-9);
}

#[test]
fn human_measure_output_is_a_number() {
    let output = run(&["distance", &fixture("water.xyz"), "0", "1"]);
    assert!(output.status.success());
    let text = String::from_utf8(output.stdout).expect("stdout is UTF-8");
    assert_eq!(text.trim(), "0.960000");
}

#[test]
fn quiet_hides_analysis_headers() {
    let output = run(&["--quiet", "rings", &fixture("hexagon.xyz")]);
    assert!(output.status.success());
    let text = String::from_utf8(output.stdout).expect("stdout is UTF-8");
    assert!(!text.starts_with('#'), "got: {text}");
    assert!(text.starts_with("6 "), "got: {text}");
}

#[test]
fn info_reads_poscar_and_xyz() {
    let output = run(&["--json", "info", &fixture("silicon.poscar")]);
    let value = stdout_json(&output);
    assert_eq!(value["formula"], "Si2");
    assert_eq!(value["n_sites"], 2);

    let output = run(&["--json", "info", &fixture("water.xyz")]);
    let value = stdout_json(&output);
    assert_eq!(value["n_sites"], 3);
    assert_eq!(value["formula"], "O H2");
}

#[test]
fn info_reads_cif_from_standard_input() {
    let cif = include_str!("../../../tests/corpus/simple/Si.cif");
    let output = run_with_stdin(&["--json", "info", "-"], cif);
    let value = stdout_json(&output);
    assert_eq!(value["n_sites"], 1);
}

#[test]
fn convert_reads_xyz_from_standard_input() {
    let xyz = include_str!("fixtures/water.xyz");
    let output = run_with_stdin(&["convert", "-", "--format", "cif"], xyz);
    assert!(output.status.success());
    let text = String::from_utf8(output.stdout).expect("stdout is UTF-8");
    assert!(text.starts_with("data_"));
    assert!(text.contains("_atom_site_fract_x"));
}

#[test]
fn invalid_ring_size_exits_with_code_three_in_json() {
    let output = run(&["--json", "rings", &fixture("water.xyz"), "--max-size", "2"]);
    assert_eq!(output.status.code(), Some(3));
    let value: serde_json::Value =
        serde_json::from_slice(&output.stderr).expect("stderr must be JSON");
    assert_eq!(value["error"]["code"], "invalid_input");
    assert!(
        value["error"]["message"]
            .as_str()
            .expect("message")
            .contains("at least 3")
    );
}

#[test]
fn out_of_range_atom_index_exits_with_code_three() {
    let output = run(&["distance", &fixture("water.xyz"), "0", "9"]);
    assert_eq!(output.status.code(), Some(3));
    let text = String::from_utf8(output.stderr).expect("stderr is UTF-8");
    assert!(text.contains("outside the P1 cell"));
}

#[test]
fn global_flags_work_after_the_subcommand() {
    let output = run(&[
        "supercell",
        &fixture("water.xyz"),
        "1",
        "1",
        "1",
        "--json",
        "--limits",
        "generous",
        "--threads",
        "2",
    ]);
    let value = stdout_json(&output);
    assert_eq!(value["schema_version"], 1);
}

#[test]
fn zero_threads_exits_with_code_three() {
    let output = run(&["--threads", "0", "version"]);
    assert_eq!(output.status.code(), Some(3));
}

#[test]
fn help_lists_the_phase_one_commands() {
    let output = run(&["--help"]);
    assert!(output.status.success());
    let text = String::from_utf8(output.stdout).expect("stdout is UTF-8");
    for command in [
        "convert",
        "supercell",
        "slab",
        "bonds",
        "rings",
        "groups",
        "neighbors",
        "distance",
        "angle",
        "torsion",
    ] {
        assert!(text.contains(command), "help is missing {command}");
    }
}
