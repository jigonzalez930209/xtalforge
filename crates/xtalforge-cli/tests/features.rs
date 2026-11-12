//! Integration tests for the Phase 9 CLI commands.

#![allow(
    clippy::unwrap_used,
    clippy::expect_used,
    clippy::panic,
    clippy::indexing_slicing
)]

use std::process::{Command, Output};

const SILICON: &str = concat!(env!("CARGO_MANIFEST_DIR"), "/tests/fixtures/silicon.poscar");
const SILICON_CIF: &str = concat!(
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
fn compare_reports_a_match() {
    let value = json(run(&["--json", "compare", SILICON, SILICON]));
    assert_eq!(value["matched"], true);
    assert!(value["rms"].as_f64().unwrap_or(1.0) < 1e-6);
}

#[test]
fn group_puts_equal_structures_in_one_class() {
    let value = json(run(&["--json", "group", SILICON, SILICON]));
    assert_eq!(value["n_groups"], 1);
    assert_eq!(value["n_structures"], 2);
}

#[test]
fn deform_scales_the_cell() {
    let value = json(run(&["--json", "deform", SILICON, "--volume", "200"]));
    assert_eq!(value["schema_version"], 1);
    assert_eq!(
        value["structure"]["sites"].as_array().map(Vec::len),
        Some(2)
    );
}

#[test]
fn features_miller_lists_the_cubic_families() {
    let value = json(run(&[
        "--json",
        "features",
        "miller",
        SILICON_CIF,
        "--max-index",
        "1",
    ]));
    let indices = value.as_array().expect("list");
    assert_eq!(indices.len(), 3);
    assert!(
        indices
            .iter()
            .any(|index| index["hkl"] == serde_json::json!([1, 0, 0]))
    );
}

#[test]
fn features_eos_fits_the_points() {
    let value = json(run(&[
        "--json", "features", "eos", "--point", "18,1", "--point", "19,0.4", "--point", "20,0",
        "--point", "21,0.5", "--point", "22,1.6",
    ]));
    let v0 = value["v0"].as_f64().expect("v0");
    assert!((18.0..22.0).contains(&v0));
}
