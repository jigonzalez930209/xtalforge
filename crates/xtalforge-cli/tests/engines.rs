//! Integration tests for the Phase 18 engine commands.

#![allow(
    clippy::unwrap_used,
    clippy::expect_used,
    clippy::panic,
    clippy::indexing_slicing
)]

use std::process::{Command, Output};

const NACL: &str = concat!(
    env!("CARGO_MANIFEST_DIR"),
    "/../../tests/corpus/simple/NaCl.cif"
);
const EIGENVAL: &str = concat!(
    env!("CARGO_MANIFEST_DIR"),
    "/../xtalforge-engines/tests/fixtures/engines/vasp/EIGENVAL"
);
const DOSCAR: &str = concat!(
    env!("CARGO_MANIFEST_DIR"),
    "/../xtalforge-engines/tests/fixtures/engines/vasp/DOSCAR"
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
fn presets_and_failures_are_listed() {
    let presets = json(run(&["--json", "engines", "presets"]));
    assert!(presets.as_array().is_some_and(|items| items.len() >= 20));
    let failures = json(run(&["--json", "engines", "failures"]));
    assert!(
        failures
            .as_array()
            .is_some_and(|items| items.iter().any(|item| item["code"] == "scf_not_converged"))
    );
}

#[test]
fn inputs_and_scripts_are_written() {
    let directory = std::env::temp_dir().join(format!(
        "xtalforge-engines-cli-{}-{}",
        std::process::id(),
        std::time::SystemTime::now()
            .duration_since(std::time::UNIX_EPOCH)
            .unwrap()
            .as_nanos()
    ));
    let out = directory.to_string_lossy().into_owned();
    let value = json(run(&[
        "--json",
        "engines",
        "inputs",
        "vasp/relax",
        NACL,
        "--out",
        &out,
    ]));
    assert_eq!(value["preset"], "vasp/relax");
    for name in ["INCAR", "KPOINTS", "POSCAR", "POTCAR.spec"] {
        assert!(directory.join(name).is_file(), "missing {name}");
    }
    let script = json(run(&[
        "--json",
        "engines",
        "script",
        "slurm",
        "xtalforge info a.cif",
        "--out",
        &out,
        "--dry-run",
    ]));
    assert!(script["path"].as_str().is_some());
    let _ = std::fs::remove_dir_all(&directory);
}

#[test]
fn parsers_and_electronic_outputs_are_typed() {
    let value = json(run(&[
        "--json",
        "engines",
        "parse",
        "vasp-eigenval",
        EIGENVAL,
    ]));
    assert!(value["bands"].is_object());
    let bands = json(run(&[
        "--json",
        "engines",
        "bands",
        "vasp-eigenval",
        EIGENVAL,
    ]));
    assert!(bands["bands"]["bands"].as_array().is_some());
    assert!(bands["gap"]["energy"].as_f64().is_some());
    let dos = json(run(&["--json", "engines", "dos", "vasp-doscar", DOSCAR]));
    assert!(dos["smeared"].as_array().is_some());
}

#[test]
fn unknown_codes_and_schedulers_are_errors() {
    let output = run(&["--json", "engines", "parse", "nope", EIGENVAL]);
    assert!(!output.status.success());
    let output = run(&["--json", "engines", "script", "nope", "true"]);
    assert!(!output.status.success());
}
