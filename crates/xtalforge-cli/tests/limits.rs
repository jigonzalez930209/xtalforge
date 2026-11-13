//! Integration tests for the `--limits` presets (Task 15.2.5).

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

fn fixture(name: &str) -> String {
    Path::new(env!("CARGO_MANIFEST_DIR"))
        .join("tests/fixtures")
        .join(name)
        .to_string_lossy()
        .into_owned()
}

fn stdout_json(output: &Output) -> serde_json::Value {
    assert!(output.status.success());
    serde_json::from_slice(&output.stdout).expect("stdout must be JSON")
}

#[test]
fn the_huge_limits_preset_is_accepted() {
    let output = run(&["--json", "--limits", "huge", "info", &fixture("water.xyz")]);
    assert!(output.status.success());
    let value = stdout_json(&output);
    assert_eq!(value["n_sites"], 3);
    let help = run(&["--help"]);
    let text = String::from_utf8(help.stdout).expect("stdout is UTF-8");
    assert!(text.contains("huge"), "help does not list the huge preset");
}
