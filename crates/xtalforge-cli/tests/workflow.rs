//! `xtalforge run`: declarative workflows from the CLI (Task 16.1.2).

#![allow(
    clippy::unwrap_used,
    clippy::expect_used,
    clippy::panic,
    clippy::indexing_slicing
)]

use std::path::{Path, PathBuf};
use std::process::{Command, Output};

const NACL: &str = include_str!("../../../tests/corpus/simple/NaCl.cif");
const SI: &str = include_str!("../../../tests/corpus/simple/Si.cif");
const FIXTURE: &str = r#"{
  "workflow_version": 1,
  "name": "inspect-copy",
  "steps": [
    {"kind": "read", "path": "{{input}}"},
    {"kind": "analyse", "name": "summary", "verb": "info"},
    {"kind": "write", "name": "copy", "path": "{{dir}}/out/{{stem}}.cif"}
  ]
}
"#;

fn run(args: &[&str]) -> Output {
    Command::new(env!("CARGO_BIN_EXE_xtalforge"))
        .args(args)
        .output()
        .expect("the CLI binary must run")
}

fn temp_dir(tag: &str) -> PathBuf {
    let path = std::env::temp_dir().join(format!("xtalforge-cli-wf-{}-{tag}", std::process::id()));
    let _ = std::fs::remove_dir_all(&path);
    std::fs::create_dir_all(&path).expect("temp dir");
    path
}

fn write_workflow(directory: &Path) -> String {
    let path = directory.join("wf.json");
    std::fs::write(&path, FIXTURE).expect("workflow");
    path.to_string_lossy().into_owned()
}

fn stdout_json(output: &Output) -> serde_json::Value {
    serde_json::from_slice(&output.stdout).expect("stdout must be one JSON value")
}

#[test]
fn the_schema_is_printed_without_a_workflow() {
    let output = run(&["run", "--schema"]);
    assert!(
        output.status.success(),
        "stderr: {}",
        String::from_utf8_lossy(&output.stderr)
    );
    let schema: serde_json::Value = serde_json::from_slice(&output.stdout).expect("schema");
    assert_eq!(
        schema["$id"],
        xtalforge::agent::workflow::WORKFLOW_SCHEMA_ID
    );
    assert_eq!(schema["properties"]["workflow_version"]["const"], 1);
}

#[test]
fn validate_prints_the_steps_and_needs_no_inputs() {
    let directory = temp_dir("validate");
    let workflow = write_workflow(&directory);
    let output = run(&["--json", "run", &workflow, "--validate"]);
    assert!(output.status.success());
    let summary = stdout_json(&output);
    assert_eq!(summary["name"], "inspect-copy");
    assert_eq!(summary["workflow_version"], 1);
    assert_eq!(summary["steps"].as_array().expect("steps").len(), 3);
    let _ = std::fs::remove_dir_all(&directory);
}

#[test]
fn a_glob_runs_every_input_in_order_and_writes_the_report() {
    let directory = temp_dir("glob");
    let workflow = write_workflow(&directory);
    std::fs::write(directory.join("NaCl.cif"), NACL).expect("write");
    std::fs::write(directory.join("Si.cif"), SI).expect("write");
    let pattern = format!("{}/*.cif", directory.display());
    let report_path = directory.join("report.json");
    let output = run(&[
        "--json",
        "run",
        &workflow,
        "--inputs",
        &pattern,
        "--jobs",
        "2",
        "--out",
        &report_path.to_string_lossy(),
    ]);
    assert!(
        output.status.success(),
        "stderr: {}",
        String::from_utf8_lossy(&output.stderr)
    );
    let report = stdout_json(&output);
    assert_eq!(report["name"], "inspect-copy");
    assert_eq!(report["inputs"], 2);
    assert_eq!(report["succeeded"], 2);
    assert_eq!(report["failed"], 0);
    assert_eq!(report["jobs"], 2);
    assert!(
        report["items"][0]["input"]
            .as_str()
            .expect("input")
            .ends_with("NaCl.cif")
    );
    assert!(
        report["items"][1]["input"]
            .as_str()
            .expect("input")
            .ends_with("Si.cif")
    );
    assert!(directory.join("out/NaCl.cif").is_file());
    assert!(directory.join("out/Si.cif").is_file());
    let written: serde_json::Value =
        serde_json::from_str(&std::fs::read_to_string(&report_path).expect("report"))
            .expect("report JSON");
    assert_eq!(written["succeeded"], 2);
    let _ = std::fs::remove_dir_all(&directory);
}

#[test]
fn a_bad_workflow_exits_with_code_three_and_a_typed_error() {
    let directory = temp_dir("bad");
    let path = directory.join("bad.json");
    std::fs::write(
        &path,
        r#"{"workflow_version": 1, "name": "x", "steps": [{"kind": "write", "path": "out.cif"}]}"#,
    )
    .expect("write");
    let output = run(&[
        "--json",
        "run",
        &path.to_string_lossy(),
        "--inputs",
        "*.cif",
    ]);
    assert_eq!(output.status.code(), Some(3));
    let error: serde_json::Value = serde_json::from_slice(&output.stderr).expect("stderr JSON");
    assert_eq!(error["error"]["code"], "workflow_invalid");
    assert_eq!(error["error"]["category"], "input");
    assert!(
        error["error"]["message"]
            .as_str()
            .expect("message")
            .contains("steps[0].kind")
    );
    let _ = std::fs::remove_dir_all(&directory);
}

#[test]
fn a_failed_item_is_isolated_and_sets_the_exit_code() {
    let directory = temp_dir("isolated");
    let workflow = write_workflow(&directory);
    std::fs::write(directory.join("NaCl.cif"), NACL).expect("write");
    std::fs::write(directory.join("broken.cif"), "not a structure").expect("write");
    let pattern = format!("{}/*.cif", directory.display());
    let output = run(&["--json", "run", &workflow, "--inputs", &pattern]);
    assert_eq!(output.status.code(), Some(3));
    let report = stdout_json(&output);
    assert_eq!(report["succeeded"], 1);
    assert_eq!(report["failed"], 1);
    assert_eq!(report["inputs"], 2);
    let error: serde_json::Value = serde_json::from_slice(&output.stderr).expect("stderr JSON");
    assert_eq!(error["error"]["code"], "workflow_failed");
    assert_eq!(error["error"]["category"], "workflow");
    assert!(directory.join("out/NaCl.cif").is_file());
    let _ = std::fs::remove_dir_all(&directory);
}

#[test]
fn a_checkpoint_resumes_failed_items() {
    let directory = temp_dir("resume");
    let workflow = write_workflow(&directory);
    let checkpoint = directory.join("run.checkpoint.json");
    std::fs::write(directory.join("NaCl.cif"), NACL).expect("write");
    std::fs::write(directory.join("broken.cif"), "not a structure").expect("write");
    let pattern = format!("{}/*.cif", directory.display());
    let checkpoint_text = checkpoint.to_string_lossy().into_owned();
    let first = run(&[
        "--json",
        "run",
        &workflow,
        "--inputs",
        &pattern,
        "--checkpoint",
        &checkpoint_text,
    ]);
    assert_eq!(first.status.code(), Some(3));
    let report = stdout_json(&first);
    assert_eq!(report["failed"], 1);
    assert_eq!(report["resumed"], 0);
    assert!(checkpoint.is_file());

    std::fs::write(directory.join("broken.cif"), SI).expect("fix");
    let second = run(&[
        "--json",
        "run",
        &workflow,
        "--inputs",
        &pattern,
        "--checkpoint",
        &checkpoint_text,
    ]);
    assert!(
        second.status.success(),
        "stderr: {}",
        String::from_utf8_lossy(&second.stderr)
    );
    let report = stdout_json(&second);
    assert_eq!(report["succeeded"], 2);
    assert_eq!(report["failed"], 0);
    assert_eq!(report["resumed"], 1);

    let third = run(&[
        "--json",
        "run",
        &workflow,
        "--inputs",
        &pattern,
        "--checkpoint",
        &checkpoint_text,
        "--fresh",
    ]);
    assert!(third.status.success());
    let report = stdout_json(&third);
    assert_eq!(report["resumed"], 0);
    let _ = std::fs::remove_dir_all(&directory);
}

#[test]
fn missing_inputs_and_bad_jobs_are_input_errors() {
    let directory = temp_dir("arguments");
    let workflow = write_workflow(&directory);
    let output = run(&["--json", "run", &workflow]);
    assert_eq!(output.status.code(), Some(3));
    let error: serde_json::Value = serde_json::from_slice(&output.stderr).expect("stderr JSON");
    assert_eq!(error["error"]["category"], "input");

    let empty = format!("{}/none/*.cif", directory.display());
    let output = run(&[
        "--json", "run", &workflow, "--inputs", &empty, "--jobs", "0",
    ]);
    assert_eq!(output.status.code(), Some(3));
    let _ = std::fs::remove_dir_all(&directory);
}

#[test]
fn run_help_lists_the_flags() {
    let output = run(&["run", "--help"]);
    assert!(output.status.success());
    let text = String::from_utf8(output.stdout).expect("UTF-8");
    for flag in [
        "--inputs",
        "--jobs",
        "--checkpoint",
        "--out",
        "--fresh",
        "--validate",
        "--schema",
    ] {
        assert!(text.contains(flag), "missing {flag}");
    }
}
