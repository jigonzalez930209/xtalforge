//! In-memory workflow runner (Task 17.2.5).

#![allow(
    clippy::unwrap_used,
    clippy::expect_used,
    clippy::panic,
    clippy::indexing_slicing
)]

use std::collections::BTreeMap;

use xtalforge_agent::workflow::{ItemStatus, Workflow, run_workflow_memory};

const NACL: &str = include_str!("../../../tests/corpus/simple/NaCl.cif");

fn workflow() -> Workflow {
    Workflow::from_json(
        r#"{
          "workflow_version": 1,
          "name": "memory",
          "steps": [
            {"kind": "read", "path": "{{input}}"},
            {"kind": "analyse", "verb": "info", "name": "summary"},
            {"kind": "write", "path": "out/{{stem}}.cif"}
          ]
        }"#,
    )
    .unwrap()
}

#[test]
fn runs_without_touching_the_filesystem() {
    let mut files = BTreeMap::new();
    files.insert("NaCl.cif".to_owned(), NACL.to_owned());
    let result = run_workflow_memory(&workflow(), &files).unwrap();
    assert_eq!(result.report.name, "memory");
    assert_eq!(result.report.inputs, 1);
    assert_eq!(result.report.succeeded, 1);
    assert_eq!(result.report.failed, 0);
    assert_eq!(result.report.jobs, 1);
    let item = result.report.items.first().unwrap();
    assert_eq!(item.status, ItemStatus::Ok);
    assert_eq!(item.steps.len(), 3);
    assert_eq!(item.steps[1].step, "summary");
    let text = result.outputs.get("out/NaCl.cif").expect("output");
    assert!(text.contains("_cell_length_a"));
    assert!(!std::path::Path::new("out/NaCl.cif").exists());
}

#[test]
fn a_missing_memory_input_is_reported_per_item() {
    let mut files = BTreeMap::new();
    files.insert("missing.cif".to_owned(), NACL.to_owned());
    let mut workflow = workflow();
    workflow.steps[0].path = Some("other.cif".to_owned());
    let result = run_workflow_memory(&workflow, &files).unwrap();
    assert_eq!(result.report.failed, 1);
    let item = result.report.items.first().unwrap();
    assert_eq!(item.status, ItemStatus::Failed);
    let error = item.error.as_ref().unwrap();
    assert_eq!(error.code, "io");
    assert!(error.message.contains("other.cif"));
}

#[test]
fn an_unknown_output_format_is_an_error() {
    let mut files = BTreeMap::new();
    files.insert("NaCl.cif".to_owned(), NACL.to_owned());
    let mut workflow = workflow();
    workflow.steps[2].path = Some("out/data.unknown".to_owned());
    let result = run_workflow_memory(&workflow, &files).unwrap();
    assert_eq!(result.report.failed, 1);
    assert!(result.outputs.is_empty());
}
