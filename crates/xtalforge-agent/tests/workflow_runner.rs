//! Workflow runner: order, isolation, resume and checkpoints (Task 16.1.2).

#![allow(
    clippy::unwrap_used,
    clippy::expect_used,
    clippy::panic,
    clippy::indexing_slicing
)]

use std::path::{Path, PathBuf};

use xtalforge_agent::workflow::{
    ItemStatus, RunOptions, Workflow, WorkflowCheckpoint, expand_inputs, run_workflow,
};

const NACL: &str = include_str!("../../../tests/corpus/simple/NaCl.cif");
const SI: &str = include_str!("../../../tests/corpus/simple/Si.cif");
const FIXTURE: &str = include_str!("fixtures/workflow.json");

fn setup(tag: &str) -> PathBuf {
    let path = std::env::temp_dir().join(format!("xtalforge-wf-run-{}-{tag}", std::process::id()));
    let _ = std::fs::remove_dir_all(&path);
    std::fs::create_dir_all(&path).expect("temp dir");
    std::fs::write(path.join("NaCl.cif"), NACL).expect("write");
    std::fs::write(path.join("Si.cif"), SI).expect("write");
    path
}

fn workflow() -> Workflow {
    Workflow::from_json(FIXTURE).expect("fixture")
}

fn options(directory: &Path, jobs: usize) -> RunOptions {
    let pattern = format!("{}/*.cif", directory.display());
    RunOptions {
        inputs: expand_inputs([pattern.as_str()]).expect("expand"),
        jobs,
        checkpoint: None,
        resume: true,
    }
}

#[test]
fn every_input_runs_in_deterministic_order() {
    let directory = setup("order");
    let report = run_workflow(&workflow(), &options(&directory, 1)).expect("run");
    assert_eq!(report.name, "inspect-copy");
    assert_eq!(report.inputs, 2);
    assert_eq!(report.succeeded, 2);
    assert_eq!(report.failed, 0);
    assert_eq!(report.resumed, 0);
    assert_eq!(report.jobs, 1);
    assert_eq!(report.items.len(), 2);
    assert!(report.items[0].input.ends_with("NaCl.cif"));
    assert!(report.items[1].input.ends_with("Si.cif"));
    for item in &report.items {
        assert_eq!(item.status, ItemStatus::Ok);
        assert_eq!(item.steps.len(), 3);
        assert_eq!(item.steps[1].step, "summary");
        assert_eq!(
            item.steps[1].answer.as_ref().expect("answer")["kind"],
            "info"
        );
        assert!(!item.resumed);
    }
    assert!(directory.join("out/NaCl.cif").is_file());
    assert!(directory.join("out/Si.cif").is_file());
    let _ = std::fs::remove_dir_all(&directory);
}

#[test]
fn parallel_jobs_keep_the_same_report_and_output() {
    let sequential = setup("seq");
    let parallel = setup("par");
    let first = run_workflow(&workflow(), &options(&sequential, 1)).expect("run");
    let second = run_workflow(&workflow(), &options(&parallel, 4)).expect("run");
    assert_eq!(second.jobs, 2);
    let order = |report: &xtalforge_agent::workflow::WorkflowReport| {
        report
            .items
            .iter()
            .map(|item| {
                let name = Path::new(&item.input)
                    .file_name()
                    .expect("name")
                    .to_string_lossy()
                    .into_owned();
                (name, item.status)
            })
            .collect::<Vec<_>>()
    };
    assert_eq!(order(&first), order(&second));
    assert_eq!(first.succeeded, second.succeeded);
    for name in ["NaCl.cif", "Si.cif"] {
        let a = std::fs::read_to_string(sequential.join("out").join(name)).expect("first output");
        let b = std::fs::read_to_string(parallel.join("out").join(name)).expect("second output");
        assert_eq!(a, b);
    }
    let _ = std::fs::remove_dir_all(&sequential);
    let _ = std::fs::remove_dir_all(&parallel);
}

#[test]
fn one_bad_item_does_not_stop_the_others() {
    let directory = setup("isolation");
    std::fs::write(directory.join("broken.cif"), "this is not a structure").expect("write");
    let report = run_workflow(&workflow(), &options(&directory, 3)).expect("run");
    assert_eq!(report.inputs, 3);
    assert_eq!(report.succeeded, 2);
    assert_eq!(report.failed, 1);
    let failed = report
        .items
        .iter()
        .find(|item| item.status == ItemStatus::Failed)
        .expect("failed item");
    assert!(failed.input.ends_with("broken.cif"));
    assert!(failed.steps.is_empty());
    let problem = failed.error.as_ref().expect("problem");
    assert!(!problem.code.is_empty());
    assert_eq!(problem.step.as_deref(), Some("step-1"));
    assert!(directory.join("out/NaCl.cif").is_file());
    let _ = std::fs::remove_dir_all(&directory);
}

#[test]
fn a_checkpoint_resumes_only_the_failed_items() {
    let directory = setup("resume");
    let checkpoint = directory.join("run.checkpoint.json");
    std::fs::write(directory.join("broken.cif"), "this is not a structure").expect("write");
    let mut first = options(&directory, 2);
    first.checkpoint = Some(checkpoint.clone());
    let report = run_workflow(&workflow(), &first).expect("run");
    assert_eq!(report.failed, 1);
    assert_eq!(report.resumed, 0);
    let saved = WorkflowCheckpoint::read(&checkpoint).expect("checkpoint");
    assert_eq!(saved.items.len(), 3);
    assert_eq!(saved.name, "inspect-copy");
    assert!(saved.fingerprint.starts_with("fnv1a64:"));
    assert!(saved.items.iter().all(|item| !item.resumed));

    // A resumed run must not rewrite the output of an item it skipped.
    std::fs::remove_file(directory.join("out/NaCl.cif")).expect("remove output");
    std::fs::write(directory.join("broken.cif"), SI).expect("fix input");
    let report = run_workflow(&workflow(), &first).expect("resume");
    assert_eq!(report.succeeded, 3);
    assert_eq!(report.failed, 0);
    assert_eq!(report.resumed, 2);
    assert!(!directory.join("out/NaCl.cif").is_file());
    assert!(directory.join("out/broken.cif").is_file());
    let resumed: Vec<bool> = report.items.iter().map(|item| item.resumed).collect();
    assert_eq!(resumed, vec![true, true, false]);

    // `--fresh` ignores the checkpoint and runs everything again.
    let mut fresh = first.clone();
    fresh.resume = false;
    let report = run_workflow(&workflow(), &fresh).expect("fresh");
    assert_eq!(report.resumed, 0);
    assert_eq!(report.succeeded, 3);
    assert!(directory.join("out/NaCl.cif").is_file());
    let _ = std::fs::remove_dir_all(&directory);
}

#[test]
fn a_checkpoint_from_another_workflow_is_refused() {
    let directory = setup("foreign");
    let checkpoint = directory.join("run.checkpoint.json");
    let mut first = options(&directory, 1);
    first.checkpoint = Some(checkpoint.clone());
    run_workflow(&workflow(), &first).expect("run");

    let other = Workflow::from_value(serde_json::json!({
        "workflow_version": 1,
        "name": "another",
        "steps": [{"kind": "read"}]
    }))
    .expect("workflow");
    let mut second = options(&directory, 1);
    second.checkpoint = Some(checkpoint);
    let error = run_workflow(&other, &second).expect_err("must refuse");
    assert_eq!(error.code(), "checkpoint");
    assert!(error.to_string().contains("another"));
    let _ = std::fs::remove_dir_all(&directory);
}

#[test]
fn an_empty_input_list_is_an_empty_report() {
    let report = run_workflow(
        &workflow(),
        &RunOptions {
            inputs: Vec::new(),
            jobs: 4,
            checkpoint: None,
            resume: true,
        },
    )
    .expect("run");
    assert_eq!(report.inputs, 0);
    assert_eq!(report.items.len(), 0);
    assert_eq!(report.jobs, 1);
}

#[test]
fn jobs_zero_is_sequential() {
    let directory = setup("jobs-zero");
    let report = run_workflow(&workflow(), &options(&directory, 0)).expect("run");
    assert_eq!(report.jobs, 1);
    assert_eq!(report.succeeded, 2);
    let _ = std::fs::remove_dir_all(&directory);
}

#[test]
fn write_paths_create_directories_and_use_one_based_indexes() {
    let directory = setup("template");
    let workflow = Workflow::from_value(serde_json::json!({
        "workflow_version": 1,
        "name": "indexed",
        "steps": [
            {"kind": "read"},
            {"kind": "write", "path": "{{dir}}/deep/nested/{{index}}-{{stem}}.{{ext}}"}
        ]
    }))
    .expect("workflow");
    let report = run_workflow(&workflow, &options(&directory, 2)).expect("run");
    assert_eq!(report.succeeded, 2);
    assert!(directory.join("deep/nested/1-NaCl.cif").is_file());
    assert!(directory.join("deep/nested/2-Si.cif").is_file());
    let _ = std::fs::remove_dir_all(&directory);
}

#[test]
fn duplicates_in_the_input_list_are_removed() {
    let directory = setup("duplicates");
    let file = directory.join("NaCl.cif");
    let options = RunOptions {
        inputs: vec![file.clone(), file.clone(), directory.join("Si.cif")],
        jobs: 2,
        checkpoint: None,
        resume: true,
    };
    let report = run_workflow(&workflow(), &options).expect("run");
    assert_eq!(report.inputs, 2);
    assert_eq!(report.succeeded, 2);
    let _ = std::fs::remove_dir_all(&directory);
}
