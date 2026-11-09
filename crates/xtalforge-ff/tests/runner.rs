//! Process runner tests with a simulated executable (Task 4.9).

#![cfg(unix)]
#![allow(
    clippy::unwrap_used,
    clippy::expect_used,
    clippy::panic,
    clippy::indexing_slicing
)]

use std::fs;
use std::os::unix::fs::PermissionsExt;
use std::path::PathBuf;
use std::sync::Arc;
use std::sync::atomic::{AtomicBool, Ordering};
use std::time::Duration;

use xtalforge_ff::{RunSpec, RunStatus, run_in_scratch};

fn script(name: &str, body: &str) -> PathBuf {
    let path =
        std::env::temp_dir().join(format!("xtalforge-runner-{}-{name}.sh", std::process::id()));
    fs::write(&path, format!("#!/bin/sh\n{body}\n")).unwrap();
    let mut permissions = fs::metadata(&path).unwrap().permissions();
    permissions.set_mode(0o755);
    fs::set_permissions(&path, permissions).unwrap();
    path
}

#[test]
fn captures_logs_and_inputs_in_the_scratch_folder() {
    let program = script("ok", "echo hello\ncat input.txt\necho warning >&2\nexit 0");
    let mut spec = RunSpec::new(program.to_string_lossy(), Vec::new());
    spec.input_files
        .insert("input.txt".to_owned(), "payload".to_owned());
    spec.timeout = Duration::from_secs(5);
    spec.keep_scratch = true;
    let outcome = run_in_scratch(&spec, &AtomicBool::new(false)).unwrap();
    assert_eq!(outcome.status, RunStatus::Success);
    assert_eq!(outcome.exit_code, Some(0));
    assert!(outcome.stdout.contains("hello"));
    assert!(outcome.stdout.contains("payload"));
    assert!(outcome.stderr.contains("warning"));
    assert!(!outcome.stdout_truncated);
    assert!(outcome.scratch.join("input.txt").exists());
    fs::remove_dir_all(&outcome.scratch).unwrap();
    fs::remove_file(program).unwrap();
}

#[test]
fn timeout_kills_the_process() {
    let program = script("timeout", "sleep 5");
    let mut spec = RunSpec::new(program.to_string_lossy(), Vec::new());
    spec.timeout = Duration::from_millis(100);
    let outcome = run_in_scratch(&spec, &AtomicBool::new(false)).unwrap();
    assert_eq!(outcome.status, RunStatus::TimedOut);
    assert!(outcome.duration < Duration::from_secs(2));
    fs::remove_file(program).unwrap();
}

#[test]
fn cancellation_kills_the_process() {
    let program = script("cancel", "sleep 5");
    let mut spec = RunSpec::new(program.to_string_lossy(), Vec::new());
    spec.timeout = Duration::from_secs(5);
    let cancel = Arc::new(AtomicBool::new(false));
    let flag = Arc::clone(&cancel);
    std::thread::spawn(move || {
        std::thread::sleep(Duration::from_millis(100));
        flag.store(true, Ordering::Relaxed);
    });
    let outcome = run_in_scratch(&spec, &cancel).unwrap();
    assert_eq!(outcome.status, RunStatus::Cancelled);
    assert!(outcome.duration < Duration::from_secs(2));
    fs::remove_file(program).unwrap();
}

#[test]
fn logs_are_capped() {
    let program = script("logs", "printf 'abcdefghij'");
    let mut spec = RunSpec::new(program.to_string_lossy(), Vec::new());
    spec.max_log_bytes = 4;
    let outcome = run_in_scratch(&spec, &AtomicBool::new(false)).unwrap();
    assert_eq!(outcome.stdout, "abcd");
    assert!(outcome.stdout_truncated);
    fs::remove_file(program).unwrap();
}

#[test]
fn unsafe_input_paths_are_refused() {
    let mut spec = RunSpec::new("echo", vec!["x".to_owned()]);
    spec.input_files
        .insert("../outside.txt".to_owned(), "x".to_owned());
    assert!(run_in_scratch(&spec, &AtomicBool::new(false)).is_err());
}

#[test]
fn a_missing_program_is_an_error() {
    let spec = RunSpec::new("xtalforge-program-that-does-not-exist", Vec::new());
    assert!(run_in_scratch(&spec, &AtomicBool::new(false)).is_err());
}

#[test]
fn scratch_is_removed_by_default() {
    let program = script("clean", "echo done");
    let spec = RunSpec::new(program.to_string_lossy(), Vec::new());
    let outcome = run_in_scratch(&spec, &AtomicBool::new(false)).unwrap();
    assert_eq!(outcome.status, RunStatus::Success);
    assert!(!outcome.scratch.exists());
    fs::remove_file(program).unwrap();
}
