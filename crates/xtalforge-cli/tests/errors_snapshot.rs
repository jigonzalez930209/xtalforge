//! Snapshot tests for the JSON error rendering of the CLI (Task 17.1.5).
//!
//! Each fixture in `tests/fixtures/errors/` records one stable failure shape:
//! the command arguments, the exit code and the expected JSON object keys,
//! error `code` and category. Volatile text (absolute paths) is skipped; the
//! `message_contains` substrings keep the wording reviewable without freezing
//! a filesystem path.

#![allow(
    clippy::unwrap_used,
    clippy::expect_used,
    clippy::panic,
    clippy::indexing_slicing
)]

use std::collections::BTreeSet;
use std::fs;
use std::path::{Path, PathBuf};
use std::process::Command;
use std::time::{SystemTime, UNIX_EPOCH};

use serde::Deserialize;

const NACL: &str = include_str!("../../../tests/corpus/simple/NaCl.cif");

/// Every error snapshot, paired with the fixture file name.
const SNAPSHOTS: &[(&str, &str)] = &[
    (
        "unreadable_input",
        include_str!("fixtures/errors/unreadable_input.json"),
    ),
    (
        "unknown_format",
        include_str!("fixtures/errors/unknown_format.json"),
    ),
    (
        "invalid_index",
        include_str!("fixtures/errors/invalid_index.json"),
    ),
    (
        "invalid_numeric",
        include_str!("fixtures/errors/invalid_numeric.json"),
    ),
    (
        "mof_failure",
        include_str!("fixtures/errors/mof_failure.json"),
    ),
];

#[derive(Debug, Deserialize)]
struct Snapshot {
    description: String,
    args: Vec<String>,
    exit_code: i32,
    #[serde(default)]
    write_file: Option<WriteFile>,
    #[serde(default)]
    message_contains: Vec<String>,
    error: serde_json::Value,
}

#[derive(Debug, Deserialize)]
struct WriteFile {
    name: String,
    content: String,
}

/// Temporary directory removed when the test ends.
struct Scratch(PathBuf);

impl Scratch {
    fn new() -> Self {
        let nanos = SystemTime::now()
            .duration_since(UNIX_EPOCH)
            .expect("the clock must be after the epoch")
            .as_nanos();
        let path = std::env::temp_dir().join(format!(
            "xtalforge-error-snapshots-{}-{nanos}",
            std::process::id()
        ));
        fs::create_dir_all(&path).expect("the snapshot directory must be creatable");
        Self(path)
    }

    fn path(&self) -> &Path {
        &self.0
    }
}

impl Drop for Scratch {
    fn drop(&mut self) {
        let _ = fs::remove_dir_all(&self.0);
    }
}

/// Replaces `a.cif` with the absolute corpus path.
fn resolve_fixtures(args: &mut [String], scratch: &Path) {
    for argument in args {
        if argument == "a.cif" {
            *argument = scratch.join("a.cif").display().to_string();
        }
    }
}

fn run(args: &[String], scratch: &Path) -> std::process::Output {
    Command::new(env!("CARGO_BIN_EXE_xtalforge"))
        .args(args)
        .current_dir(scratch)
        .output()
        .expect("the CLI binary must run")
}

#[test]
fn every_error_snapshot_matches() {
    let scratch = Scratch::new();
    fs::write(scratch.path().join("a.cif"), NACL).expect("NaCl fixture");
    for (name, text) in SNAPSHOTS {
        let snapshot: Snapshot = serde_json::from_str(text)
            .unwrap_or_else(|error| panic!("fixture {name} is not valid JSON: {error}"));
        if let Some(file) = &snapshot.write_file {
            fs::write(scratch.path().join(&file.name), &file.content).unwrap_or_else(|error| {
                panic!("fixture {name} cannot write {}: {error}", file.name)
            });
        }
        let mut args = vec!["--json".to_owned()];
        args.extend(snapshot.args.iter().cloned());
        resolve_fixtures(&mut args, scratch.path());
        let output = run(&args, scratch.path());
        let stdout = String::from_utf8_lossy(&output.stdout);
        let stderr = String::from_utf8_lossy(&output.stderr);
        assert_eq!(
            output.status.code(),
            Some(snapshot.exit_code),
            "fixture {name} ({}) exited with the wrong code\nstdout:\n{stdout}\nstderr:\n{stderr}",
            snapshot.description
        );
        assert!(
            stdout.trim().is_empty(),
            "fixture {name} wrote to stdout: {stdout}"
        );
        assert!(
            !stderr.contains("panicked at"),
            "fixture {name} panicked: {stderr}"
        );
        let value: serde_json::Value = serde_json::from_slice(&output.stderr)
            .unwrap_or_else(|error| panic!("fixture {name} stderr is not JSON: {error}\n{stderr}"));
        let top = value.as_object().expect("the error document is an object");
        let top_keys: BTreeSet<&str> = top.keys().map(String::as_str).collect();
        assert_eq!(
            top_keys,
            BTreeSet::from(["error"]),
            "fixture {name} has unexpected top-level keys"
        );
        let error = value["error"]
            .as_object()
            .expect("the error field is an object");
        let expected = snapshot
            .error
            .as_object()
            .unwrap_or_else(|| panic!("fixture {name} has no error object"));
        let actual_keys: BTreeSet<&str> = error.keys().map(String::as_str).collect();
        let expected_keys: BTreeSet<&str> = expected.keys().map(String::as_str).collect();
        assert_eq!(
            actual_keys, expected_keys,
            "fixture {name} has unexpected error keys"
        );
        assert_eq!(
            error["code"], expected["code"],
            "fixture {name} code drifted"
        );
        assert_eq!(
            error["category"], expected["category"],
            "fixture {name} category drifted"
        );
        let message = error["message"]
            .as_str()
            .expect("the message field is a string");
        assert!(!message.is_empty(), "fixture {name} has an empty message");
        for fragment in &snapshot.message_contains {
            assert!(
                message.contains(fragment),
                "fixture {name} message `{message}` lost `{fragment}`"
            );
        }
    }
}
