//! Integration tests for the Phase 6 build commands.

#![allow(
    clippy::unwrap_used,
    clippy::expect_used,
    clippy::panic,
    clippy::indexing_slicing
)]

use std::process::{Command, Output};

fn run(args: &[&str]) -> Output {
    Command::new(env!("CARGO_BIN_EXE_xtalforge"))
        .args(args)
        .output()
        .expect("the CLI binary must run")
}

fn json(output: &Output) -> serde_json::Value {
    assert!(output.status.success(), "{output:?}");
    serde_json::from_slice(&output.stdout).expect("JSON")
}

#[test]
fn fullerenes_and_twists_build() {
    let value = json(&run(&["--json", "build", "c20"]));
    assert_eq!(
        value["structure"]["sites"].as_array().expect("sites").len(),
        20
    );

    let value = json(&run(&[
        "--json", "build", "twist", "2", "2", "0", "--layers", "2",
    ]));
    assert_eq!(
        value["structure"]["sites"].as_array().expect("sites").len(),
        16
    );
}

#[test]
fn ribbons_surfaces_and_polymers_build() {
    let value = json(&run(&[
        "--json", "build", "ribbon", "2", "3", "--vacuum", "12",
    ]));
    assert!(
        !value["structure"]["sites"]
            .as_array()
            .expect("sites")
            .is_empty()
    );

    let value = json(&run(&[
        "--json",
        "build",
        "ztc",
        "--samples",
        "8",
        "--cell",
        "12",
    ]));
    assert!(
        !value["structure"]["sites"]
            .as_array()
            .expect("sites")
            .is_empty()
    );

    let value = json(&run(&[
        "--json",
        "build",
        "schwarzite",
        "--samples",
        "8",
        "--cell",
        "12",
    ]));
    assert!(
        !value["structure"]["sites"]
            .as_array()
            .expect("sites")
            .is_empty()
    );

    let value = json(&run(&["--json", "build", "polymer", "C=C", "3"]));
    assert!(!value["steps"].as_array().expect("steps").is_empty());
}
