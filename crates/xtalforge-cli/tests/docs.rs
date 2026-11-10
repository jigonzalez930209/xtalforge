//! Completions and the manual page (Task 8.5.4).

#![allow(clippy::unwrap_used, clippy::expect_used, clippy::panic)]

use std::process::Command;

fn run(args: &[&str]) -> std::process::Output {
    Command::new(env!("CARGO_BIN_EXE_xtalforge"))
        .args(args)
        .output()
        .expect("the CLI binary must run")
}

#[test]
fn every_shell_completion_is_generated() {
    for shell in ["bash", "zsh", "fish", "powershell"] {
        let output = run(&["completions", shell]);
        assert!(output.status.success(), "{shell}");
        let text = String::from_utf8(output.stdout).unwrap();
        assert!(text.contains("xtalforge"), "{shell} output");
        assert!(text.len() > 100, "{shell} output is too short");
    }
}

#[test]
fn the_manual_page_has_a_header() {
    let output = run(&["man"]);
    assert!(output.status.success());
    let text = String::from_utf8(output.stdout).unwrap();
    assert!(text.contains(".TH"), "roff header");
    assert!(text.contains("xtalforge"));
}
