//! CLI test for `build carbon` (Task 6.8).

#![allow(
    clippy::unwrap_used,
    clippy::expect_used,
    clippy::panic,
    clippy::indexing_slicing
)]

use std::process::Command;

#[test]
fn build_carbon_reports_a_percolating_recipe() {
    let directory =
        std::env::temp_dir().join(format!("xtalforge-carbon-{}-{}", std::process::id(), 1));
    std::fs::create_dir_all(&directory).unwrap();
    let path = directory.join("carbon.cif");
    let output = Command::new(env!("CARGO_BIN_EXE_xtalforge"))
        .args([
            "build",
            "carbon",
            "--repeat-x",
            "1",
            "--repeat-y",
            "1",
            "--repeat-z",
            "1",
            "--seed",
            "0",
            "--out",
            &path.to_string_lossy(),
        ])
        .output()
        .expect("the CLI must run");
    assert!(
        output.status.success(),
        "{}",
        String::from_utf8_lossy(&output.stderr)
    );
    let text = std::fs::read_to_string(&path).unwrap();
    assert!(text.contains("_cell_length_a"));
    assert!(text.lines().count() > 300);
    std::fs::remove_dir_all(&directory).unwrap();
}
