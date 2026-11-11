//! CLI test for `powder automatic` and its space-group screen (Task 7.10).

#![allow(
    clippy::unwrap_used,
    clippy::expect_used,
    clippy::panic,
    clippy::indexing_slicing
)]

use std::path::PathBuf;
use std::process::Command;

fn repository() -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR"))
        .parent()
        .and_then(|path| path.parent())
        .expect("the crate lives under crates/")
        .to_path_buf()
}

#[test]
fn powder_automatic_reports_space_group_candidates() {
    let root = repository();
    let structure = root.join("tests/corpus/simple/TiO2.cif");
    let pattern = root.join("tests/golden/phase-7/rietveld/rutile.xy");
    let output = Command::new(env!("CARGO_BIN_EXE_xtalforge"))
        .args([
            "--json",
            "powder",
            "automatic",
            &structure.to_string_lossy(),
            &pattern.to_string_lossy(),
        ])
        .output()
        .expect("the CLI must run");
    assert!(
        output.status.success(),
        "{}",
        String::from_utf8_lossy(&output.stderr)
    );
    let value: serde_json::Value = serde_json::from_slice(&output.stdout).expect("JSON");
    let rows = value["leading_cells"].as_array().expect("a table");
    assert!(!rows.is_empty());
    let groups = rows[0]["space_groups"].as_array().expect("a screen");
    assert!(!groups.is_empty());
    assert!(groups[0]["number"].as_u64().unwrap() > 0);
    assert!(!groups[0]["symbol"].as_str().unwrap().is_empty());
    // `best_cell` is the Pawley-refined cell, not the quick cubic-index edge
    // in `index.cell` (a = 10.26 for rutile).
    let best = value["best_cell"].as_array().expect("a best cell");
    assert_eq!(best.len(), 6);
    let a = best[0].as_f64().unwrap();
    let c = best[2].as_f64().unwrap();
    assert!((a - 4.5993).abs() < 0.01, "{a}");
    assert!((c - 2.9629).abs() < 0.01, "{c}");
    let cubic = value["index"]["cell"]["a"].as_f64().unwrap();
    assert!((cubic - 10.2619).abs() < 0.01, "{cubic}");
}
