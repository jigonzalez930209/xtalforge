//! CLI test for `powder index --general` (Task 7.5).

#![allow(
    clippy::unwrap_used,
    clippy::expect_used,
    clippy::panic,
    clippy::indexing_slicing
)]

use std::process::Command;

fn run(args: &[&str]) -> std::process::Output {
    Command::new(env!("CARGO_BIN_EXE_xtalforge"))
        .args(args)
        .output()
        .expect("the CLI must run")
}

/// TiO2 rutile (a = 4.5937, c = 2.9587) at Cu K-alpha, as a sampled pattern.
fn rutile_pattern() -> String {
    let centers = [
        19.32, 27.46, 30.21, 36.11, 39.22, 41.27, 44.08, 54.37, 56.68, 62.81, 69.07, 72.49,
    ];
    let mut text = String::new();
    for step in 0..3000 {
        let two_theta = 15.0 + step as f64 * 0.02;
        let mut intensity = 0.0;
        for (index, center) in centers.iter().enumerate() {
            let weight = 1.0 / (1.0 + index as f64);
            let delta = (two_theta - center) / 0.12;
            intensity += weight * (-delta * delta).exp();
        }
        text.push_str(&format!("{two_theta:.3} {intensity:.6}\n"));
    }
    text
}

#[test]
fn powder_index_reports_general_candidates() {
    let directory =
        std::env::temp_dir().join(format!("xtalforge-powder-{}-{}", std::process::id(), 1));
    std::fs::create_dir_all(&directory).unwrap();
    let path = directory.join("rutile.xy");
    std::fs::write(&path, rutile_pattern()).unwrap();
    let output = run(&[
        "--json",
        "powder",
        "index",
        &path.to_string_lossy(),
        "--general",
        "--max-index",
        "3",
        "--tolerance",
        "0.05",
    ]);
    assert!(
        output.status.success(),
        "{}",
        String::from_utf8_lossy(&output.stderr)
    );
    let value: serde_json::Value = serde_json::from_slice(&output.stdout).expect("JSON");
    let candidates = value["candidates"]["candidates"]
        .as_array()
        .expect("candidate list");
    assert!(!candidates.is_empty());
    let best = &candidates[0];
    assert_eq!(best["system"], "tetragonal");
    assert!((best["cell"]["a"].as_f64().unwrap() - 4.5937).abs() < 0.05);
    assert!((best["cell"]["c"].as_f64().unwrap() - 2.9587).abs() < 0.05);
    assert!(best["m20"].as_f64().unwrap() > 0.0);
    std::fs::remove_dir_all(&directory).unwrap();
}
