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

#[test]
fn powder_pawley_recovers_the_golden_pattern() {
    let golden: serde_json::Value = serde_json::from_str(include_str!(
        "../../../tests/golden/phase-7/pawley/one-atom-cubic.json"
    ))
    .unwrap();
    let cell = golden["cell"].as_array().unwrap();
    let value = |index: usize| cell[index].as_f64().unwrap();
    let directory =
        std::env::temp_dir().join(format!("xtalforge-pawley-{}-{}", std::process::id(), 2));
    std::fs::create_dir_all(&directory).unwrap();
    let structure_path = directory.join("one.cif");
    std::fs::write(
        &structure_path,
        format!(
            "data_one\n_cell_length_a {}\n_cell_length_b {}\n_cell_length_c {}\n\
             _cell_angle_alpha {}\n_cell_angle_beta {}\n_cell_angle_gamma {}\n\
             loop_\n_atom_site_label\n_atom_site_type_symbol\n_atom_site_fract_x\n\
             _atom_site_fract_y\n_atom_site_fract_z\n_atom_site_occupancy\nC1 C 0 0 0 1\n",
            value(0),
            value(1),
            value(2),
            value(3),
            value(4),
            value(5)
        ),
    )
    .unwrap();
    let pattern_path = directory.join("one.xy");
    let mut pattern = String::new();
    for peak in golden["peaks"].as_array().unwrap() {
        pattern.push_str(&format!(
            "{} {}\n",
            peak[0].as_f64().unwrap(),
            peak[1].as_f64().unwrap()
        ));
    }
    std::fs::write(&pattern_path, pattern).unwrap();
    let output = run(&[
        "--json",
        "powder",
        "pawley",
        &structure_path.to_string_lossy(),
        &pattern_path.to_string_lossy(),
        "--profile",
        "gaussian",
    ]);
    assert!(
        output.status.success(),
        "{}",
        String::from_utf8_lossy(&output.stderr)
    );
    let value: serde_json::Value = serde_json::from_slice(&output.stdout).expect("JSON");
    assert_eq!(value["reflections"].as_array().unwrap().len(), 80);
    assert!((value["cell"][0].as_f64().unwrap() - 4.0).abs() < 5e-3);
    assert!((value["zero_shift"].as_f64().unwrap() - 0.02).abs() < 5e-3);
    assert!(value["rwp"].as_f64().unwrap() < 0.02);
    std::fs::remove_dir_all(&directory).unwrap();
}
