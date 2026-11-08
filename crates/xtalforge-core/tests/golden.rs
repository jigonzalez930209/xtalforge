//! Integration tests that read the Phase 1 golden data.

#![allow(
    clippy::unwrap_used,
    clippy::expect_used,
    clippy::panic,
    clippy::indexing_slicing
)]

use std::fs;
use std::path::Path;

use serde_json::Value;
use xtalforge_core::cell_volume;

fn golden_dir() -> std::path::PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR")).join("../../tests/golden/phase-1")
}

#[test]
fn phase1_golden_cells_match_the_cell_volume_api() {
    let simple = golden_dir().join("simple");
    let mut checked = 0;
    for entry in fs::read_dir(&simple).expect("golden simple directory must exist") {
        let path = entry.expect("directory entry").path();
        if path.extension().and_then(|extension| extension.to_str()) != Some("json") {
            continue;
        }
        let data: Value =
            serde_json::from_str(&fs::read_to_string(&path).expect("golden file must be readable"))
                .expect("golden file must be valid JSON");
        let parameters = data["cell"]["parameters"]
            .as_array()
            .expect("cell parameters must be an array");
        let expected = data["cell"]["volume"]
            .as_f64()
            .expect("cell volume must be a number");
        let value = |index: usize| {
            parameters[index]
                .as_f64()
                .expect("cell parameter must be a number")
        };
        let got = cell_volume(value(0), value(1), value(2), value(3), value(4), value(5))
            .expect("golden cell must be valid");
        assert!(
            (got - expected).abs() <= 1e-8 * expected.abs().max(1.0),
            "{}: volume {got} != {expected}",
            path.display()
        );
        assert!(
            data["p1"]["n_atoms"].as_u64().expect("P1 atom count")
                >= data["n_sites"].as_u64().expect("site count")
        );
        checked += 1;
    }
    assert!(checked >= 10, "expected at least ten simple golden files");
}

#[test]
fn malformed_golden_case_records_a_normalized_error() {
    let path = golden_dir().join("edge/malformed.json");
    let data: Value = serde_json::from_str(
        &fs::read_to_string(path).expect("malformed golden file must be readable"),
    )
    .expect("malformed golden file must be valid JSON");
    let message = data["error"]["message"]
        .as_str()
        .expect("an error message must be recorded");
    assert!(message.contains("edge/malformed.cif"));
    assert!(!message.contains("/home/"));
}
