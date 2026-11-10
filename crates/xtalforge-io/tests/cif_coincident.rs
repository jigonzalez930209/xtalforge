//! A CIF that lists the whole cell under a non-P1 group (Task 1.6.2).

#![allow(
    clippy::unwrap_used,
    clippy::expect_used,
    clippy::panic,
    clippy::indexing_slicing
)]

use std::path::PathBuf;

use xtalforge_io::read_cif;

const WHOLE_CELL: &str = include_str!("fixtures/whole_cell_under_fm3m.cif");

fn warning(structure: &xtalforge_core::Structure) -> Option<&str> {
    structure
        .metadata()
        .get("coincident_atoms_warning")
        .and_then(serde_json::Value::as_str)
}

#[test]
fn a_whole_cell_under_a_space_group_is_reported() {
    let structure = read_cif(WHOLE_CELL).unwrap();
    let message = warning(&structure).expect("the repeated orbits must be reported");
    assert!(message.contains("whole cell"), "{message}");
}

#[test]
fn the_curated_corpus_has_no_false_positives() {
    let root = PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("../../tests/corpus");
    let mut flagged = Vec::new();
    let mut total = 0;
    let mut stack = vec![root];
    while let Some(directory) = stack.pop() {
        for entry in std::fs::read_dir(&directory).unwrap() {
            let path = entry.unwrap().path();
            if path.is_dir() {
                stack.push(path);
            } else if path.extension().is_some_and(|extension| extension == "cif") {
                let text = std::fs::read_to_string(&path).unwrap();
                if let Ok(structure) = read_cif(&text) {
                    total += 1;
                    if warning(&structure).is_some() {
                        flagged.push(path.file_name().unwrap().to_string_lossy().into_owned());
                    }
                }
            }
        }
    }
    assert!(total > 60, "{total} files read");
    flagged.sort();
    assert!(flagged.is_empty(), "false positives: {flagged:?}");
}
