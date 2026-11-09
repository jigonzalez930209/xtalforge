//! Golden ring comparison (Task 1.12).
//!
//! Counts must match the Crystal-Builder Franzblau search for every structure
//! with ring data. When the golden list is complete, the periodic paths are
//! compared after translating each ring home and forgetting its direction.

#![allow(
    clippy::unwrap_used,
    clippy::expect_used,
    clippy::panic,
    clippy::indexing_slicing
)]

use std::collections::BTreeSet;
use std::fs;
use std::path::{Path, PathBuf};

use serde_json::Value;
use xtalforge_core::{Ring, rings};
use xtalforge_io::read_cif;

#[test]
fn golden_ring_counts_and_paths_match() {
    let root = Path::new(env!("CARGO_MANIFEST_DIR")).join("../../tests");
    let mut checked = 0;
    let mut full_paths = 0;
    for path in golden_files(&root.join("golden/phase-1")) {
        let golden: Value = serde_json::from_str(&fs::read_to_string(&path).unwrap()).unwrap();
        let Some(source) = golden["source"].as_str() else {
            continue;
        };
        let Some(expected) = golden.get("rings") else {
            continue;
        };
        if expected.is_null() {
            continue;
        }
        let structure =
            read_cif(&fs::read_to_string(root.join("corpus").join(source)).unwrap()).unwrap();
        if let Some(error) = expected.get("error") {
            assert!(
                rings(&structure, 8).is_err(),
                "{}: golden {} should also be refused",
                path.display(),
                error["type"].as_str().unwrap_or("error")
            );
            checked += 1;
            continue;
        }
        let count = expected["count"].as_u64().unwrap() as usize;
        let found = rings(&structure, 8).unwrap();
        assert_eq!(found.len(), count, "{}: ring count", path.display());
        let truncated = expected["truncated"].as_bool().unwrap_or(false);
        if !truncated {
            let p1 = structure.p1().unwrap();
            let mapping = map_atoms(&golden, &p1);
            assert!(
                mapping.iter().all(Option::is_some),
                "{}: every P1 atom must map",
                path.display()
            );
            let expected_paths = expected["paths"]
                .as_array()
                .unwrap()
                .iter()
                .map(golden_key)
                .collect::<BTreeSet<_>>();
            let found_paths = found
                .iter()
                .map(|ring| ring_key(ring, &mapping))
                .collect::<BTreeSet<_>>();
            assert_eq!(
                found_paths,
                expected_paths,
                "{}: ring paths",
                path.display()
            );
            full_paths += 1;
        }
        checked += 1;
    }
    assert!(checked >= 40, "expected the ring corpus, got {checked}");
    assert!(
        full_paths >= 30,
        "expected complete golden paths, got {full_paths}"
    );
}

type Step = (usize, [i32; 3]);
type RingKey = BTreeSet<Step>;

fn ring_key(ring: &Ring, mapping: &[Option<usize>]) -> RingKey {
    normalize(
        ring.path
            .iter()
            .map(|step| (mapping[step.atom].unwrap(), step.image))
            .collect::<Vec<_>>(),
    )
}

fn golden_key(path: &Value) -> RingKey {
    let steps = path
        .as_array()
        .unwrap()
        .iter()
        .map(|step| {
            let atom = step[0].as_u64().unwrap() as usize;
            let image = [
                step[1][0].as_i64().unwrap() as i32,
                step[1][1].as_i64().unwrap() as i32,
                step[1][2].as_i64().unwrap() as i32,
            ];
            (atom, image)
        })
        .collect::<Vec<_>>();
    normalize(steps)
}

fn normalize(steps: Vec<Step>) -> RingKey {
    let minimum = steps.iter().min().copied().unwrap_or((0, [0, 0, 0]));
    steps
        .into_iter()
        .map(|(atom, image)| {
            (
                atom,
                [
                    image[0] - minimum.1[0],
                    image[1] - minimum.1[1],
                    image[2] - minimum.1[2],
                ],
            )
        })
        .collect()
}

fn golden_files(root: &Path) -> Vec<PathBuf> {
    let mut files = Vec::new();
    let mut pending = vec![root.to_path_buf()];
    while let Some(directory) = pending.pop() {
        for entry in fs::read_dir(directory).unwrap() {
            let path = entry.unwrap().path();
            if path.is_dir() {
                pending.push(path);
            } else if path.extension().and_then(|value| value.to_str()) == Some("json") {
                files.push(path);
            }
        }
    }
    files.sort();
    files
}

fn map_atoms(golden: &Value, p1: &xtalforge_core::P1Cell) -> Vec<Option<usize>> {
    let elements = golden["p1"]["elements"].as_array().unwrap();
    let fracs = golden["p1"]["frac"].as_array().unwrap();
    let mut used = vec![false; elements.len()];
    let mut mapping = Vec::with_capacity(p1.atoms.len());
    for atom in &p1.atoms {
        let mut best: Option<(f64, usize)> = None;
        for (index, element) in elements.iter().enumerate() {
            if used[index] || element.as_str() != Some(atom.element.as_str()) {
                continue;
            }
            let frac = &fracs[index];
            let target = [
                frac[0].as_f64().unwrap(),
                frac[1].as_f64().unwrap(),
                frac[2].as_f64().unwrap(),
            ];
            let distance = p1.lattice.distance(atom.frac, target);
            if best.is_none_or(|(value, _)| distance < value) {
                best = Some((distance, index));
            }
        }
        match best {
            Some((distance, index)) if distance < 0.05 => {
                used[index] = true;
                mapping.push(Some(index));
            }
            _ => mapping.push(None),
        }
    }
    mapping
}
