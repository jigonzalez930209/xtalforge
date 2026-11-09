//! Golden bond-set and coordination comparison (Task 1.11.5).
//!
//! The golden P1 cells come from Crystal-Builder and are not required to use
//! the same atom order as ours; atoms are matched by element and periodic
//! distance before the bonds are compared. Six prepared structures above the
//! golden generator's atom limit have no bond list and are skipped.

#![allow(
    clippy::unwrap_used,
    clippy::expect_used,
    clippy::panic,
    clippy::indexing_slicing
)]

use std::collections::BTreeMap;
use std::fs;
use std::path::{Path, PathBuf};

use serde_json::Value;
use xtalforge_core::P1Cell;
use xtalforge_io::read_cif;

/// P1 positions may differ because the golden side standardizes some sites
/// (MIL-100 O27 is 0.021 Angstrom away). Bond lengths inherit that shift.
const POSITION_TOLERANCE: f64 = 0.05;
const DISTANCE_TOLERANCE: f64 = 0.05;
const GOLDEN_WITHOUT_BONDS: usize = 6;

type BondMap = BTreeMap<(usize, usize, [i32; 3]), f64>;

#[test]
fn full_corpus_bond_sets_and_coordination_match_the_golden_data() {
    let root = Path::new(env!("CARGO_MANIFEST_DIR")).join("../../tests");
    let mut checked = 0;
    let mut without_bonds = 0;
    for path in golden_files(&root.join("golden/phase-1")) {
        let golden: Value = serde_json::from_str(&fs::read_to_string(&path).unwrap()).unwrap();
        if golden.get("error").is_some() {
            continue;
        }
        let Some(source) = golden["source"].as_str() else {
            continue;
        };
        let corpus = root.join("corpus").join(source);
        let structure = read_cif(&fs::read_to_string(&corpus).unwrap()).unwrap();
        let p1 = structure.p1().unwrap();
        let expected_atoms = golden["p1"]["n_atoms"].as_u64().unwrap() as usize;
        assert_eq!(
            p1.n_atoms(),
            expected_atoms,
            "{}: P1 atom count",
            path.display()
        );
        let Some(entries) = golden["bonds"].as_array() else {
            without_bonds += 1;
            continue;
        };

        let mapping = map_atoms(&golden, &p1);
        assert!(
            mapping.iter().all(Option::is_some),
            "{}: every P1 atom must map to a golden atom",
            path.display()
        );
        let expected = expected_bonds(entries);
        let (found, found_coordination) = found_bonds(&structure, &mapping);
        if found.len() != expected.len() {
            let extra: Vec<_> = found
                .keys()
                .filter(|key| !expected.contains_key(*key))
                .collect();
            let missing: Vec<_> = expected
                .keys()
                .filter(|key| !found.contains_key(*key))
                .collect();
            eprintln!("{} extra {:?} missing {:?}", path.display(), extra, missing);
        }
        if found.len() != expected.len() {
            let extra: Vec<_> = found
                .keys()
                .filter(|key| !expected.contains_key(*key))
                .collect();
            eprintln!("{} extra {:?}", path.display(), extra);
        }
        assert_eq!(
            found.len(),
            expected.len(),
            "{}: bond count",
            path.display()
        );
        for (key, distance) in &expected {
            let got = found.get(key).unwrap_or_else(|| {
                panic!("{}: bond {key:?} is missing", path.display());
            });
            assert!(
                (got - distance).abs() <= DISTANCE_TOLERANCE,
                "{}: bond {key:?} length {got} != {distance}",
                path.display()
            );
        }

        let mut expected_coordination = vec![0_usize; p1.n_atoms()];
        for (i, j, _) in expected.keys() {
            expected_coordination[*i] += 1;
            expected_coordination[*j] += 1;
        }
        assert_eq!(
            found_coordination,
            expected_coordination,
            "{}: coordination numbers",
            path.display()
        );
        checked += 1;
    }
    assert!(
        checked >= 60,
        "expected at least 60 golden structures, got {checked}"
    );
    assert_eq!(without_bonds, GOLDEN_WITHOUT_BONDS);
}

fn expected_bonds(entries: &[Value]) -> BTreeMap<(usize, usize, [i32; 3]), f64> {
    let mut bonds = BTreeMap::new();
    for bond in entries {
        let i = bond["i"].as_u64().unwrap() as usize;
        let j = bond["j"].as_u64().unwrap() as usize;
        let image = [
            bond["image"][0].as_i64().unwrap() as i32,
            bond["image"][1].as_i64().unwrap() as i32,
            bond["image"][2].as_i64().unwrap() as i32,
        ];
        let (i, j, image) = normalize(i, j, image);
        bonds.insert((i, j, image), bond["distance"].as_f64().unwrap());
    }
    bonds
}

fn found_bonds(
    structure: &xtalforge_core::Structure,
    mapping: &[Option<usize>],
) -> (BondMap, Vec<usize>) {
    let graph = structure.bond_graph().unwrap();
    let mut bonds = BTreeMap::new();
    for bond in &graph.bonds {
        let (Some(Some(first)), Some(Some(second))) = (mapping.get(bond.i), mapping.get(bond.j))
        else {
            continue;
        };
        let (i, j, image) = normalize(*first, *second, bond.image);
        bonds.insert((i, j, image), bond.distance);
    }
    let coordination = (0..graph.n_atoms())
        .map(|index| graph.degree(index))
        .collect();
    (bonds, coordination)
}

/// One-to-one match by element and periodic distance.
fn map_atoms(golden: &Value, p1: &P1Cell) -> Vec<Option<usize>> {
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
            Some((distance, index)) if distance < POSITION_TOLERANCE => {
                used[index] = true;
                mapping.push(Some(index));
            }
            _ => mapping.push(None),
        }
    }
    mapping
}

/// Bond identity: the smaller of `(i, j, image)` and `(j, i, -image)`.
fn normalize(first: usize, second: usize, image: [i32; 3]) -> (usize, usize, [i32; 3]) {
    let flipped = (second, first, [-image[0], -image[1], -image[2]]);
    if (second, first, flipped.2) < (first, second, image) {
        flipped
    } else {
        (first, second, image)
    }
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
