//! Wyckoff generator tests (Task 20.1.1).

#![allow(
    clippy::unwrap_used,
    clippy::expect_used,
    clippy::panic,
    clippy::indexing_slicing
)]

use xtalforge_core::StructureDocument;
use xtalforge_generation::{
    CrystalOptions, generate_crystal, generate_crystal_range, generation_statistics,
};
use xtalforge_symmetry::{SymmetryOptions, detect_symmetry};

fn nacl(seed: u64) -> CrystalOptions {
    CrystalOptions {
        space_group: 225,
        composition: vec![("Na".to_owned(), 4), ("Cl".to_owned(), 4)],
        volume_per_atom: 22.0,
        min_distance: 1.2,
        max_attempts: 200,
        seed,
    }
}

#[test]
fn generates_rock_salt_and_redetects_the_group() {
    let crystal = generate_crystal(&nacl(7)).expect("crystal");
    assert_eq!(crystal.space_group, 225);
    assert_eq!(crystal.structure.sites().len(), 2);
    let p1 = crystal.structure.p1().expect("p1");
    assert_eq!(p1.atoms.len(), 8);
    let dataset = detect_symmetry(&crystal.structure, SymmetryOptions::default()).expect("detect");
    assert_eq!(dataset.number, 225);
    let sodium = p1.atoms.iter().filter(|atom| atom.element == "Na").count();
    let chlorine = p1.atoms.iter().filter(|atom| atom.element == "Cl").count();
    assert_eq!(sodium, 4);
    assert_eq!(chlorine, 4);
}

#[test]
fn same_seed_gives_identical_structure() {
    let first = generate_crystal(&nacl(11)).expect("first");
    let second = generate_crystal(&nacl(11)).expect("second");
    assert_eq!(first.structure, second.structure);
    assert_eq!(first.wyckoffs, second.wyckoffs);
    assert_eq!(first.attempt, second.attempt);
}

#[test]
fn different_seeds_differ_for_a_free_parameter() {
    let options = |seed| CrystalOptions {
        space_group: 186,
        composition: vec![("Zn".to_owned(), 2), ("O".to_owned(), 2)],
        volume_per_atom: 24.0,
        min_distance: 1.0,
        max_attempts: 200,
        seed,
    };
    let first = generate_crystal(&options(1)).expect("first");
    let second = generate_crystal(&options(2)).expect("second");
    assert_ne!(first.structure.sites(), second.structure.sites());
}

#[test]
fn hexagonal_metal_is_generated() {
    let options = CrystalOptions {
        space_group: 194,
        composition: vec![("Mg".to_owned(), 2)],
        volume_per_atom: 22.0,
        min_distance: 1.5,
        max_attempts: 200,
        seed: 3,
    };
    let crystal = generate_crystal(&options).expect("hexagonal");
    let dataset = detect_symmetry(&crystal.structure, SymmetryOptions::default()).expect("detect");
    assert_eq!(dataset.number, 194);
}

#[test]
fn range_reports_every_group() {
    let options = nacl(5);
    let groups = [216u32, 225];
    let results = generate_crystal_range(&groups, &options).expect("range");
    assert_eq!(results.len(), groups.len());
    for (group, result) in groups.iter().zip(results) {
        let crystal = result.expect("group");
        assert_eq!(crystal.space_group, *group);
    }
}

#[test]
fn statistics_count_successes() {
    let stats = generation_statistics(&nacl(100), 4).expect("stats");
    assert_eq!(stats.seeds, 4);
    assert_eq!(stats.succeeded + stats.failed, 4);
    assert!(
        stats.succeeded >= 3,
        "expected most seeds to work: {stats:?}"
    );
}

#[test]
fn invalid_options_are_rejected() {
    let mut options = nacl(1);
    options.space_group = 231;
    assert!(generate_crystal(&options).is_err());
    options = nacl(1);
    options.composition.clear();
    assert!(generate_crystal(&options).is_err());
    options = nacl(1);
    options.volume_per_atom = 0.0;
    assert!(generate_crystal(&options).is_err());
    options = nacl(1);
    options.max_attempts = 0;
    assert!(generate_crystal(&options).is_err());
    options = nacl(1);
    options.composition = vec![("Zz".to_owned(), 4)];
    assert!(generate_crystal(&options).is_err());
}

#[test]
fn json_round_trip() {
    let crystal = generate_crystal(&nacl(9)).expect("crystal");
    let crystal_text = crystal.to_json().expect("json");
    let document = StructureDocument::new(crystal.structure.clone());
    let text = document.to_json().expect("document json");
    let parsed = StructureDocument::from_json(&text).expect("parse");
    assert_eq!(parsed.structure, crystal.structure);
    assert!(crystal_text.contains("\"space_group\":225"));
}
