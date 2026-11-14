//! Layer and rod generator tests (Task 20.1.3).

#![allow(
    clippy::unwrap_used,
    clippy::expect_used,
    clippy::panic,
    clippy::indexing_slicing
)]

use xtalforge_generation::{
    LayerCrystalOptions, RodCrystalOptions, generate_layer_crystal, generate_rod_crystal,
    layer_groups, rod_groups,
};
use xtalforge_symmetry::{SymmetryOptions, detect_layer_symmetry};

#[test]
fn layer_catalogue_has_eighty_groups() {
    let groups = layer_groups().expect("groups");
    assert_eq!(groups.len(), 80);
    assert_eq!(groups[0].number, 1);
    assert_eq!(groups[79].number, 80);
}

#[test]
fn generates_a_layer_in_p1() {
    let options = LayerCrystalOptions::new(1, vec![("C".to_owned(), 1), ("N".to_owned(), 1)], 7);
    let crystal = generate_layer_crystal(&options).expect("layer");
    assert_eq!(crystal.layer_group, 1);
    let dataset = detect_layer_symmetry(&crystal.structure, SymmetryOptions::default().symprec)
        .expect("detect");
    assert_eq!(dataset.number, 1);
    let c = crystal.structure.lattice().parameters()[2];
    assert!((c - xtalforge_generation::LAYER_VACUUM).abs() < 1e-9);
}

#[test]
fn generates_an_inversion_layer() {
    let options = LayerCrystalOptions::new(2, vec![("C".to_owned(), 2)], 3);
    let crystal = generate_layer_crystal(&options).expect("layer");
    let p1 = crystal.structure.p1().expect("p1");
    assert_eq!(p1.atoms.len(), 2);
    let dataset = detect_layer_symmetry(&crystal.structure, SymmetryOptions::default().symprec)
        .expect("detect");
    assert_eq!(dataset.number, 2);
}

#[test]
fn layer_and_rod_options_are_validated() {
    let mut options = LayerCrystalOptions::new(0, vec![("C".to_owned(), 1)], 1);
    assert!(generate_layer_crystal(&options).is_err());
    options = LayerCrystalOptions::new(1, Vec::new(), 1);
    assert!(generate_layer_crystal(&options).is_err());
    let mut rod = RodCrystalOptions::new(0, vec![("C".to_owned(), 1)], 1);
    assert!(generate_rod_crystal(&rod).is_err());
    rod = RodCrystalOptions::new(1, vec![("C".to_owned(), 1)], 1);
    rod.repeat = 0.0;
    assert!(generate_rod_crystal(&rod).is_err());
}

#[test]
fn rod_catalogue_has_fourteen_groups() {
    let groups = rod_groups().expect("groups");
    assert_eq!(groups.len(), 14);
    assert!(groups.iter().any(|group| group.symbol == "p6mm"));
}

#[test]
fn generates_an_inversion_rod() {
    let options = RodCrystalOptions {
        rod_group: 2,
        composition: vec![("C".to_owned(), 2)],
        repeat: 4.0,
        min_distance: 0.5,
        max_attempts: 64,
        seed: 5,
    };
    let crystal = generate_rod_crystal(&options).expect("rod");
    let p1 = crystal.structure.p1().expect("p1");
    assert!(!p1.atoms.is_empty());
    assert!(p1.atoms.len().is_multiple_of(2));
}
