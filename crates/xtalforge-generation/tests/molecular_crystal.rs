//! Molecular crystal tests (Task 20.1.2).

#![allow(
    clippy::unwrap_used,
    clippy::expect_used,
    clippy::panic,
    clippy::indexing_slicing
)]

use xtalforge_core::{Lattice, Site, Structure, StructureBuilder};
use xtalforge_generation::{MolecularCrystalOptions, generate_molecular_crystal};
use xtalforge_symmetry::{SymmetryOptions, detect_symmetry};

fn water() -> Structure {
    let lattice = Lattice::from_parameters(10.0, 10.0, 10.0, 90.0, 90.0, 90.0).expect("lattice");
    StructureBuilder::new("water", lattice)
        .site(Site::new("O1", "O", [0.5, 0.5, 0.5]).expect("O"))
        .site(Site::new("H2", "H", [0.596, 0.5, 0.5]).expect("H"))
        .site(Site::new("H3", "H", [0.476, 0.593, 0.5]).expect("H"))
        .build()
        .expect("water")
}

#[test]
fn places_a_water_molecule_in_p1() {
    let options = MolecularCrystalOptions {
        space_group: 1,
        molecule: water(),
        copies: 1,
        volume_per_atom: 20.0,
        min_distance: 0.5,
        tolerance: 0.25,
        max_attempts: 32,
        seed: 4,
    };
    let crystal = generate_molecular_crystal(&options).expect("crystal");
    let p1 = crystal.structure.p1().expect("p1");
    assert_eq!(p1.atoms.len(), 3);
    assert_eq!(crystal.molecules_per_cell, 1);
    assert_eq!(crystal.copies, 1);
}

#[test]
fn places_a_molecule_on_a_general_position_of_p_minus_1() {
    let options = MolecularCrystalOptions {
        space_group: 2,
        molecule: water(),
        copies: 1,
        volume_per_atom: 20.0,
        min_distance: 0.5,
        tolerance: 0.25,
        max_attempts: 64,
        seed: 9,
    };
    let crystal = generate_molecular_crystal(&options).expect("crystal");
    let p1 = crystal.structure.p1().expect("p1");
    assert_eq!(p1.atoms.len(), 6);
    assert_eq!(crystal.molecules_per_cell, 2);
    let dataset = detect_symmetry(&crystal.structure, SymmetryOptions::default()).expect("detect");
    assert_eq!(dataset.number, 2);
}

#[test]
fn same_seed_is_deterministic() {
    let options = MolecularCrystalOptions {
        space_group: 2,
        molecule: water(),
        copies: 1,
        volume_per_atom: 20.0,
        min_distance: 0.5,
        tolerance: 0.25,
        max_attempts: 64,
        seed: 12,
    };
    let first = generate_molecular_crystal(&options).expect("first");
    let second = generate_molecular_crystal(&options).expect("second");
    assert_eq!(first.structure, second.structure);
}

#[test]
fn empty_molecule_and_invalid_options_fail() {
    let lattice = Lattice::from_parameters(10.0, 10.0, 10.0, 90.0, 90.0, 90.0).expect("lattice");
    let empty = StructureBuilder::new("empty", lattice)
        .build()
        .expect("empty structure");
    let mut options = MolecularCrystalOptions::new(1, empty, 1);
    assert!(generate_molecular_crystal(&options).is_err());
    options = MolecularCrystalOptions::new(1, water(), 1);
    options.copies = 0;
    assert!(generate_molecular_crystal(&options).is_err());
    options = MolecularCrystalOptions::new(1, water(), 1);
    options.tolerance = 0.0;
    assert!(generate_molecular_crystal(&options).is_err());
}
