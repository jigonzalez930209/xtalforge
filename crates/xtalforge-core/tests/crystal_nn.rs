//! Integration tests for the CrystalNN-like coordination model (Task 17.3.4).
//!
//! The expected coordination numbers follow the pymatgen `CrystalNN` reference
//! behaviour for the same idealised structures. No external package is used.

#![allow(
    clippy::unwrap_used,
    clippy::expect_used,
    clippy::panic,
    clippy::indexing_slicing
)]

use xtalforge_core::{
    CoordinationShell, CoordinationStrategy, Lattice, Site, Structure, StructureBuilder,
    coordination,
};

fn build(name: &str, lattice: Lattice, sites: &[(&str, &str, [f64; 3])]) -> Structure {
    let mut builder = StructureBuilder::new(name, lattice);
    for (label, element, frac) in sites {
        builder = builder.site(Site::new(*label, element, *frac).expect("valid site"));
    }
    builder.build().expect("valid structure")
}

fn crystal_nn(structure: &Structure, index: usize, cutoff: f64) -> CoordinationShell {
    coordination(structure, index, CoordinationStrategy::CrystalNn { cutoff })
        .expect("CrystalNN coordination")
}

fn elements_of(structure: &Structure, shell: &CoordinationShell) -> Vec<String> {
    let p1 = structure.p1().expect("P1 expansion");
    shell
        .neighbors
        .iter()
        .map(|neighbor| {
            p1.atoms
                .get(neighbor.index)
                .expect("neighbor index")
                .element
                .clone()
        })
        .collect()
}

fn assert_single_element(structure: &Structure, shell: &CoordinationShell, element: &str) {
    let elements = elements_of(structure, shell);
    assert_eq!(elements.len(), shell.number);
    assert!(
        elements.iter().all(|found| found == element),
        "expected only {element}, got {elements:?}"
    );
}

fn assert_distances(shell: &CoordinationShell, expected: f64, tolerance: f64) {
    for neighbor in &shell.neighbors {
        assert!(
            (neighbor.distance - expected).abs() < tolerance,
            "expected {expected}, got {}",
            neighbor.distance
        );
        assert!(neighbor.weight > 0.0);
    }
}

fn p1_elements(structure: &Structure) -> Vec<String> {
    structure
        .p1()
        .expect("P1 expansion")
        .atoms
        .iter()
        .map(|atom| atom.element.clone())
        .collect()
}

fn assert_equal_weights(shell: &CoordinationShell, tolerance: f64) {
    let Some(first) = shell.neighbors.first() else {
        return;
    };
    for neighbor in &shell.neighbors {
        assert!(
            (neighbor.weight - first.weight).abs() < tolerance,
            "weights differ: {} vs {}",
            neighbor.weight,
            first.weight
        );
    }
}

fn rocksalt() -> Structure {
    build(
        "NaCl",
        Lattice::cubic(5.64).unwrap(),
        &[
            ("Na1", "Na", [0.0, 0.0, 0.0]),
            ("Na2", "Na", [0.0, 0.5, 0.5]),
            ("Na3", "Na", [0.5, 0.0, 0.5]),
            ("Na4", "Na", [0.5, 0.5, 0.0]),
            ("Cl1", "Cl", [0.5, 0.0, 0.0]),
            ("Cl2", "Cl", [0.0, 0.5, 0.0]),
            ("Cl3", "Cl", [0.0, 0.0, 0.5]),
            ("Cl4", "Cl", [0.5, 0.5, 0.5]),
        ],
    )
}

fn diamond() -> Structure {
    build(
        "C",
        Lattice::cubic(3.567).unwrap(),
        &[
            ("C1", "C", [0.0, 0.0, 0.0]),
            ("C2", "C", [0.0, 0.5, 0.5]),
            ("C3", "C", [0.5, 0.0, 0.5]),
            ("C4", "C", [0.5, 0.5, 0.0]),
            ("C5", "C", [0.25, 0.25, 0.25]),
            ("C6", "C", [0.25, 0.75, 0.75]),
            ("C7", "C", [0.75, 0.25, 0.75]),
            ("C8", "C", [0.75, 0.75, 0.25]),
        ],
    )
}

fn cesium_chloride() -> Structure {
    build(
        "CsCl",
        Lattice::cubic(4.11).unwrap(),
        &[
            ("Cl1", "Cl", [0.0, 0.0, 0.0]),
            ("Cs1", "Cs", [0.5, 0.5, 0.5]),
        ],
    )
}

fn bcc() -> Structure {
    build(
        "bcc",
        Lattice::cubic(4.0).unwrap(),
        &[
            ("Fe1", "Fe", [0.0, 0.0, 0.0]),
            ("Fe2", "Fe", [0.5, 0.5, 0.5]),
        ],
    )
}

fn polonium() -> Structure {
    let a = 3.35;
    build(
        "Po",
        Lattice::orthorhombic(2.0 * a, a, a).unwrap(),
        &[
            ("Po1", "Po", [0.0, 0.0, 0.0]),
            ("Po2", "Po", [0.5, 0.0, 0.0]),
        ],
    )
}

/// A synthetic L1_2 cell with Na at the corners and Cl at one face centre.
///
/// This is the Na3Cl arrangement from Cl corners and Na face centres,
/// translated so that a Na sits at the origin. Each Na then has four Cl and
/// eight Na neighbours at the same distance, which isolates the cation-anion
/// preference from the distance and solid-angle terms. The two species have a
/// large Pauling difference (Na 0.93, Cl 3.16).
fn sodium_chloride_l12() -> Structure {
    build(
        "Na3Cl",
        Lattice::cubic(4.0).unwrap(),
        &[
            ("Na1", "Na", [0.0, 0.0, 0.0]),
            ("Na2", "Na", [0.0, 0.5, 0.5]),
            ("Na3", "Na", [0.5, 0.0, 0.5]),
            ("Cl1", "Cl", [0.5, 0.5, 0.0]),
        ],
    )
}

#[test]
fn rock_salt_has_six_coordination() {
    let structure = rocksalt();
    let shell = crystal_nn(&structure, 0, 4.5);
    assert_eq!(shell.number, 6);
    assert_single_element(&structure, &shell, "Cl");
    assert_distances(&shell, 2.82, 1e-6);
    assert_equal_weights(&shell, 1e-9);
}

#[test]
fn diamond_has_four_coordination() {
    let structure = diamond();
    let shell = crystal_nn(&structure, 0, 1.7);
    assert_eq!(shell.number, 4);
    assert_single_element(&structure, &shell, "C");
    assert_distances(&shell, 1.544_556, 1e-5);
    assert_equal_weights(&shell, 1e-9);
}

#[test]
fn cesium_chloride_has_eight_coordination() {
    let structure = cesium_chloride();
    let shell = crystal_nn(&structure, 0, 4.0);
    assert_eq!(shell.number, 8);
    assert_single_element(&structure, &shell, "Cs");
    assert_distances(&shell, 3.559_364, 1e-5);
    assert_equal_weights(&shell, 1e-9);
}

#[test]
fn simple_cubic_polonium_has_six_coordination() {
    let structure = polonium();
    let shell = crystal_nn(&structure, 0, 3.5);
    assert_eq!(shell.number, 6);
    assert_single_element(&structure, &shell, "Po");
    assert_distances(&shell, 3.35, 1e-6);
    assert_equal_weights(&shell, 1e-9);
}

#[test]
fn bcc_distance_weight_keeps_only_the_first_shell() {
    let structure = bcc();
    let shell = crystal_nn(&structure, 0, 4.1);
    assert_eq!(shell.number, 8);
    assert_single_element(&structure, &shell, "Fe");
    assert_distances(&shell, 3.464_102, 1e-5);
    assert_equal_weights(&shell, 1e-9);
}

#[test]
fn every_rock_salt_site_is_six_coordinate() {
    let structure = rocksalt();
    let elements = p1_elements(&structure);
    for (index, element) in elements.iter().enumerate() {
        let shell = crystal_nn(&structure, index, 4.5);
        assert_eq!(shell.number, 6, "atom {index} ({element})");
        let expected = if element == "Na" { "Cl" } else { "Na" };
        assert_single_element(&structure, &shell, expected);
        assert_distances(&shell, 2.82, 1e-6);
    }
}

#[test]
fn every_diamond_site_is_four_coordinate() {
    let structure = diamond();
    let elements = p1_elements(&structure);
    for (index, element) in elements.iter().enumerate() {
        let shell = crystal_nn(&structure, index, 1.7);
        assert_eq!(shell.number, 4, "atom {index} ({element})");
        assert_single_element(&structure, &shell, "C");
        assert_distances(&shell, 1.544_556, 1e-5);
    }
}

#[test]
fn every_cesium_chloride_site_is_eight_coordinate() {
    let structure = cesium_chloride();
    let elements = p1_elements(&structure);
    for (index, element) in elements.iter().enumerate() {
        let shell = crystal_nn(&structure, index, 4.0);
        assert_eq!(shell.number, 8, "atom {index} ({element})");
        let expected = if element == "Cl" { "Cs" } else { "Cl" };
        assert_single_element(&structure, &shell, expected);
        assert_distances(&shell, 3.559_364, 1e-5);
    }
}

#[test]
fn every_polonium_site_is_six_coordinate() {
    let structure = polonium();
    for index in 0..p1_elements(&structure).len() {
        let shell = crystal_nn(&structure, index, 3.5);
        assert_eq!(shell.number, 6, "atom {index}");
        assert_single_element(&structure, &shell, "Po");
        assert_distances(&shell, 3.35, 1e-6);
    }
}

#[test]
fn every_bcc_site_is_eight_coordinate() {
    let structure = bcc();
    for index in 0..p1_elements(&structure).len() {
        let shell = crystal_nn(&structure, index, 4.1);
        assert_eq!(shell.number, 8, "atom {index}");
        assert_single_element(&structure, &shell, "Fe");
        assert_distances(&shell, 3.464_102, 1e-5);
    }
}

#[test]
fn every_l12_site_sees_its_expected_species() {
    let structure = sodium_chloride_l12();
    let elements = p1_elements(&structure);
    for (index, element) in elements.iter().enumerate() {
        let shell = crystal_nn(&structure, index, 3.0);
        if element == "Na" {
            assert_eq!(shell.number, 4, "Na atom {index}");
            assert_single_element(&structure, &shell, "Cl");
            assert_distances(&shell, 2.828_427, 1e-5);
        } else {
            assert_eq!(shell.number, 12, "Cl atom {index}");
            assert_single_element(&structure, &shell, "Na");
            assert_distances(&shell, 2.828_427, 1e-5);
        }
    }
}

#[test]
fn naive_distance_cutoff_overcounts_rock_salt() {
    let structure = rocksalt();
    let cutoff = coordination(&structure, 0, CoordinationStrategy::Cutoff(4.5))
        .expect("cutoff coordination");
    assert_eq!(cutoff.number, 18);

    let shell = crystal_nn(&structure, 0, 4.5);
    assert_eq!(shell.number, 6);
    assert_distances(&shell, 2.82, 1e-6);
}

#[test]
fn crystal_nn_prefers_cation_anion_pairs() {
    let structure = sodium_chloride_l12();
    let cutoff = coordination(&structure, 0, CoordinationStrategy::Cutoff(3.0))
        .expect("cutoff coordination");
    assert_eq!(cutoff.number, 12);
    assert!(elements_of(&structure, &cutoff).contains(&"Cl".to_owned()));
    assert!(elements_of(&structure, &cutoff).contains(&"Na".to_owned()));

    let shell = crystal_nn(&structure, 0, 3.0);
    assert_eq!(shell.number, 4);
    assert_single_element(&structure, &shell, "Cl");
    assert_distances(&shell, 2.828_427, 1e-5);
    assert_equal_weights(&shell, 1e-9);
    assert!(
        shell.effective_number < 0.5,
        "the same-class Na partners should be suppressed, effective {}",
        shell.effective_number
    );
}

#[test]
fn crystal_nn_tolerates_missing_electronegativity() {
    let structure = build(
        "XC",
        Lattice::cubic(4.0).unwrap(),
        &[("X1", "X", [0.0, 0.0, 0.0]), ("C1", "C", [0.5, 0.5, 0.5])],
    );
    let shell = crystal_nn(&structure, 0, 4.1);
    assert_eq!(shell.number, 8);
    assert_single_element(&structure, &shell, "C");
}

#[test]
fn crystal_nn_rejects_invalid_input() {
    let structure = rocksalt();
    assert!(
        coordination(
            &structure,
            0,
            CoordinationStrategy::CrystalNn { cutoff: 0.0 }
        )
        .is_err()
    );
    assert!(
        coordination(
            &structure,
            99,
            CoordinationStrategy::CrystalNn { cutoff: 4.5 }
        )
        .is_err()
    );
}
