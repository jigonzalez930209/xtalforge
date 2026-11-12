//! Integration tests for bond valence, surfaces, defects, orderings, SQS and
//! interfaces (Tasks 9.3-9.5).

#![allow(
    clippy::unwrap_used,
    clippy::expect_used,
    clippy::panic,
    clippy::indexing_slicing
)]

use xtalforge_analysis::{
    BvsOptions, InterfaceOptions, OrderingOptions, SqsOptions, bond_valence_sums, interstitials,
    lattice_matches, orderings, oxidation_states, slab_terminations, special_quasirandom,
    substitutions, unique_miller_indices, vacancies,
};
use xtalforge_core::{Lattice, Site, Structure, StructureBuilder};

fn nacl() -> Structure {
    let a = 5.64;
    let positions = [
        ("Na", [0.0, 0.0, 0.0]),
        ("Na", [0.0, 0.5, 0.5]),
        ("Na", [0.5, 0.0, 0.5]),
        ("Na", [0.5, 0.5, 0.0]),
        ("Cl", [0.5, 0.5, 0.5]),
        ("Cl", [0.5, 0.0, 0.0]),
        ("Cl", [0.0, 0.5, 0.0]),
        ("Cl", [0.0, 0.0, 0.5]),
    ];
    let mut builder = StructureBuilder::new("NaCl", Lattice::cubic(a).unwrap());
    for (index, (element, position)) in positions.iter().enumerate() {
        builder = builder.site(Site::new(format!("{element}{index}"), element, *position).unwrap());
    }
    builder.build().unwrap()
}

fn cscl() -> Structure {
    StructureBuilder::new("CsCl", Lattice::cubic(4.2).unwrap())
        .site(Site::new("Cs1", "Cs", [0.0, 0.0, 0.0]).unwrap())
        .site(Site::new("Cl1", "Cl", [0.5, 0.5, 0.5]).unwrap())
        .build()
        .unwrap()
}

#[test]
fn bond_valence_sums_of_nacl_are_about_one() {
    let report = bond_valence_sums(&nacl(), &BvsOptions::default()).unwrap();
    let sodium = report
        .sites
        .iter()
        .find(|site| site.element == "Na")
        .unwrap();
    let chlorine = report
        .sites
        .iter()
        .find(|site| site.element == "Cl")
        .unwrap();
    assert!((sodium.bvs - 1.0).abs() < 0.2, "Na BVS {}", sodium.bvs);
    assert!((chlorine.bvs + 1.0).abs() < 0.2, "Cl BVS {}", chlorine.bvs);
}

#[test]
fn oxidation_states_of_nacl_are_neutral() {
    let report = oxidation_states(&nacl(), &BvsOptions::default()).unwrap();
    assert_eq!(report.by_element.get("Na"), Some(&1));
    assert_eq!(report.by_element.get("Cl"), Some(&-1));
    assert!(report.neutral, "charge {}", report.total_charge);
}

#[test]
fn unique_miller_indices_of_a_cubic_crystal() {
    let indices = unique_miller_indices(&cscl(), 1, 1e-5).unwrap();
    assert_eq!(indices.len(), 3);
    let hundred = indices.iter().find(|index| index.hkl == [1, 0, 0]).unwrap();
    assert_eq!(hundred.multiplicity, 6);
    let eleven = indices.iter().find(|index| index.hkl == [1, 1, 0]).unwrap();
    assert_eq!(eleven.multiplicity, 12);
    let oneoneone = indices.iter().find(|index| index.hkl == [1, 1, 1]).unwrap();
    assert_eq!(oneoneone.multiplicity, 8);
}

#[test]
fn a_two_species_slab_has_two_terminations() {
    let report = slab_terminations(&cscl(), [1, 0, 0], 1, 10.0, false).unwrap();
    assert_eq!(report.n_terminations, 2);
    assert!(report.terminations.iter().all(|entry| entry.n_atoms > 0));
}

#[test]
fn vacancy_and_substitution_are_symmetry_unique() {
    let vacancies = vacancies(&cscl(), 8.0).unwrap();
    assert_eq!(vacancies.len(), 2);
    assert!(vacancies.iter().all(|defect| defect.image_distance >= 8.0));
    let substitutions = substitutions(&cscl(), "Na", 6.0).unwrap();
    // Both the Cs and the Cl orbits accept a Na substitution.
    assert_eq!(substitutions.len(), 2);
    assert!(interstitials(&cscl(), "H", 6.0, 4.0, 0.8).is_ok());
}

fn disordered() -> Structure {
    StructureBuilder::new("ZnO", Lattice::cubic(4.0).unwrap())
        .site(
            Site::new("Zn1", "Zn", [0.0, 0.0, 0.0])
                .unwrap()
                .with_occupancy(0.5)
                .unwrap(),
        )
        .site(
            Site::new("O1", "O", [0.0, 0.0, 0.0])
                .unwrap()
                .with_occupancy(0.5)
                .unwrap(),
        )
        .build()
        .unwrap()
}

#[test]
fn orderings_enumerate_and_rank_decorations() {
    let options = OrderingOptions {
        supercell: [2, 1, 1],
        max_orderings: 32,
        ..OrderingOptions::default()
    };
    let orderings = orderings(&disordered(), &options).unwrap();
    assert!(!orderings.is_empty());
    assert!(orderings.iter().any(|ordering| ordering.energy.is_some()));
    assert!(
        orderings
            .iter()
            .all(|ordering| ordering.assignment.len() == 2)
    );
}

#[test]
fn sqs_reports_pair_and_triplet_correlations() {
    let options = SqsOptions {
        supercell: [2, 2, 2],
        seed: 11,
        cutoff: 5.0,
        iterations: 2000,
        ..SqsOptions::default()
    };
    let result = special_quasirandom(&disordered(), &options).unwrap();
    assert_eq!(result.assignment.len(), 8);
    assert!(!result.shells.is_empty());
    assert!(result.objective.is_finite());
    assert!(result.objective >= 0.0);
}

#[test]
fn identical_slabs_have_a_zero_strain_match() {
    let options = InterfaceOptions::default();
    let first = StructureBuilder::new("substrate", Lattice::orthorhombic(4.0, 4.0, 20.0).unwrap())
        .site(Site::new("Si1", "Si", [0.0, 0.0, 0.0]).unwrap())
        .build()
        .unwrap();
    let second = StructureBuilder::new("film", Lattice::orthorhombic(4.0, 4.0, 8.0).unwrap())
        .site(Site::new("Si1", "Si", [0.0, 0.0, 0.0]).unwrap())
        .build()
        .unwrap();
    let matches = lattice_matches(&second, &first, &options).unwrap();
    assert!(!matches.is_empty());
    assert!(matches[0].strain < 1e-9);
}
