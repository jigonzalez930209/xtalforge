//! Substitution tests (Task 20.3.3).

#![allow(
    clippy::unwrap_used,
    clippy::expect_used,
    clippy::panic,
    clippy::indexing_slicing
)]

use xtalforge_generation::{
    PrototypeOptions, SubstitutionOptions, build_prototype, ionic_radii, substituted_structures,
    substitution_candidates, substitution_report,
};

#[test]
fn radii_table_loads() {
    let radii = ionic_radii().expect("radii");
    assert!(radii.len() >= 40);
    assert!(
        radii
            .iter()
            .any(|ion| ion.element == "Na" && ion.charge == 1)
    );
    assert!(
        radii
            .iter()
            .any(|ion| ion.element == "O" && ion.charge == -2)
    );
}

#[test]
fn proposes_alkali_and_halide_substitutions_for_rock_salt() {
    let structure = build_prototype("B1", &PrototypeOptions::default()).expect("B1");
    let options = SubstitutionOptions {
        radius_tolerance: 0.4,
        ..SubstitutionOptions::default()
    };
    let candidates = substitution_candidates(&structure, &options).expect("candidates");
    assert!(
        candidates
            .iter()
            .any(|candidate| candidate.from == "Na" && candidate.to == "K")
    );
    assert!(
        candidates
            .iter()
            .any(|candidate| candidate.from == "Cl" && candidate.to == "F")
    );
    assert!(candidates.iter().all(|candidate| candidate.mismatch <= 0.4));
}

#[test]
fn report_groups_compositions_and_keeps_charge() {
    let structure = build_prototype("B1", &PrototypeOptions::default()).expect("B1");
    let report = substitution_report(
        &structure,
        &SubstitutionOptions {
            radius_tolerance: 0.4,
            ..SubstitutionOptions::default()
        },
    )
    .expect("report");
    assert!(report.neutral);
    assert!(report.compositions.iter().any(|entry| entry == "Na->K"));
}

#[test]
fn applying_a_substitution_replaces_every_site() {
    let structure = build_prototype("B1", &PrototypeOptions::default()).expect("B1");
    let substituted =
        xtalforge_generation::apply_substitution(&structure, "Na", "K").expect("substitute");
    let p1 = substituted.p1().expect("p1");
    assert_eq!(
        p1.atoms.iter().filter(|atom| atom.element == "K").count(),
        4
    );
    assert_eq!(
        p1.atoms.iter().filter(|atom| atom.element == "Na").count(),
        0
    );
    assert_eq!(substituted.operations().len(), structure.operations().len());
}

#[test]
fn structures_are_produced_for_candidates() {
    let structure = build_prototype("B1", &PrototypeOptions::default()).expect("B1");
    let variants = substituted_structures(
        &structure,
        &SubstitutionOptions {
            limit: 4,
            ..SubstitutionOptions::default()
        },
    )
    .expect("variants");
    assert!(!variants.is_empty());
    assert!(
        variants
            .iter()
            .all(|(_, variant)| !variant.sites().is_empty())
    );
}

#[test]
fn invalid_options_are_rejected() {
    let structure = build_prototype("B1", &PrototypeOptions::default()).expect("B1");
    assert!(
        substitution_candidates(
            &structure,
            &SubstitutionOptions {
                radius_tolerance: -1.0,
                ..SubstitutionOptions::default()
            }
        )
        .is_err()
    );
    assert!(
        substitution_candidates(
            &structure,
            &SubstitutionOptions {
                limit: 0,
                ..SubstitutionOptions::default()
            }
        )
        .is_err()
    );
    assert!(xtalforge_generation::apply_substitution(&structure, "Na", "Zz").is_err());
}
