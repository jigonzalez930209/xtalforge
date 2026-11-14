//! Prototype library tests (Task 20.3).

#![allow(
    clippy::unwrap_used,
    clippy::expect_used,
    clippy::panic,
    clippy::indexing_slicing
)]

use xtalforge_generation::{
    PrototypeOptions, build_prototype, identify_prototype, prototype_library, prototype_names,
};
use xtalforge_symmetry::{SymmetryOptions, detect_symmetry};

#[test]
fn library_loads_and_has_the_expected_labels() {
    let library = prototype_library().expect("library");
    assert!(library.len() >= 18);
    let names = prototype_names().expect("names");
    for label in [
        "A1", "A2", "A3", "A4", "B1", "B2", "B3", "B4", "C1", "E1", "L12",
    ] {
        assert!(names.iter().any(|name| name == label), "missing {label}");
    }
}

#[test]
fn every_prototype_builds_and_redetects_its_group() {
    for prototype in prototype_library().expect("library") {
        let structure =
            build_prototype(&prototype.label, &PrototypeOptions::default()).expect("build");
        let dataset = detect_symmetry(&structure, SymmetryOptions::default()).expect("detect");
        assert_eq!(
            dataset.number, prototype.space_group,
            "prototype {} detected as {}",
            prototype.label, dataset.number
        );
    }
}

#[test]
fn rock_salt_has_the_expected_composition() {
    let structure = build_prototype("B1", &PrototypeOptions::default()).expect("B1");
    let p1 = structure.p1().expect("p1");
    assert_eq!(p1.atoms.len(), 8);
    assert_eq!(
        p1.atoms.iter().filter(|atom| atom.element == "Na").count(),
        4
    );
    assert_eq!(
        p1.atoms.iter().filter(|atom| atom.element == "Cl").count(),
        4
    );
}

#[test]
fn scale_changes_only_the_cell() {
    let normal = build_prototype("A1", &PrototypeOptions::default()).expect("A1");
    let scaled = build_prototype(
        "A1",
        &PrototypeOptions {
            scale: 2.0,
            ..PrototypeOptions::default()
        },
    )
    .expect("scaled");
    assert!(
        (scaled.lattice().parameters()[0] - 2.0 * normal.lattice().parameters()[0]).abs() < 1e-9
    );
    assert_eq!(normal.sites()[0].frac(), scaled.sites()[0].frac());
}

#[test]
fn identifies_the_built_prototypes() {
    for label in ["A1", "B1", "B2", "E1", "L12"] {
        let structure = build_prototype(label, &PrototypeOptions::default()).expect("build");
        let identification =
            identify_prototype(&structure, &PrototypeOptions::default()).expect("identify");
        assert_eq!(
            identification.label.as_deref(),
            Some(label),
            "identification: {identification:?}"
        );
    }
}

#[test]
fn identifies_a_rock_salt_with_other_elements() {
    let structure = build_prototype("B1", &PrototypeOptions::default()).expect("B1");
    let substituted = xtalforge_generation::apply_substitution(&structure, "Na", "K").expect("K");
    let substituted =
        xtalforge_generation::apply_substitution(&substituted, "Cl", "Br").expect("Br");
    let identification =
        identify_prototype(&substituted, &PrototypeOptions::default()).expect("identify");
    assert_eq!(identification.label.as_deref(), Some("B1"));
}

#[test]
fn unknown_labels_and_scales_are_rejected() {
    assert!(build_prototype("ZZ", &PrototypeOptions::default()).is_err());
    assert!(
        build_prototype(
            "A1",
            &PrototypeOptions {
                scale: 0.0,
                ..PrototypeOptions::default()
            }
        )
        .is_err()
    );
}
