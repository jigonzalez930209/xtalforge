//! Integration tests for functional groups (Task 1.13).

#![allow(
    clippy::unwrap_used,
    clippy::expect_used,
    clippy::panic,
    clippy::indexing_slicing
)]

use xtalforge_core::{GroupKind, Lattice, Site, StructureBuilder, find_groups, group_counts};

fn build(sites: &[(&str, [f64; 3])]) -> xtalforge_core::Structure {
    let mut builder = StructureBuilder::new("groups", Lattice::cubic(10.0).unwrap());
    for (index, (element, frac)) in sites.iter().enumerate() {
        builder = builder.site(Site::new(format!("{element}{index}"), element, *frac).unwrap());
    }
    builder.build().unwrap()
}

#[test]
fn finds_hydroxyl_in_methanol() {
    let structure = build(&[
        ("C", [0.5, 0.5, 0.5]),
        ("O", [0.56, 0.5, 0.5]),
        ("H", [0.656, 0.5, 0.5]),
    ]);
    let groups = find_groups(&structure).unwrap();
    assert!(groups.iter().any(|group| group.kind == GroupKind::Hydroxyl));
    assert_eq!(
        groups
            .iter()
            .filter(|group| group.kind == GroupKind::Hydroxyl)
            .count(),
        1
    );
}

#[test]
fn finds_carboxylate_and_fluoride() {
    let carboxylate = build(&[
        ("C", [0.5, 0.5, 0.5]),
        ("O", [0.56, 0.5, 0.5]),
        ("O", [0.44, 0.5, 0.5]),
    ]);
    assert_eq!(
        find_groups(&carboxylate)
            .unwrap()
            .iter()
            .filter(|group| group.kind == GroupKind::Carboxylate)
            .count(),
        1
    );

    let fluoride = build(&[("C", [0.5, 0.5, 0.5]), ("F", [0.56, 0.5, 0.5])]);
    assert!(
        find_groups(&fluoride)
            .unwrap()
            .iter()
            .any(|group| group.kind == GroupKind::Fluoride)
    );
}

#[test]
fn counts_groups() {
    let structure = build(&[
        ("C", [0.5, 0.5, 0.5]),
        ("O", [0.56, 0.5, 0.5]),
        ("H", [0.656, 0.5, 0.5]),
    ]);
    let counts = group_counts(&structure).unwrap();
    assert!(counts.contains(&(GroupKind::Hydroxyl, 1)));
}
