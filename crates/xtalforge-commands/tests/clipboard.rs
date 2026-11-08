//! Integration tests for clipboard commands (Task 3.10.6).

#![allow(
    clippy::unwrap_used,
    clippy::expect_used,
    clippy::panic,
    clippy::indexing_slicing
)]

use std::collections::BTreeSet;

use xtalforge_commands::clipboard::{Fragment, PasteFragment, copy_selection};
use xtalforge_commands::command::{Command, CommandStack, Origin};
use xtalforge_commands::document::Document;
use xtalforge_core::{Lattice, Site, Structure, StructureBuilder};

fn assert_frac(site: &xtalforge_core::Site, expected: [f64; 3]) {
    let actual = site.frac();
    for (actual, expected) in actual.iter().zip(expected.iter()) {
        assert!((actual - expected).abs() < 1e-12);
    }
}

fn base() -> Structure {
    StructureBuilder::new("clipboard", Lattice::cubic(4.0).unwrap())
        .site(Site::new("C1", "C", [0.0, 0.0, 0.0]).unwrap())
        .site(Site::new("O1", "O", [0.5, 0.0, 0.0]).unwrap())
        .build()
        .unwrap()
}

#[test]
fn copy_selection_returns_recentred_cartesian_points() {
    let mut document = Document::new(base());
    assert!(copy_selection(&document).is_empty());
    document.set_selection(vec![0, 1]).unwrap();
    let fragment = copy_selection(&document);
    assert_eq!(fragment.len(), 2);
    assert_eq!(fragment.points[0].0, "C");
    assert_eq!(fragment.points[1].0, "O");
    assert_eq!(fragment.points[0].1, [-1.0, 0.0, 0.0]);
    assert_eq!(fragment.points[1].1, [1.0, 0.0, 0.0]);
    assert_eq!(fragment.centroid(), [0.0, 0.0, 0.0]);
}

#[test]
fn paste_fragment_adds_sites_and_round_trips() {
    let mut document = Document::new(base());
    document.set_selection(vec![0, 1]).unwrap();
    let fragment = copy_selection(&document);

    let mut structure = base();
    let mut stack = CommandStack::new();
    let command = PasteFragment::new(fragment, [1.0, 1.0, 0.0]);
    assert_eq!(command.fragment().len(), 2);
    assert_eq!(command.offset(), [1.0, 1.0, 0.0]);
    stack
        .push(&mut structure, Box::new(command), Origin::User)
        .unwrap();
    assert_eq!(structure.sites().len(), 4);
    assert_frac(&structure.sites()[2], [0.0, 0.25, 0.0]);
    assert_frac(&structure.sites()[3], [0.5, 0.25, 0.0]);
    let labels: BTreeSet<&str> = structure.sites().iter().map(|site| site.label()).collect();
    assert_eq!(labels.len(), 4);

    assert!(stack.undo(&mut structure));
    assert_eq!(structure.sites().len(), 2);
    assert!(stack.redo(&mut structure));
    assert_eq!(structure.sites().len(), 4);
}

#[test]
fn paste_rejects_unknown_elements_and_bad_offsets() {
    let structure = base();
    let unknown = Fragment::new(vec![("Qq".to_owned(), [0.0, 0.0, 0.0])]);
    let error = PasteFragment::new(unknown, [0.0, 0.0, 0.0])
        .apply(&structure)
        .unwrap_err();
    assert_eq!(error.code(), "invalid_input");

    let fragment = Fragment::new(vec![("C".to_owned(), [0.0, 0.0, 0.0])]);
    assert!(
        PasteFragment::new(fragment, [f64::INFINITY, 0.0, 0.0])
            .apply(&structure)
            .is_err()
    );
}

#[test]
fn fragment_serializes_and_keeps_its_centroid() {
    let fragment = Fragment::new(vec![
        ("C".to_owned(), [0.0, 0.0, 0.0]),
        ("C".to_owned(), [2.0, 0.0, 0.0]),
    ]);
    assert_eq!(fragment.centroid(), [1.0, 0.0, 0.0]);
    let json = serde_json::to_string(&fragment).unwrap();
    let decoded: Fragment = serde_json::from_str(&json).unwrap();
    assert_eq!(decoded, fragment);
    assert_eq!(Fragment::default().centroid(), [0.0, 0.0, 0.0]);
}
