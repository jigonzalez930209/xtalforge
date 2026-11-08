//! Integration tests for selection and tracking (Task 3.11).

#![allow(
    clippy::unwrap_used,
    clippy::expect_used,
    clippy::panic,
    clippy::indexing_slicing
)]

use std::collections::BTreeMap;

use xtalforge_core::{Selection, TrackingTable};

#[test]
fn selection_supports_sets_and_toggles() {
    let mut selection = Selection::new();
    selection.add_atoms([3, 1, 2]);
    assert_eq!(selection.len(), 3);
    assert!(selection.contains(2));
    selection.toggle_atom(2);
    assert!(!selection.contains(2));
    selection.add_bond(7, 4);
    assert!(selection.bonds.contains(&(4, 7)));
    let mut other = Selection::new();
    other.add_atom(1);
    other.add_atom(9);
    selection.union(&other);
    assert!(selection.contains(9));
    selection.intersect(&other);
    assert_eq!(selection.atoms.len(), 2);
    selection.clear();
    assert!(selection.is_empty());
}

#[test]
fn tracking_assigns_and_remaps_ids() {
    let mut table = TrackingTable::new();
    let first = table.assign(0);
    let second = table.assign(1);
    assert_eq!(table.site_of(first), Some(0));
    assert_eq!(table.id_of(1), Some(second));
    table.remap(&BTreeMap::from([(0, 5)]));
    assert_eq!(table.site_of(first), Some(5));
    assert_eq!(table.site_of(second), None);
    table.remove(first);
    assert_eq!(table.site_of(first), None);
}
