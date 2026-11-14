//! Layer-group catalogue tests (Phase 20 support).

#![allow(
    clippy::unwrap_used,
    clippy::expect_used,
    clippy::panic,
    clippy::indexing_slicing
)]

use xtalforge_symmetry::{LayerGroup, LayerWyckoffPosition, layer_groups, layer_number};

#[test]
fn catalogue_has_eighty_default_settings() {
    let groups = layer_groups().expect("groups");
    assert_eq!(groups.len(), 80);
    for (index, group) in groups.iter().enumerate() {
        assert_eq!(group.number as usize, index + 1);
        assert_eq!(layer_number(group.hall_number), Some(group.number));
    }
}

#[test]
fn p1_has_one_operation_and_one_wyckoff_position() {
    let group = LayerGroup::from_number(1).expect("group");
    assert_eq!(group.hm_symbol, "p 1");
    assert_eq!(group.operations.len(), 1);
    let positions = LayerWyckoffPosition::list(group.hall_number).expect("positions");
    assert_eq!(positions.len(), 1);
    assert_eq!(positions[0].multiplicity, 1);
    assert_eq!(positions[0].coordinates, "x,y,z");
}

#[test]
fn centered_layer_has_two_operations_per_generator() {
    let group = LayerGroup::from_number(10).expect("group");
    assert_eq!(group.centering, 'c');
    assert_eq!(group.bravais, "oc");
    assert!(group.operations.len() >= 4);
}

#[test]
fn invalid_numbers_are_rejected() {
    assert!(LayerGroup::from_number(0).is_err());
    assert!(LayerGroup::from_number(81).is_err());
    assert!(LayerGroup::from_hall_number(0).is_err());
    assert!(LayerGroup::from_hall_number(117).is_err());
    assert!(LayerWyckoffPosition::list(0).is_err());
    assert!(LayerWyckoffPosition::list(117).is_err());
}
