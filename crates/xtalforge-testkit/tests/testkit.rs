//! Integration tests for the golden-data comparison helpers.

#![allow(
    clippy::unwrap_used,
    clippy::expect_used,
    clippy::panic,
    clippy::indexing_slicing
)]

use xtalforge_testkit::{
    Atom, Bond, Peak, assert_atom_sets_equal, assert_bond_sets_equal, assert_close,
    assert_peaks_equal, compare_atom_sets, compare_close, compare_peaks,
};

#[test]
fn close_values_respect_both_tolerances() {
    assert_close(1.0, 1.0, 1e-12, 0.0);
    assert_close(100.0, 100.1, 0.0, 1e-3);
    assert!(compare_close(1.0, 1.01, 1e-6, 0.0).is_err());
    assert!(compare_close(f64::NAN, 1.0, 1e-6, 0.0).is_err());
}

#[test]
fn atom_sets_ignore_order_and_translation() {
    let expected = vec![
        Atom::new("Si", [0.0, 0.0, 0.0]),
        Atom::new("Si", [0.25, 0.25, 0.25]),
    ];
    let actual = vec![
        Atom::new("Si", [0.75, 0.75, 0.75]),
        Atom::new("Si", [0.0, 0.0, 0.0]),
    ];
    assert_atom_sets_equal(&expected, &actual, 1e-9);
    assert!(compare_atom_sets(&expected, &actual[..1], 1e-9).is_err());
}

#[test]
fn atom_sets_detect_element_changes() {
    let expected = vec![Atom::new("C", [0.0, 0.0, 0.0])];
    let actual = vec![Atom::new("N", [0.0, 0.0, 0.0])];
    assert!(compare_atom_sets(&expected, &actual, 1e-9).is_err());
}

#[test]
fn bond_sets_ignore_order_and_images() {
    let expected = vec![
        Bond::new(0, 1, [0, 0, 0], 1.0),
        Bond::new(1, 2, [-1, 0, 0], 1.5),
    ];
    let actual = vec![
        Bond::new(1, 2, [-1, 0, 0], 1.5),
        Bond::new(0, 1, [0, 0, 0], 1.0),
    ];
    assert_bond_sets_equal(&expected, &actual, 1e-6);
}

#[test]
fn bond_sets_detect_order_changes() {
    let expected = vec![Bond::new(0, 1, [0, 0, 0], 1.0)];
    let actual = vec![Bond::new(0, 1, [0, 0, 0], 1.5)];
    assert!(xtalforge_testkit::compare_bond_sets(&expected, &actual, 1e-6).is_err());
}

#[test]
fn peaks_are_normalized_before_comparison() {
    let expected = vec![
        Peak {
            position: 10.0,
            intensity: 100.0,
        },
        Peak {
            position: 20.0,
            intensity: 50.0,
        },
    ];
    let actual = vec![
        Peak {
            position: 20.0,
            intensity: 0.5,
        },
        Peak {
            position: 10.0,
            intensity: 1.0,
        },
    ];
    assert_peaks_equal(&expected, &actual, 1e-6, 1e-6);
    assert!(compare_peaks(&expected, &actual[..1], 1e-6, 1e-6).is_err());
}
