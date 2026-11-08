//! Integration tests for periodic honeycomb nanoribbons (Task 6.8.3).

#![allow(
    clippy::unwrap_used,
    clippy::expect_used,
    clippy::panic,
    clippy::indexing_slicing
)]

use xtalforge_build::nanoribbon;

#[test]
fn ribbon_has_two_atoms_per_cell_and_a_finite_cell() {
    let ribbon = nanoribbon(2, 3, 12.0).unwrap();
    assert_eq!(ribbon.sites().len(), 2 * 2 * 3);
    assert!(ribbon.is_p1());
    let lengths = ribbon.lattice().lengths();
    assert!(
        lengths
            .iter()
            .all(|value| value.is_finite() && *value > 0.0)
    );
    assert!(
        ribbon
            .cart_sites()
            .iter()
            .flatten()
            .all(|value| value.is_finite())
    );
}

#[test]
fn vacuum_increases_the_b_axis() {
    let tight = nanoribbon(2, 3, 0.0).unwrap();
    let wide = nanoribbon(2, 3, 12.0).unwrap();
    let growth = wide.lattice().lengths()[1] - tight.lattice().lengths()[1];
    assert!((growth - 12.0).abs() < 1e-9, "b grew by {growth}");
    assert_eq!(tight.sites(), wide.sites());
}

#[test]
fn zero_size_and_bad_vacuum_error() {
    assert!(nanoribbon(0, 3, 12.0).is_err());
    assert!(nanoribbon(2, 0, 12.0).is_err());
    assert!(nanoribbon(2, 3, -1.0).is_err());
    assert!(nanoribbon(2, 3, f64::NAN).is_err());
}
