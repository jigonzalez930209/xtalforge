//! Integration tests for guest placement (Task 6.4).

#![allow(
    clippy::unwrap_used,
    clippy::expect_used,
    clippy::panic,
    clippy::indexing_slicing
)]

use xtalforge_build::place_guest;
use xtalforge_core::{Lattice, Site, StructureBuilder};

#[test]
fn guest_is_placed_when_clear() {
    let host = StructureBuilder::new("host", Lattice::cubic(20.0).unwrap())
        .site(Site::new("C1", "C", [0.0, 0.0, 0.0]).unwrap())
        .build()
        .unwrap();
    let guest = vec![
        ("O".to_owned(), [0.0, 0.0, 0.0]),
        ("H".to_owned(), [1.0, 0.0, 0.0]),
    ];
    let filled = place_guest(&host, &guest, [5.0, 5.0, 5.0], 1.0).unwrap();
    assert_eq!(filled.sites().len(), 3);
    assert!(place_guest(&host, &guest, [0.0, 0.0, 0.0], 1.0).is_err());
    assert!(place_guest(&host, &[], [0.0, 0.0, 0.0], 1.0).is_err());
}
