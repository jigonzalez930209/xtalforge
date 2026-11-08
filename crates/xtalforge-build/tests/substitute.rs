//! Integration tests for clearance and substitution (Task 6.3).

#![allow(
    clippy::unwrap_used,
    clippy::expect_used,
    clippy::panic,
    clippy::indexing_slicing
)]

use xtalforge_build::{clearance, substitute};
use xtalforge_core::{Lattice, Site, StructureBuilder};

#[test]
fn clearance_and_substitution() {
    let structure = StructureBuilder::new("test", Lattice::cubic(10.0).unwrap())
        .site(Site::new("C1", "C", [0.0, 0.0, 0.0]).unwrap())
        .site(Site::new("C2", "C", [0.5, 0.5, 0.5]).unwrap())
        .build()
        .unwrap();
    assert!(clearance(&structure, [0.5, 0.0, 0.0]).unwrap() > 4.0);
    let substituted = substitute(&structure, 0, "N").unwrap();
    assert_eq!(substituted.sites()[0].element(), "N");
    assert_eq!(substituted.sites()[1].element(), "C");
    assert!(substitute(&structure, 9, "N").is_err());
    assert!(substitute(&structure, 0, "Zz").is_err());
    assert!(clearance(&structure, [f64::NAN, 0.0, 0.0]).is_err());
}
