//! Integration tests for primitive ring perception (Task 1.12).

#![allow(
    clippy::unwrap_used,
    clippy::expect_used,
    clippy::panic,
    clippy::indexing_slicing
)]

use xtalforge_core::{Lattice, Site, StructureBuilder, ring_counts, rings};

fn hexagon() -> xtalforge_core::Structure {
    let mut builder = StructureBuilder::new("hexagon", Lattice::cubic(20.0).unwrap());
    for index in 0..6 {
        let angle = std::f64::consts::TAU * index as f64 / 6.0;
        let x = 8.0 + 1.4 * angle.cos();
        let y = 8.0 + 1.4 * angle.sin();
        builder =
            builder.site(Site::new(format!("C{index}"), "C", [x / 20.0, y / 20.0, 0.5]).unwrap());
    }
    builder.build().unwrap()
}

#[test]
fn finds_one_hexagon() {
    let structure = hexagon();
    let found = rings(&structure, 8).unwrap();
    assert_eq!(found.len(), 1);
    assert_eq!(found[0].size(), 6);
    assert_eq!(found[0].atoms.len(), 6);
}

#[test]
fn counts_rings_by_size() {
    let structure = hexagon();
    assert_eq!(ring_counts(&structure, 8).unwrap(), vec![(6, 1)]);
}

#[test]
fn triangle_is_a_ring() {
    let structure = StructureBuilder::new("triangle", Lattice::cubic(10.0).unwrap())
        .site(Site::new("C0", "C", [0.5, 0.5, 0.5]).unwrap())
        .site(Site::new("C1", "C", [0.56, 0.5, 0.5]).unwrap())
        .site(Site::new("C2", "C", [0.53, 0.554, 0.5]).unwrap())
        .build()
        .unwrap();
    let found = rings(&structure, 4).unwrap();
    assert!(found.iter().any(|ring| ring.size() == 3));
}

#[test]
fn size_limit_is_validated() {
    let structure = hexagon();
    assert!(rings(&structure, 2).is_err());
}
