//! Integration tests for the periodic neighbor search (Task 1.10).

#![allow(
    clippy::unwrap_used,
    clippy::expect_used,
    clippy::panic,
    clippy::indexing_slicing
)]

use xtalforge_core::{Lattice, Site, StructureBuilder, neighbors, neighbors_of};

fn structure(sites: &[(&str, [f64; 3])]) -> xtalforge_core::Structure {
    let mut builder = StructureBuilder::new("test", Lattice::cubic(4.0).unwrap());
    for (index, (element, frac)) in sites.iter().enumerate() {
        builder = builder.site(Site::new(format!("{element}{index}"), element, *frac).unwrap());
    }
    builder.build().unwrap()
}

#[test]
fn nearest_pair_and_image() {
    let structure = structure(&[("C", [0.0, 0.0, 0.0]), ("C", [0.5, 0.0, 0.0])]);
    let found = neighbors(&structure, 2.1).unwrap();
    assert_eq!(found.len(), 2);
    assert!(
        found
            .iter()
            .all(|neighbor| neighbor.i == 0 && neighbor.j == 1)
    );
    assert!(found.iter().any(|neighbor| neighbor.image == [0, 0, 0]));
    assert!(found.iter().any(|neighbor| neighbor.image == [-1, 0, 0]));
    assert!(
        found
            .iter()
            .all(|neighbor| (neighbor.distance - 2.0).abs() < 1e-12)
    );
}

#[test]
fn self_images_are_reported_once() {
    let structure = structure(&[("C", [0.0, 0.0, 0.0])]);
    let found = neighbors(&structure, 4.1).unwrap();
    assert_eq!(found.len(), 6);
    for neighbor in &found {
        assert_eq!(neighbor.i, neighbor.j);
        assert!((neighbor.distance - 4.0).abs() < 1e-12);
    }
}

#[test]
fn periodic_image_between_boundary_atoms() {
    let structure = structure(&[("C", [0.9, 0.0, 0.0]), ("C", [0.1, 0.0, 0.0])]);
    let found = neighbors(&structure, 1.0).unwrap();
    assert_eq!(found.len(), 1);
    assert_eq!(found[0].image, [1, 0, 0]);
    assert!((found[0].distance - 0.8).abs() < 1e-12);
}

#[test]
fn results_are_sorted_and_filterable() {
    let structure = structure(&[
        ("C", [0.0, 0.0, 0.0]),
        ("C", [0.25, 0.0, 0.0]),
        ("C", [0.5, 0.0, 0.0]),
    ]);
    let found = neighbors(&structure, 1.1).unwrap();
    assert!(found.windows(2).all(|pair| {
        (pair[0].i, pair[0].j, pair[0].image) <= (pair[1].i, pair[1].j, pair[1].image)
    }));
    let first = neighbors_of(&structure, 0, 1.1).unwrap();
    assert!(
        first
            .iter()
            .all(|neighbor| neighbor.i == 0 || neighbor.j == 0)
    );
}

#[test]
fn cutoff_is_validated() {
    let structure = structure(&[("C", [0.0, 0.0, 0.0])]);
    assert!(neighbors(&structure, 0.0).is_err());
    assert!(neighbors(&structure, f64::NAN).is_err());
}
