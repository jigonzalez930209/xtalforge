//! Integration tests for measures (Task 1.14.1).

#![allow(
    clippy::unwrap_used,
    clippy::expect_used,
    clippy::panic,
    clippy::indexing_slicing
)]

use xtalforge_core::{Lattice, Site, StructureBuilder, measure};

fn points(sites: &[[f64; 3]]) -> xtalforge_core::Structure {
    let mut builder = StructureBuilder::new("measure", Lattice::cubic(10.0).unwrap());
    for (index, frac) in sites.iter().enumerate() {
        builder = builder.site(Site::new(format!("C{index}"), "C", *frac).unwrap());
    }
    builder.build().unwrap()
}

#[test]
fn minimum_image_distance() {
    let structure = points(&[[0.9, 0.0, 0.0], [0.1, 0.0, 0.0]]);
    assert!((measure::distance(&structure, 0, 1).unwrap() - 2.0).abs() < 1e-12);
}

#[test]
fn angle_and_torsion() {
    let structure = points(&[[0.6, 0.5, 0.5], [0.5, 0.5, 0.5], [0.5, 0.6, 0.5]]);
    assert!((measure::angle(&structure, 0, 1, 2).unwrap() - 90.0).abs() < 1e-9);

    let structure = points(&[
        [0.1, 0.1, 0.1],
        [0.2, 0.1, 0.1],
        [0.2, 0.2, 0.1],
        [0.2, 0.2, 0.2],
    ]);
    let value = measure::torsion(&structure, 0, 1, 2, 3).unwrap();
    assert!((value.abs() - 90.0).abs() < 1e-9);
}

#[test]
fn centroid_and_plane_angle() {
    let centroid = measure::centroid(&[[0.0, 0.0, 0.0], [2.0, 0.0, 0.0]]).unwrap();
    assert!((centroid[0] - 1.0).abs() < 1e-12);
    assert!(
        (measure::interplanar_angle([1.0, 0.0, 0.0], [0.0, 1.0, 0.0]).unwrap() - 90.0).abs()
            < 1e-12
    );
}

#[test]
fn indices_and_degenerate_inputs_are_validated() {
    let structure = points(&[[0.0, 0.0, 0.0]]);
    assert!(measure::distance(&structure, 0, 9).is_err());
    assert!(measure::centroid(&[]).is_err());
    assert!(measure::interplanar_angle([0.0, 0.0, 0.0], [1.0, 0.0, 0.0]).is_err());
}
