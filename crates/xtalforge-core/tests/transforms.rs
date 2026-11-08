//! Integration tests for point-set transforms, finite cells and ellipsoids.

#![allow(
    clippy::unwrap_used,
    clippy::expect_used,
    clippy::panic,
    clippy::indexing_slicing
)]

use xtalforge_core::transforms::{
    align_vector, displacement_ellipsoid, fit_plane, invert, mirror, rotate, translate,
};
use xtalforge_core::{Lattice, Plane, Site, StructureBuilder};

#[test]
fn translation_rotation_mirror_and_inversion() {
    let points = [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0]];
    assert_eq!(
        translate(&points[..1], [1.0, 2.0, 3.0]),
        vec![[2.0, 2.0, 3.0]]
    );
    let rotated = rotate(&points[..1], [0.0, 0.0, 1.0], 90.0, [0.0, 0.0, 0.0]).unwrap();
    assert!((rotated[0][0]).abs() < 1e-12);
    assert!((rotated[0][1] - 1.0).abs() < 1e-12);
    let plane = Plane {
        normal: [0.0, 0.0, 1.0],
        point: [0.0, 0.0, 0.0],
    };
    let mirrored = mirror(&[[1.0, 2.0, 3.0]], plane);
    assert_eq!(mirrored, vec![[1.0, 2.0, -3.0]]);
    let inverted = invert(&[[1.0, 2.0, 3.0]], [0.0, 0.0, 0.0]);
    assert_eq!(inverted, vec![[-1.0, -2.0, -3.0]]);
}

#[test]
fn alignment_and_plane_fit() {
    let aligned = align_vector(
        &[[1.0, 0.0, 0.0]],
        [1.0, 0.0, 0.0],
        [0.0, 1.0, 0.0],
        [0.0, 0.0, 0.0],
    )
    .unwrap();
    assert!(aligned[0][0].abs() < 1e-12);
    assert!((aligned[0][1] - 1.0).abs() < 1e-12);

    let plane = fit_plane(&[[0.0, 0.0, 2.0], [1.0, 0.0, 2.0], [0.0, 1.0, 2.0]]).unwrap();
    assert!((plane.normal[2].abs() - 1.0).abs() < 1e-12);
    assert!(plane.normal[0].abs() < 1e-12);
}

#[test]
fn displacement_ellipsoid_has_expected_axes() {
    let ellipsoid = displacement_ellipsoid([1.0, 4.0, 9.0, 0.0, 0.0, 0.0], 0.5).unwrap();
    assert!((ellipsoid.lengths[0] - 1.538).abs() < 1e-9);
    assert!((ellipsoid.lengths[1] - 3.076).abs() < 1e-9);
    assert!((ellipsoid.lengths[2] - 4.614).abs() < 1e-9);
    for axis in ellipsoid.axes {
        assert!((axis[0] * axis[0] + axis[1] * axis[1] + axis[2] * axis[2] - 1.0).abs() < 1e-9);
    }
    assert!(displacement_ellipsoid([1.0, 1.0, 1.0, 2.0, 0.0, 0.0], 0.5).is_err());
}

#[test]
fn finite_cell_duplicates_corner_atoms() {
    let lattice = Lattice::cubic(4.0).unwrap();
    let corner = StructureBuilder::new("corner", lattice)
        .site(Site::new("C1", "C", [0.0, 0.0, 0.0]).unwrap())
        .build()
        .unwrap();
    assert_eq!(corner.finite_cell(1e-6).unwrap().sites().len(), 8);
    let center = StructureBuilder::new("center", lattice)
        .site(Site::new("C1", "C", [0.5, 0.5, 0.5]).unwrap())
        .build()
        .unwrap();
    assert_eq!(center.finite_cell(1e-6).unwrap().sites().len(), 1);
}
