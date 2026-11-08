//! Integration tests for surface bases and slabs (Task 1.8).

#![allow(
    clippy::unwrap_used,
    clippy::expect_used,
    clippy::panic,
    clippy::indexing_slicing
)]

use xtalforge_core::{Lattice, Site, StructureBuilder, make_slab, plane_angle, surface_basis};

fn dot(first: [i32; 3], second: [i32; 3]) -> i32 {
    first[0] * second[0] + first[1] * second[1] + first[2] * second[2]
}

fn determinant(matrix: [[i32; 3]; 3]) -> i32 {
    let [[a, b, c], [d, e, f], [g, h, i]] = matrix;
    a * (e * i - f * h) - b * (d * i - f * g) + c * (d * h - e * g)
}

#[test]
fn surface_basis_is_unimodular_and_in_the_plane() {
    for hkl in [
        [1, 0, 0],
        [0, 1, 0],
        [0, 0, 1],
        [1, 1, 0],
        [1, 1, 1],
        [2, 1, 0],
    ] {
        let basis = surface_basis(hkl).unwrap();
        assert_eq!(determinant(basis), 1, "hkl={hkl:?}");
        assert_eq!(dot(hkl, basis[0]), 0);
        assert_eq!(dot(hkl, basis[1]), 0);
        assert_eq!(dot(hkl, basis[2]), 1);
    }
    assert!(surface_basis([0, 0, 0]).is_err());
}

#[test]
fn slab_has_the_requested_layers_and_vacuum() {
    let lattice = Lattice::cubic(4.0).unwrap();
    let structure = StructureBuilder::new("C", lattice)
        .site(Site::new("C1", "C", [0.0, 0.0, 0.0]).unwrap())
        .build()
        .unwrap();
    let slab = make_slab(&structure, [1, 0, 0], 2, 10.0, 0.0).unwrap();
    assert_eq!(slab.layers, 2);
    assert_eq!(slab.hkl, [1, 0, 0]);
    assert_eq!(slab.structure.sites().len(), 2);
    assert!(slab.structure.is_p1());
    assert!((slab.thickness - 8.0).abs() < 1e-9);
    assert!((slab.structure.lattice().lengths()[2] - 18.0).abs() < 1e-9);
    assert!((slab.vacuum - 10.0).abs() < 1e-9);
}

#[test]
fn slab_preserves_one_cell_per_layer_for_common_surfaces() {
    let fcc_lattice = Lattice::cubic(4.0).unwrap();
    let fcc = StructureBuilder::new("fcc", fcc_lattice)
        .site(Site::new("Al1", "Al", [0.0, 0.0, 0.0]).unwrap())
        .site(Site::new("Al2", "Al", [0.5, 0.5, 0.0]).unwrap())
        .site(Site::new("Al3", "Al", [0.5, 0.0, 0.5]).unwrap())
        .site(Site::new("Al4", "Al", [0.0, 0.5, 0.5]).unwrap())
        .build()
        .unwrap();
    let slab = make_slab(&fcc, [1, 1, 1], 1, 10.0, 0.0).unwrap();
    assert_eq!(slab.structure.sites().len(), 4);

    let hex_lattice = Lattice::from_parameters(3.2, 3.2, 5.2, 90.0, 90.0, 120.0).unwrap();
    let wurtzite = StructureBuilder::new("wurtzite", hex_lattice)
        .site(Site::new("Zn1", "Zn", [0.0, 0.0, 0.0]).unwrap())
        .site(Site::new("Zn2", "Zn", [1.0 / 3.0, 2.0 / 3.0, 0.5]).unwrap())
        .site(Site::new("O1", "O", [0.0, 0.0, 0.375]).unwrap())
        .site(Site::new("O2", "O", [1.0 / 3.0, 2.0 / 3.0, 0.875]).unwrap())
        .build()
        .unwrap();
    let slab = make_slab(&wurtzite, [0, 0, 1], 1, 10.0, 0.0).unwrap();
    assert_eq!(slab.structure.sites().len(), 4);
}

#[test]
fn slab_rejects_bad_parameters() {
    let lattice = Lattice::cubic(4.0).unwrap();
    let structure = StructureBuilder::new("C", lattice)
        .site(Site::new("C1", "C", [0.0, 0.0, 0.0]).unwrap())
        .build()
        .unwrap();
    assert!(make_slab(&structure, [0, 0, 0], 1, 10.0, 0.0).is_err());
    assert!(make_slab(&structure, [1, 0, 0], 0, 10.0, 0.0).is_err());
    assert!(make_slab(&structure, [1, 0, 0], 1, -1.0, 0.0).is_err());
}

#[test]
fn plane_angles_are_orthogonal_for_axis_planes() {
    assert!((plane_angle([1, 0, 0], [0, 1, 0]) - 90.0).abs() < 1e-12);
    assert!((plane_angle([1, 0, 0], [1, 0, 0])).abs() < 1e-12);
}
