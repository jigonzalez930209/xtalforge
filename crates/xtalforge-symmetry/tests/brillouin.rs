//! Brillouin-zone geometry tests (Task 2.10.1).

#![allow(
    clippy::unwrap_used,
    clippy::expect_used,
    clippy::panic,
    clippy::indexing_slicing
)]

use xtalforge_core::Lattice;
use xtalforge_symmetry::brillouin_zone;

#[test]
fn simple_cubic_zone_is_a_cube() {
    let zone = brillouin_zone(&Lattice::cubic(4.0).unwrap()).unwrap();
    assert_eq!(
        (zone.n_vertices(), zone.n_faces(), zone.edges.len()),
        (8, 6, 12)
    );
    for vertex in &zone.vertices {
        for value in vertex {
            assert!((value.abs() - std::f64::consts::PI / 4.0).abs() < 1e-9);
        }
    }
}

#[test]
fn orthorhombic_zone_is_a_rectangular_box() {
    let lattice = Lattice::new([[5.0, 0.0, 0.0], [0.0, 6.0, 0.0], [0.0, 0.0, 7.0]]).unwrap();
    let zone = brillouin_zone(&lattice).unwrap();
    assert_eq!(zone.n_vertices(), 8);
    assert_eq!(zone.n_faces(), 6);
    assert_eq!(zone.edges.len(), 12);
}

#[test]
fn anisotropic_hexagonal_cell_gives_a_hexagonal_prism() {
    let lattice = Lattice::new([
        [2.46, 0.0, 0.0],
        [-1.23, 2.130_422_493_309_721, 0.0],
        [0.0, 0.0, 20.0],
    ])
    .unwrap();
    let zone = brillouin_zone(&lattice).unwrap();
    assert_eq!(
        (zone.n_vertices(), zone.n_faces(), zone.edges.len()),
        (12, 8, 18)
    );
}

#[test]
fn fcc_primitive_lattice_gives_a_truncated_octahedron() {
    let a = 4.0;
    let lattice = Lattice::new([
        [0.0, a / 2.0, a / 2.0],
        [a / 2.0, 0.0, a / 2.0],
        [a / 2.0, a / 2.0, 0.0],
    ])
    .unwrap();
    let zone = brillouin_zone(&lattice).unwrap();
    assert_eq!(zone.n_vertices(), 24);
    assert_eq!(zone.n_faces(), 14);
    assert_eq!(zone.edges.len(), 36);
}

#[test]
fn bcc_primitive_lattice_gives_a_rhombic_dodecahedron() {
    let a = 4.0;
    let lattice = Lattice::new([
        [-a / 2.0, a / 2.0, a / 2.0],
        [a / 2.0, -a / 2.0, a / 2.0],
        [a / 2.0, a / 2.0, -a / 2.0],
    ])
    .unwrap();
    let zone = brillouin_zone(&lattice).unwrap();
    assert_eq!(zone.n_vertices(), 14);
    assert_eq!(zone.n_faces(), 12);
    assert_eq!(zone.edges.len(), 24);
}

#[test]
fn negatively_oriented_orthorhombic_cell_still_has_six_faces() {
    let lattice = Lattice::new([
        [-20.0511, 0.0, 0.0],
        [0.0, -19.8757, 0.0],
        [0.0, 0.0, 13.36823],
    ])
    .unwrap();
    let zone = brillouin_zone(&lattice).unwrap();
    println!("faces {:?}", zone.faces);
    println!("vertices {:?}", zone.vertices);
    assert_eq!(
        (zone.n_vertices(), zone.n_faces(), zone.edges.len()),
        (8, 6, 12)
    );
}
