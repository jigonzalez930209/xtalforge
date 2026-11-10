//! Integration tests for graphene and nanotube builders (Tasks 6.5, 6.6).

#![allow(
    clippy::unwrap_used,
    clippy::expect_used,
    clippy::panic,
    clippy::indexing_slicing
)]

use xtalforge_build::{graphene, nanotube};

#[test]
fn graphene_sheet_has_two_atoms_per_cell() {
    let sheet = graphene(3, 2, 1, 3.35, 15.0).unwrap();
    assert_eq!(sheet.sites().len(), 12);
    assert!(sheet.is_p1());
    let [a, _b, c] = sheet.lattice().lengths();
    assert!((a - 3.0 * 2.4595).abs() < 1e-2);
    assert!((c - 18.35).abs() < 1e-9);
    let graph = sheet.bond_graph().unwrap();
    assert!(!graph.bonds.is_empty());
}

#[test]
fn graphene_geometry_is_a_honeycomb_for_every_repeat() {
    for (nx, ny) in [(1, 1), (3, 3), (3, 2), (2, 4)] {
        let sheet = graphene(nx, ny, 1, 3.35, 15.0).unwrap();
        let atoms = sheet.sites().len();
        assert_eq!(atoms, (2 * nx * ny) as usize);
        let graph = sheet.bond_graph().unwrap();
        // Every carbon has three neighbours at the C-C bond length.
        assert_eq!(graph.bonds.len(), 3 * atoms / 2, "{nx}x{ny}");
        for bond in &graph.bonds {
            assert!((bond.distance - 1.42).abs() < 1e-3, "{nx}x{ny}");
        }
        let hexagons = xtalforge_core::rings(&sheet, 8).unwrap();
        assert_eq!(hexagons.len(), (nx * ny) as usize, "{nx}x{ny}");
        assert!(hexagons.iter().all(|ring| ring.size() == 6));
    }
}

#[test]
fn stacked_graphene_has_layers() {
    let sheet = graphene(1, 1, 3, 3.35, 10.0).unwrap();
    assert_eq!(sheet.sites().len(), 6);
    assert!((sheet.lattice().lengths()[2] - 20.05).abs() < 1e-9);
}

#[test]
fn nanotube_atom_count_and_bonds_follow_the_chirality() {
    // N = 2 (n^2 + n m + m^2) / d_R hexagons per translational cell, 2 N atoms.
    for (n, m, length) in [(5, 5, 1), (9, 0, 1), (6, 4, 1), (5, 5, 2)] {
        let tube = nanotube(n, m, length, 10.0).unwrap();
        let d_r = gcd(2 * n + m, 2 * m + n);
        let atoms = 4 * (n * n + n * m + m * m) / d_r * length as i32;
        assert_eq!(tube.sites().len(), atoms as usize, "({n},{m}) x {length}");
        let graph = tube.bond_graph().unwrap();
        // Every carbon has three neighbours; the chord is a little shorter
        // than the 1.42 Angstrom arc.
        assert_eq!(graph.bonds.len(), 3 * atoms as usize / 2, "({n},{m})");
        for bond in &graph.bonds {
            assert!((bond.distance - 1.42).abs() < 0.05, "({n},{m})");
        }
    }
}

fn gcd(first: i32, second: i32) -> i32 {
    if second == 0 {
        first.abs()
    } else {
        gcd(second, first % second)
    }
}

#[test]
fn nanotube_rolls_a_chiral_sheet() {
    let tube = nanotube(5, 5, 2, 10.0).unwrap();
    assert!(tube.sites().len() > 10);
    let [a, _b, c] = tube.lattice().lengths();
    assert!(a > 0.0);
    assert!(c > 10.0);
    assert!(nanotube(0, 0, 1, 10.0).is_err());
    assert!(nanotube(5, 5, 0, 10.0).is_err());
}
