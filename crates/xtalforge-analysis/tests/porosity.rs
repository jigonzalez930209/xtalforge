//! Integration tests for distance grids and porosity (Tasks 5.2 and 5.3).

#![allow(
    clippy::unwrap_used,
    clippy::expect_used,
    clippy::panic,
    clippy::indexing_slicing
)]

use std::collections::BTreeMap;
use std::f64::consts::PI;

use xtalforge_analysis::grid::{default_radii, distance_grid};
use xtalforge_core::elements::element_by_symbol;
use xtalforge_core::{Lattice, Site, Structure, StructureBuilder};

fn one_atom(frac: [f64; 3]) -> Structure {
    StructureBuilder::new("Na", Lattice::cubic(5.0).unwrap())
        .site(Site::new("Na1", "Na", frac).unwrap())
        .build()
        .unwrap()
}

fn radii(element: &str, radius: f64) -> BTreeMap<String, f64> {
    let mut map = BTreeMap::new();
    map.insert(element.to_owned(), radius);
    map
}

#[test]
fn distances_use_minimum_images_and_are_negative_inside_atoms() {
    let structure = one_atom([0.5, 0.5, 0.5]);
    let grid = distance_grid(&structure, &radii("Na", 1.0), 0.5).unwrap();
    assert_eq!(grid.grid.shape, [10, 10, 10]);
    assert!(grid.grid.periodic);
    assert_eq!(grid.len(), 1000);
    assert!((grid.grid.get(5, 5, 5).unwrap() + 1.0).abs() < 1e-12);
    assert!(grid.grid.get(0, 0, 0).unwrap() > 0.0);

    // An atom just across the cell face is 0.25 Angstrom away, not 4.75.
    let near_face = one_atom([0.05, 0.5, 0.5]);
    let grid = distance_grid(&near_face, &radii("Na", 1.0), 0.5).unwrap();
    let value = grid.grid.get(0, 5, 5).unwrap();
    assert!(
        (value + 0.75).abs() < 1e-12,
        "wrapped surface distance {value}"
    );
}

#[test]
fn volumes_are_monotonic_in_the_probe_radius() {
    let structure = one_atom([0.5, 0.5, 0.5]);
    let grid = distance_grid(&structure, &radii("Na", 1.0), 0.5).unwrap();
    assert!((grid.cell_volume() - 125.0).abs() < 1e-9);
    let free = grid.accessible_volume(0.0).unwrap();
    assert!(free > 0.0 && free < 125.0);
    assert_eq!(grid.pore_volume(), free);
    let mut previous = free;
    for probe in [0.25, 0.5, 0.75, 1.0, 2.0] {
        let volume = grid.accessible_volume(probe).unwrap();
        assert!(volume <= previous + 1e-12);
        previous = volume;
    }
    assert_eq!(grid.accessible_volume(10.0).unwrap(), 0.0);
    assert!(grid.accessible_volume(-1.0).unwrap() == free);
    assert!(grid.accessible_volume(f64::NAN).is_err());
}

#[test]
fn surface_area_is_positive_and_close_to_the_sphere() {
    let structure = one_atom([0.5, 0.5, 0.5]);
    let grid = distance_grid(&structure, &radii("Na", 1.0), 0.5).unwrap();
    let area = grid.surface_area();
    let sphere = 4.0 * PI;
    assert!(area > 0.5 * sphere && area < 2.0 * sphere, "area {area}");
    let surface = grid.contact_surface();
    assert!(!surface.is_empty());
    for point in surface {
        assert!(point.iter().all(|value| value.is_finite()));
    }
}

#[test]
fn empty_cell_has_the_whole_cell_accessible() {
    let structure = StructureBuilder::new("empty", Lattice::cubic(5.0).unwrap())
        .build()
        .unwrap();
    let grid = distance_grid(&structure, &BTreeMap::new(), 0.5).unwrap();
    assert!(grid.grid.data.iter().all(|value| value.is_infinite()));
    assert_eq!(grid.accessible_volume(3.0).unwrap(), 125.0);
    assert_eq!(grid.pore_volume(), 125.0);
    assert_eq!(grid.surface_area(), 0.0);
    assert!(grid.contact_surface().is_empty());
    assert!(grid.pore_diameter(0.5).is_err());
    assert!(grid.largest_included_diameter().is_infinite());
}

#[test]
fn invalid_spacing_and_radii_are_rejected() {
    let structure = one_atom([0.5, 0.5, 0.5]);
    for spacing in [0.0, -0.5, f64::NAN, f64::INFINITY] {
        let error = distance_grid(&structure, &radii("Na", 1.0), spacing).unwrap_err();
        assert_eq!(error.code(), "invalid_input");
    }
    let error = distance_grid(&structure, &BTreeMap::new(), 0.5).unwrap_err();
    assert_eq!(error.code(), "invalid_input");
    for radius in [0.0, -1.0, f64::NAN] {
        let error = distance_grid(&structure, &radii("Na", radius), 0.5).unwrap_err();
        assert_eq!(error.code(), "invalid_input");
    }
    // Too fine for the hard sample limit; rejected before allocation.
    let error = distance_grid(&structure, &radii("Na", 1.0), 0.01).unwrap_err();
    assert_eq!(error.code(), "invalid_input");
}

#[test]
fn default_radii_come_from_the_element_table() {
    let structure = one_atom([0.5, 0.5, 0.5]);
    let table = default_radii(&structure).unwrap();
    let expected = element_by_symbol("Na").unwrap().vdw_radius;
    assert_eq!(table.get("Na").copied(), Some(expected));
    let grid = distance_grid(&structure, &table, 0.5).unwrap();
    assert_eq!(grid.radii, table);
}

#[test]
fn pore_size_histogram_accounts_for_every_free_sample() {
    let structure = one_atom([0.5, 0.5, 0.5]);
    let grid = distance_grid(&structure, &radii("Na", 1.0), 0.5).unwrap();
    let bins = grid.pore_diameter(0.5).unwrap();
    assert!(!bins.is_empty());
    let total: usize = bins.iter().map(|bin| bin.count).sum();
    let expected = grid.grid.data.iter().filter(|value| **value >= 0.0).count();
    assert_eq!(total, expected);
    for pair in bins.windows(2) {
        assert!(pair[0].diameter_max <= pair[1].diameter_min + 1e-12);
    }
    let volume: f64 = bins.iter().map(|bin| bin.volume).sum();
    assert!((volume - grid.pore_volume()).abs() < 1e-9);
    assert!(grid.largest_included_diameter() > 0.0);
    for invalid in [0.0, -0.5, f64::NAN] {
        assert!(grid.pore_diameter(invalid).is_err());
    }
}

#[test]
fn distance_grid_is_deterministic_and_samples_by_length() {
    let structure = StructureBuilder::new("box", Lattice::orthorhombic(10.0, 5.0, 5.0).unwrap())
        .site(Site::new("Na1", "Na", [0.5, 0.5, 0.5]).unwrap())
        .build()
        .unwrap();
    let table = radii("Na", 2.27);
    let first = distance_grid(&structure, &table, 1.0).unwrap();
    let second = distance_grid(&structure, &table, 1.0).unwrap();
    assert_eq!(first, second);
    assert_eq!(first.grid.shape, [10, 5, 5]);
    assert_eq!(first.grid.data.len(), 250);
}
