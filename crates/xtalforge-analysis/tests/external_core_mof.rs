//! Porosity and PXRD on real CoRE MOF 2019 frameworks (Task 0.9).
//!
//! HKUST-1 and ZIF-8 are the two most studied CoRE MOF entries; they check
//! that the geometric analysis produces physically bounded numbers on a
//! charged Cu paddlewheel and a tetrahedral Zn imidazolate.

#![allow(
    clippy::unwrap_used,
    clippy::expect_used,
    clippy::panic,
    clippy::indexing_slicing
)]

use xtalforge_analysis::{PxrdOptions, calculate_pxrd, default_radii, distance_grid};
use xtalforge_io::read_cif;

const HKUST: &str = include_str!("../../../tests/external/core_mof/FIQCEN_clean_pacman.cif");
const ZIF8: &str = include_str!("../../../tests/external/core_mof/OFERUN_clean_pacman.cif");

fn check(name: &str, text: &str) {
    let structure = read_cif(text).unwrap_or_else(|error| panic!("{name}: {error}"));

    let pattern = calculate_pxrd(&structure, &PxrdOptions::default())
        .unwrap_or_else(|error| panic!("{name} PXRD: {error}"));
    assert!(pattern.reflections.len() > 10, "{name} reflections");
    assert!(pattern.max_intensity() > 0.0, "{name} intensity");
    assert!(
        pattern
            .points
            .iter()
            .all(|point| point.intensity.is_finite() && point.two_theta.is_finite()),
        "{name} profile"
    );

    let radii = default_radii(&structure).unwrap_or_else(|error| panic!("{name}: {error}"));
    let grid = distance_grid(&structure, &radii, 0.8)
        .unwrap_or_else(|error| panic!("{name} grid: {error}"));
    let accessible = grid
        .accessible_volume(1.2)
        .unwrap_or_else(|error| panic!("{name} accessible: {error}"));
    let fraction = accessible / grid.cell_volume();
    assert!(
        (0.0..=1.0).contains(&fraction),
        "{name} accessible fraction {fraction}"
    );
    assert!(fraction > 0.05, "{name} is not porous: {fraction}");
    assert!(grid.surface_area() > 0.0, "{name} surface area");
    assert!(
        grid.largest_included_diameter() > 0.0,
        "{name} largest included sphere"
    );
}

#[test]
fn hkust1_is_porous_and_gives_a_pattern() {
    check("HKUST-1", HKUST);
}

#[test]
fn zif8_is_porous_and_gives_a_pattern() {
    check("ZIF-8", ZIF8);
}
