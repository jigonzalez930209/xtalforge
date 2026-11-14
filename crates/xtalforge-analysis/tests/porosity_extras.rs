//! Integration tests for the Zeo++-style porosity extras (Task 17.1.2).

#![allow(
    clippy::unwrap_used,
    clippy::expect_used,
    clippy::panic,
    clippy::indexing_slicing
)]

use std::collections::BTreeMap;

use xtalforge_analysis::grid::{
    PorosityOptions, PorosityReport, default_radii, distance_grid, porosity_report,
};
use xtalforge_core::elements::element_by_symbol;
use xtalforge_core::{Lattice, Site, Structure, StructureBuilder};

const GRAMS_PER_AMU: f64 = 1.660_539_066_60e-24;

fn radii(element: &str, radius: f64) -> BTreeMap<String, f64> {
    let mut map = BTreeMap::new();
    map.insert(element.to_owned(), radius);
    map
}

fn one_atom() -> Structure {
    StructureBuilder::new("Na", Lattice::cubic(5.0).unwrap())
        .site(Site::new("Na1", "Na", [0.5, 0.5, 0.5]).unwrap())
        .build()
        .unwrap()
}

/// A simple cubic `a`-channel: a dense carbon wall in the `a`-`c` plane at
/// `b = 0` blocks the `b` axis while the two open slabs percolate along `a`
/// and `c`.
fn walled_channel() -> Structure {
    let mut builder = StructureBuilder::new("wall", Lattice::cubic(8.0).unwrap());
    let mut index = 0usize;
    for x in 0..8 {
        for z in 0..8 {
            let frac = [x as f64 / 8.0, 0.0, z as f64 / 8.0];
            builder = builder.site(Site::new(format!("C{index}"), "C", frac).unwrap());
            index += 1;
        }
    }
    builder.build().unwrap()
}

fn dense() -> Structure {
    StructureBuilder::new("dense", Lattice::cubic(4.0).unwrap())
        .site(Site::new("C1", "C", [0.5, 0.5, 0.5]).unwrap())
        .build()
        .unwrap()
}

fn report_at(probe_radius: f64) -> PorosityReport {
    porosity_report(&one_atom(), &radii("Na", 1.0), 0.5, probe_radius).unwrap()
}

#[test]
fn simple_cubic_channel_percolates_along_a_and_c_only() {
    let structure = walled_channel();
    let table = radii("C", 0.9);
    let grid = distance_grid(&structure, &table, 1.0).unwrap();
    let percolation = grid.percolation(0.5).unwrap();
    assert_eq!(percolation.directions, [true, false, true]);
    assert_eq!(percolation.dimension(), 2);
    assert!(percolation.percolates());
    assert!(!percolation.component_sizes.is_empty());
    assert!(
        percolation
            .component_sizes
            .windows(2)
            .all(|pair| pair[0] >= pair[1])
    );

    let accessible = grid.grid.data.iter().filter(|value| **value >= 0.5).count();
    let total: usize = percolation.component_sizes.iter().sum();
    assert_eq!(total, accessible, "components partition the free samples");

    let blocked = grid.percolation_radius(Some(1)).unwrap();
    assert_eq!(blocked, 0.0, "the b wall blocks every non-negative probe");
    let open = grid.percolation_radius(Some(0)).unwrap();
    assert!(open > 0.5, "the a channel opens for a probe, got {open}");
    let limiting = grid.pore_limiting_diameter().unwrap();
    assert!((limiting - 2.0 * open).abs() < 1e-12);
    assert!(grid.percolation_radius(Some(3)).is_err());

    let report = porosity_report(&structure, &table, 1.0, 0.5).unwrap();
    assert_eq!(report.percolating_directions, [true, false, true]);
    assert!(report.pore_limiting_diameter > 0.0);
}

#[test]
fn pov_and_asa_are_positive_and_monotone_in_the_probe_radius() {
    let probes = [0.25, 0.5, 0.75, 1.0];
    let reports: Vec<PorosityReport> = probes.iter().map(|probe| report_at(*probe)).collect();
    for report in &reports {
        assert!(report.accessible_volume > 0.0);
        assert!(report.probe_occupiable_volume > 0.0);
        assert!(report.pov_per_mass > 0.0);
        assert!(report.surface_area > 0.0);
        assert!(report.asa_per_mass > 0.0);
        assert!(report.pore_limiting_diameter > 0.0);
        assert_eq!(
            report.accessible_volume, report.probe_occupiable_volume,
            "POV is the probe-accessible volume in this voxel model"
        );
        assert!(
            report
                .percolating_directions
                .iter()
                .all(|direction| *direction)
        );
    }
    for pair in reports.windows(2) {
        assert!(pair[1].accessible_volume < pair[0].accessible_volume);
        assert!(pair[1].pov_per_mass < pair[0].pov_per_mass);
        assert!(pair[1].surface_area >= pair[0].surface_area);
        assert!(pair[1].asa_per_mass >= pair[0].asa_per_mass);
    }
}

#[test]
fn per_mass_values_follow_the_cell_mass() {
    let report = report_at(0.5);
    let mass = element_by_symbol("Na").unwrap().mass;
    let grams = mass * GRAMS_PER_AMU;
    let expected_pov = report.accessible_volume * 1.0e-24 / grams;
    let expected_asa = report.surface_area * 1.0e-20 / grams;
    assert!((report.pov_per_mass - expected_pov).abs() < 1e-12 * expected_pov);
    assert!((report.asa_per_mass - expected_asa).abs() < 1e-12 * expected_asa);

    // Occupancy enters the mass, so half an atom doubles both per-mass values
    // while the geometric accessible volume is unchanged.
    let half = StructureBuilder::new("half", Lattice::cubic(5.0).unwrap())
        .site(
            Site::new("Na1", "Na", [0.5, 0.5, 0.5])
                .unwrap()
                .with_occupancy(0.5)
                .unwrap(),
        )
        .build()
        .unwrap();
    let diluted = porosity_report(&half, &radii("Na", 1.0), 0.5, 0.5).unwrap();
    assert_eq!(diluted.accessible_volume, report.accessible_volume);
    assert!((diluted.pov_per_mass - 2.0 * report.pov_per_mass).abs() < 1e-12);
    assert!((diluted.asa_per_mass - 2.0 * report.asa_per_mass).abs() < 1e-12);
}

#[test]
fn dense_structure_has_no_accessible_volume() {
    let report = porosity_report(&dense(), &radii("C", 5.0), 0.5, 0.0).unwrap();
    assert_eq!(report.accessible_volume, 0.0);
    assert_eq!(report.probe_occupiable_volume, 0.0);
    assert_eq!(report.pov_per_mass, 0.0);
    assert_eq!(report.surface_area, 0.0);
    assert_eq!(report.asa_per_mass, 0.0);
    assert_eq!(report.percolating_directions, [false, false, false]);
    assert_eq!(report.pore_limiting_diameter, 0.0);
}

#[test]
fn invalid_spacing_and_probe_radius_are_rejected() {
    let structure = one_atom();
    let table = radii("Na", 1.0);
    for spacing in [0.0, -0.5, f64::NAN, f64::INFINITY] {
        let error = porosity_report(&structure, &table, spacing, 0.5).unwrap_err();
        assert_eq!(error.code(), "invalid_input");
    }
    for probe in [-0.1, f64::NAN, f64::INFINITY] {
        let error = porosity_report(&structure, &table, 0.5, probe).unwrap_err();
        assert_eq!(error.code(), "invalid_input");
    }
    let error = porosity_report(&structure, &BTreeMap::new(), 0.5, 0.5).unwrap_err();
    assert_eq!(error.code(), "invalid_input");

    let empty = StructureBuilder::new("empty", Lattice::cubic(5.0).unwrap())
        .build()
        .unwrap();
    let error = porosity_report(&empty, &BTreeMap::new(), 0.5, 0.5).unwrap_err();
    assert_eq!(error.code(), "invalid_input");

    let grid = distance_grid(&structure, &table, 0.5).unwrap();
    assert!(grid.percolation(f64::NAN).is_err());
    assert!(grid.accessible_surface_area(f64::NAN).is_err());
    assert_eq!(
        grid.percolation(-1.0).unwrap(),
        grid.percolation(0.0).unwrap()
    );
    assert_eq!(
        grid.accessible_surface_area(-1.0).unwrap(),
        grid.accessible_surface_area(0.0).unwrap()
    );
}

#[test]
fn report_is_deterministic_and_serializes() {
    let structure = one_atom();
    let table = default_radii(&structure).unwrap();
    let first = porosity_report(&structure, &table, 0.5, 0.5).unwrap();
    let second = porosity_report(&structure, &table, 0.5, 0.5).unwrap();
    assert_eq!(first, second);

    let text = serde_json::to_string(&first).unwrap();
    let parsed: PorosityReport = serde_json::from_str(&text).unwrap();
    assert_eq!(parsed, first);

    let options = PorosityOptions {
        probe_radius: 0.5,
        spacing: 0.5,
    };
    assert_eq!(options.report(&structure, &table).unwrap(), first);
    assert_eq!(PorosityOptions::default().spacing, 0.5);
    assert_eq!(PorosityOptions::default().probe_radius, 1.2);
}
