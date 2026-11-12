//! Integration tests for Voronoi, coordination, RDF and fingerprints
//! (Task 9.3).

#![allow(
    clippy::unwrap_used,
    clippy::expect_used,
    clippy::panic,
    clippy::indexing_slicing
)]

use xtalforge_core::{
    CoordinationStrategy, Lattice, RdfOptions, Site, Structure, StructureBuilder,
    bond_angle_distribution, coordination, fingerprint, radial_distribution, shape_measures,
    voronoi,
};

fn bcc() -> Structure {
    StructureBuilder::new("bcc", Lattice::cubic(4.0).unwrap())
        .site(Site::new("Fe1", "Fe", [0.0, 0.0, 0.0]).unwrap())
        .site(Site::new("Fe2", "Fe", [0.5, 0.5, 0.5]).unwrap())
        .build()
        .unwrap()
}

#[test]
fn voronoi_volumes_of_bcc_sum_to_the_cell() {
    let report = voronoi(&bcc(), 4.0).unwrap();
    assert_eq!(report.cells.len(), 2);
    for cell in &report.cells {
        assert!((cell.volume - 32.0).abs() < 0.05, "volume {}", cell.volume);
    }
    assert!((report.total_volume - 64.0).abs() < 0.1);
}

#[test]
fn bcc_coordination_is_eight() {
    let shell = coordination(
        &bcc(),
        0,
        CoordinationStrategy::Voronoi {
            cutoff: 3.8,
            solid_angle_min: 0.1,
        },
    )
    .unwrap();
    assert_eq!(shell.number, 8);
    assert!(!shell.shapes.is_empty());
    let best = &shell.shapes[0];
    assert_eq!(best.shape, "cube");
    assert!(best.cshm < 20.0, "cshm {}", best.cshm);
}

#[test]
fn octahedron_has_zero_shape_measure() {
    let points = [
        [1.0, 0.0, 0.0],
        [-1.0, 0.0, 0.0],
        [0.0, 1.0, 0.0],
        [0.0, -1.0, 0.0],
        [0.0, 0.0, 1.0],
        [0.0, 0.0, -1.0],
    ];
    let measures = shape_measures(&points);
    let octahedron = measures
        .iter()
        .find(|measure| measure.shape == "octahedral")
        .unwrap();
    assert!(octahedron.cshm < 1e-6, "cshm {}", octahedron.cshm);
}

#[test]
fn rdf_has_a_first_peak_at_the_bcc_neighbor_distance() {
    let options = RdfOptions {
        r_max: 3.6,
        bin_width: 0.02,
        smoothing: 3,
    };
    let distribution = radial_distribution(&bcc(), &options).unwrap();
    let peak = distribution
        .total
        .iter()
        .enumerate()
        .max_by(|first, second| first.1.partial_cmp(second.1).unwrap())
        .unwrap()
        .0;
    let radius = distribution.r[peak];
    assert!((radius - 3.464_101).abs() < 0.1, "peak at {radius}");
    assert!(!distribution.partial.is_empty());
}

#[test]
fn fingerprint_compares_like_structures() {
    let options = RdfOptions {
        r_max: 6.0,
        bin_width: 0.1,
        smoothing: 0,
    };
    let first = fingerprint(&bcc(), &options).unwrap();
    let second = fingerprint(&bcc(), &options).unwrap();
    assert!((first.similarity(&second) - 1.0).abs() < 1e-9);
}

#[test]
fn bond_angle_distribution_is_not_empty() {
    let angles = bond_angle_distribution(&bcc(), 4.0, 2.0).unwrap();
    assert!(angles.counts.iter().sum::<f64>() > 0.0);
}
