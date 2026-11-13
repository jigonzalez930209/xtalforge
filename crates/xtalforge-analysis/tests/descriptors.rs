//! Integration tests for the ML descriptors (Task 16.3.3).

#![allow(
    clippy::unwrap_used,
    clippy::expect_used,
    clippy::panic,
    clippy::indexing_slicing
)]

use xtalforge_analysis::{
    CoulombMatrixOptions, MAX_SOAP_L, SoapLikeOptions, composition_features, coulomb_matrix,
    rdf_descriptor, soap_like_descriptor,
};
use xtalforge_core::{Lattice, RdfOptions, Site, Structure, StructureBuilder};

fn nacl() -> Structure {
    let positions = [
        ("Na", [0.0, 0.0, 0.0]),
        ("Na", [0.0, 0.5, 0.5]),
        ("Na", [0.5, 0.0, 0.5]),
        ("Na", [0.5, 0.5, 0.0]),
        ("Cl", [0.5, 0.5, 0.5]),
        ("Cl", [0.5, 0.0, 0.0]),
        ("Cl", [0.0, 0.5, 0.0]),
        ("Cl", [0.0, 0.0, 0.5]),
    ];
    let mut builder = StructureBuilder::new("NaCl", Lattice::cubic(5.64).unwrap());
    for (index, (element, position)) in positions.iter().enumerate() {
        builder = builder.site(Site::new(format!("{element}{index}"), element, *position).unwrap());
    }
    builder.build().unwrap()
}

fn cscl() -> Structure {
    let builder = StructureBuilder::new("CsCl", Lattice::cubic(4.2).unwrap())
        .site(Site::new("Cs1", "Cs", [0.0, 0.0, 0.0]).unwrap())
        .site(Site::new("Cl1", "Cl", [0.5, 0.5, 0.5]).unwrap());
    builder.build().unwrap()
}

fn reversed_cscl() -> Structure {
    let builder = StructureBuilder::new("CsCl", Lattice::cubic(4.2).unwrap())
        .site(Site::new("Cl1", "Cl", [0.5, 0.5, 0.5]).unwrap())
        .site(Site::new("Cs1", "Cs", [0.0, 0.0, 0.0]).unwrap());
    builder.build().unwrap()
}

fn empty() -> Structure {
    StructureBuilder::new("empty", Lattice::cubic(5.0).unwrap())
        .build()
        .unwrap()
}

fn polonium() -> Structure {
    StructureBuilder::new("Po", Lattice::cubic(3.35).unwrap())
        .site(Site::new("Po1", "Po", [0.0, 0.0, 0.0]).unwrap())
        .build()
        .unwrap()
}

fn assert_close(left: &[f64], right: &[f64], tolerance: f64) {
    assert_eq!(left.len(), right.len());
    for (index, (value, expected)) in left.iter().zip(right.iter()).enumerate() {
        let scale = expected.abs().max(1.0);
        assert!(
            (value - expected).abs() <= tolerance * scale,
            "entry {index}: {value} vs {expected}"
        );
    }
}

#[test]
fn composition_features_report_atom_fractions() {
    let features = composition_features(&nacl()).unwrap();
    assert_eq!(features.elements, ["Cl", "Na"]);
    assert_eq!(features.counts, [4, 4]);
    assert_eq!(features.fractions, [0.5, 0.5]);
    assert_eq!(features.n_atoms, 8);
    assert_eq!(features.fraction("Na"), 0.5);
    assert_eq!(features.fraction("O"), 0.0);
    assert_eq!(features.count("Cl"), 4);
    assert_eq!(features.count("O"), 0);
    assert_eq!(features.to_vector(), [0.5, 0.5]);
    assert!((features.mean_atomic_number - 14.0).abs() < 1e-12);
    let expected_mass = (4.0
        * xtalforge_core::elements::element_by_symbol("Na")
            .unwrap()
            .mass
        + 4.0
            * xtalforge_core::elements::element_by_symbol("Cl")
                .unwrap()
                .mass)
        / 8.0;
    assert!((features.mean_mass - expected_mass).abs() < 1e-9);
    assert!((features.metal_fraction - 0.5).abs() < 1e-12);
    let electronegativity = features.mean_electronegativity.unwrap();
    assert!(
        (electronegativity - 2.045).abs() < 1e-6,
        "{electronegativity}"
    );

    let json = features.to_json().unwrap();
    assert!(json.contains("\"elements\""));
    let parsed: xtalforge_analysis::CompositionFeatures = serde_json::from_str(&json).unwrap();
    assert_eq!(parsed, features);
}

#[test]
fn composition_features_reject_empty_structures() {
    assert!(composition_features(&empty()).is_err());
}

#[test]
fn the_rdf_descriptor_matches_the_core_curve() {
    let options = RdfOptions {
        r_max: 8.0,
        bin_width: 0.1,
        smoothing: 0,
    };
    let descriptor = rdf_descriptor(&nacl(), &options).unwrap();
    assert_eq!(descriptor.r.len(), 80);
    assert_eq!(descriptor.total.len(), 80);
    assert_eq!(descriptor.partial_names(), ["Cl-Cl", "Cl-Na", "Na-Na"]);
    assert_eq!(descriptor.coordination.len(), 80);
    assert!(descriptor.density > 0.0);

    let first_nonzero = descriptor
        .total
        .iter()
        .position(|value| *value > 0.0)
        .unwrap();
    assert!(
        (descriptor.r[first_nonzero] - 2.82).abs() < 0.15,
        "first non-zero bin at {}",
        descriptor.r[first_nonzero]
    );
    assert_eq!(descriptor.to_vector(), descriptor.total);
    assert_eq!(descriptor.to_vector_with_partials().len(), 80 + 3 * 80);

    let json = descriptor.to_json().unwrap();
    let parsed: xtalforge_analysis::RdfDescriptor = serde_json::from_str(&json).unwrap();
    assert_eq!(parsed.partial_names(), descriptor.partial_names());
    assert_eq!(parsed.r.len(), descriptor.r.len());
    assert_close(&parsed.total, &descriptor.total, 1e-12);
    assert_close(&parsed.coordination, &descriptor.coordination, 1e-12);
    assert!((parsed.density - descriptor.density).abs() < 1e-12);
    for key in descriptor.partial.keys() {
        assert_close(&parsed.partial[key], &descriptor.partial[key], 1e-12);
    }

    assert!(rdf_descriptor(&empty(), &options).is_err());
    assert!(
        rdf_descriptor(
            &nacl(),
            &RdfOptions {
                r_max: -1.0,
                ..RdfOptions::default()
            }
        )
        .is_err()
    );
}

#[test]
fn the_coulomb_matrix_is_symmetric_and_uses_minimum_image_distances() {
    let matrix = coulomb_matrix(&nacl(), &CoulombMatrixOptions::default()).unwrap();
    assert_eq!(matrix.n_atoms(), 8);
    assert_eq!(matrix.matrix.len(), 8);
    assert_eq!(matrix.to_vector().len(), 64);
    assert_eq!(matrix.atomic_numbers[0], 11);
    assert_eq!(matrix.elements[0], "Na");
    assert_eq!(matrix.elements[4], "Cl");
    assert!((matrix.matrix[0][0] - 0.5 * 11.0_f64.powf(2.4)).abs() < 1e-12);
    assert!((matrix.matrix[4][4] - 0.5 * 17.0_f64.powf(2.4)).abs() < 1e-12);
    for row in 0..8 {
        for column in 0..8 {
            assert_eq!(matrix.matrix[row][column], matrix.matrix[column][row]);
        }
    }
    let nearest = matrix
        .positions
        .iter()
        .enumerate()
        .filter(|(index, _)| *index >= 4)
        .find(|(_, position)| **position == [0.5, 0.0, 0.0])
        .map(|(index, _)| index)
        .unwrap();
    let expected = 11.0 * 17.0 / 2.82;
    assert!(
        (matrix.matrix[0][nearest] - expected).abs() < 0.5,
        "{} vs {expected}",
        matrix.matrix[0][nearest]
    );

    let json = matrix.to_json().unwrap();
    let parsed: xtalforge_analysis::CoulombMatrix = serde_json::from_str(&json).unwrap();
    assert_eq!(parsed, matrix);
}

#[test]
fn a_single_atom_has_a_one_by_one_coulomb_matrix() {
    let matrix = coulomb_matrix(&polonium(), &CoulombMatrixOptions::default()).unwrap();
    assert_eq!(matrix.n_atoms(), 1);
    assert_eq!(matrix.matrix, vec![vec![0.5 * 84.0_f64.powf(2.4)]]);
}

#[test]
fn coulomb_matrix_error_paths() {
    assert!(coulomb_matrix(&empty(), &CoulombMatrixOptions::default()).is_err());
    assert!(coulomb_matrix(&nacl(), &CoulombMatrixOptions { max_atoms: 4 }).is_err());
    assert!(coulomb_matrix(&nacl(), &CoulombMatrixOptions { max_atoms: 0 }).is_err());
    let coincident = StructureBuilder::new("coincident", Lattice::cubic(5.0).unwrap())
        .site(Site::new("Na1", "Na", [0.0, 0.0, 0.0]).unwrap())
        .site(Site::new("Cl1", "Cl", [0.0, 0.0, 0.0]).unwrap())
        .build()
        .unwrap();
    assert!(coulomb_matrix(&coincident, &CoulombMatrixOptions::default()).is_err());
}

#[test]
fn the_soap_like_descriptor_has_the_documented_layout() {
    let descriptor = soap_like_descriptor(&cscl(), &SoapLikeOptions::default()).unwrap();
    assert_eq!(descriptor.species, ["Cl", "Cs"]);
    assert_eq!(descriptor.n_radial, 4);
    assert_eq!(descriptor.l_max, 2);
    assert_eq!(descriptor.radial_centers, [1.0, 2.0, 3.0, 4.0]);
    assert_eq!(descriptor.n_centers(), 2);
    assert_eq!(descriptor.feature_len(), 152);
    assert_eq!(descriptor.per_center.len(), 2);
    assert!(descriptor.per_center.iter().all(|row| row.len() == 152));
    assert!(descriptor.mean.iter().all(|value| value.is_finite()));
    assert!(descriptor.mean.iter().any(|value| value.abs() > 1e-12));
    assert_eq!(descriptor.to_vector(), descriptor.mean);
    assert_eq!(descriptor.center_elements, ["Cl", "Cs"]);

    let json = descriptor.to_json().unwrap();
    let parsed: xtalforge_analysis::SoapLikeDescriptor = serde_json::from_str(&json).unwrap();
    assert_eq!(parsed.species, descriptor.species);
    assert_eq!(parsed.radial_centers, descriptor.radial_centers);
    assert_eq!(parsed.feature_len(), descriptor.feature_len());
    assert_eq!(parsed.center_elements, descriptor.center_elements);
    assert_eq!(parsed.per_center.len(), descriptor.per_center.len());
    for (row, expected) in parsed.per_center.iter().zip(descriptor.per_center.iter()) {
        assert_close(row, expected, 1e-12);
    }
    assert_close(&parsed.mean, &descriptor.mean, 1e-12);
}

#[test]
fn the_soap_like_descriptor_is_order_independent() {
    let first = soap_like_descriptor(&cscl(), &SoapLikeOptions::default()).unwrap();
    let second = soap_like_descriptor(&reversed_cscl(), &SoapLikeOptions::default()).unwrap();
    assert_eq!(first, second);
}

#[test]
fn a_cutoff_below_the_nearest_neighbour_gives_zeros() {
    let options = SoapLikeOptions {
        cutoff: 1.0,
        ..SoapLikeOptions::default()
    };
    let descriptor = soap_like_descriptor(&cscl(), &options).unwrap();
    assert!(descriptor.mean.iter().all(|value| *value == 0.0));
    assert!(
        descriptor
            .per_center
            .iter()
            .all(|row| row.iter().all(|value| *value == 0.0))
    );
}

#[test]
fn periodic_self_images_are_neighbours() {
    let descriptor = soap_like_descriptor(&polonium(), &SoapLikeOptions::default()).unwrap();
    assert_eq!(descriptor.n_centers(), 1);
    assert_eq!(descriptor.feature_len(), 52);
    assert!(descriptor.mean.iter().any(|value| *value > 0.0));
}

#[test]
fn soap_like_error_paths() {
    assert!(soap_like_descriptor(&empty(), &SoapLikeOptions::default()).is_err());
    for options in [
        SoapLikeOptions {
            cutoff: 0.0,
            ..SoapLikeOptions::default()
        },
        SoapLikeOptions {
            n_radial: 0,
            ..SoapLikeOptions::default()
        },
        SoapLikeOptions {
            sigma: 0.0,
            ..SoapLikeOptions::default()
        },
        SoapLikeOptions {
            l_max: MAX_SOAP_L + 1,
            ..SoapLikeOptions::default()
        },
        SoapLikeOptions {
            max_neighbors: 0,
            ..SoapLikeOptions::default()
        },
    ] {
        assert!(soap_like_descriptor(&cscl(), &options).is_err());
    }
}
