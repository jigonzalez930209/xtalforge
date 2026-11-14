//! ACSF, MBTR and MOF feature tests (Tasks 21.1.3 and 21.1.4).

#![allow(
    clippy::unwrap_used,
    clippy::expect_used,
    clippy::panic,
    clippy::indexing_slicing
)]

use xtalforge_core::{Lattice, Site, StructureBuilder};
use xtalforge_ml::{
    AcsfOptions, MbtrOptions, MofFeaturesOptions, acsf, mbtr, mof_features, pore_geometry_vector,
};

fn nacl() -> xtalforge_core::Structure {
    let lattice = Lattice::from_parameters(5.64, 5.64, 5.64, 90.0, 90.0, 90.0).expect("cell");
    StructureBuilder::new("NaCl", lattice)
        .site(Site::new("Na1", "Na", [0.0, 0.0, 0.0]).expect("Na"))
        .site(Site::new("Cl1", "Cl", [0.5, 0.5, 0.5]).expect("Cl"))
        .build()
        .expect("structure")
}

#[test]
fn acsf_shape_and_average() {
    let structure = nacl();
    let descriptor = acsf(&structure, &AcsfOptions::default()).expect("acsf");
    assert_eq!(descriptor.rows.len(), 2);
    assert!(descriptor.n_features > 0);
    assert!(
        descriptor
            .rows
            .iter()
            .flatten()
            .all(|value| value.is_finite())
    );
    let averaged = acsf(
        &structure,
        &AcsfOptions {
            average: true,
            ..AcsfOptions::default()
        },
    )
    .expect("average");
    assert_eq!(averaged.rows.len(), 1);
    assert_eq!(averaged.n_features, descriptor.n_features);
}

#[test]
fn acsf_rejects_bad_options() {
    let structure = nacl();
    assert!(
        acsf(
            &structure,
            &AcsfOptions {
                cutoff: 0.0,
                ..AcsfOptions::default()
            }
        )
        .is_err()
    );
    assert!(
        acsf(
            &structure,
            &AcsfOptions {
                eta: Vec::new(),
                ..AcsfOptions::default()
            }
        )
        .is_err()
    );
}

#[test]
fn mbtr_shape_and_finiteness() {
    let structure = nacl();
    let options = MbtrOptions {
        cutoff: 6.0,
        n_grid: 16,
        ..MbtrOptions::default()
    };
    let descriptor = mbtr(&structure, &options).expect("mbtr");
    let elements = descriptor.elements.len();
    let pairs = elements * (elements + 1) / 2;
    let triples = elements * elements * elements;
    assert_eq!(
        descriptor.features.len(),
        (pairs + triples) * options.n_grid
    );
    assert!(descriptor.features.iter().all(|value| value.is_finite()));
    let no_angles = mbtr(
        &structure,
        &MbtrOptions {
            angular: false,
            ..options
        },
    )
    .expect("radial");
    assert_eq!(no_angles.features.len(), pairs * 16);
}

#[test]
fn mbtr_rejects_bad_options() {
    let structure = nacl();
    assert!(
        mbtr(
            &structure,
            &MbtrOptions {
                r_max: 0.0,
                ..MbtrOptions::default()
            }
        )
        .is_err()
    );
    assert!(
        mbtr(
            &structure,
            &MbtrOptions {
                n_grid: 1,
                ..MbtrOptions::default()
            }
        )
        .is_err()
    );
}

#[test]
fn mof_racs_and_pore_vector() {
    let features = mof_features(
        &nacl(),
        &MofFeaturesOptions {
            pores: false,
            ..MofFeaturesOptions::default()
        },
    )
    .expect("features");
    assert_eq!(features.properties.len(), 12);
    assert_eq!(features.features.len(), 12);
    assert!(features.pore_vector.is_empty());
    let pores = xtalforge_mof::platform::porosity::PoreDescriptors {
        pld: 1.0,
        lcd: 2.0,
        asa_m2_g: 3.0,
        asa_m2_cm3: 4.0,
        accessible_volume: 5.0,
        void_fraction: 0.6,
        channel_dimension: 3,
        density: 0.7,
        probe_radius: 1.2,
        spacing: 0.8,
    };
    let (values, names) = pore_geometry_vector(&pores);
    assert_eq!(values.len(), 9);
    assert_eq!(names.len(), 9);
    assert_eq!(values[0], 1.0);
}
