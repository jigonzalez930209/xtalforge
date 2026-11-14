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
                g4: vec![[0.5, 0.0, 1.0]],
                ..AcsfOptions::default()
            }
        )
        .is_err()
    );
}

#[test]
fn acsf_matches_the_dscribe_fixture() {
    let fixture: serde_json::Value =
        serde_json::from_str(include_str!("../data/acsf_dscribe.json")).expect("fixture");
    let lattice =
        Lattice::new([[20.0, 0.0, 0.0], [0.0, 20.0, 0.0], [0.0, 0.0, 20.0]]).expect("cell");
    let mut builder = StructureBuilder::new("acsf", lattice);
    for (index, position) in fixture["structure"]["cartesian"]
        .as_array()
        .expect("positions")
        .iter()
        .enumerate()
    {
        let values = position.as_array().expect("position");
        let cart = [
            values[0].as_f64().expect("x"),
            values[1].as_f64().expect("y"),
            values[2].as_f64().expect("z"),
        ];
        builder = builder
            .site(Site::new(format!("C{}", index + 1), "C", lattice.to_frac(cart)).expect("site"));
    }
    let structure = builder.build().expect("structure");
    let options = AcsfOptions {
        cutoff: 4.0,
        g2: vec![[0.5, 0.0]],
        g4: vec![[0.5, 1.0, 1.0]],
        ..AcsfOptions::default()
    };
    let descriptor = acsf(&structure, &options).expect("acsf");
    let expected: Vec<Vec<f64>> = fixture["rows"]
        .as_array()
        .expect("rows")
        .iter()
        .map(|row| {
            row.as_array()
                .expect("row")
                .iter()
                .map(|value| value.as_f64().expect("number"))
                .collect()
        })
        .collect();
    let mut ours = descriptor.rows.clone();
    ours.sort_by(|a, b| a[0].total_cmp(&b[0]));
    let mut reference = expected;
    reference.sort_by(|a, b| a[0].total_cmp(&b[0]));
    for (row, expected_row) in ours.iter().zip(&reference) {
        for (value, expected) in row.iter().zip(expected_row) {
            assert!(
                (value - expected).abs() < 5e-6 * expected.abs().max(1.0),
                "ACSF feature {value} differs from DScribe {expected}"
            );
        }
    }
}

#[test]
fn mbtr_matches_the_dscribe_fixture() {
    let fixture: serde_json::Value =
        serde_json::from_str(include_str!("../data/mbtr_dscribe.json")).expect("fixture");
    let lattice =
        Lattice::new([[20.0, 0.0, 0.0], [0.0, 20.0, 0.0], [0.0, 0.0, 20.0]]).expect("cell");
    let mut builder = StructureBuilder::new("mbtr", lattice);
    for (index, position) in fixture["structure"]["cartesian"]
        .as_array()
        .expect("positions")
        .iter()
        .enumerate()
    {
        let values = position.as_array().expect("position");
        let cart = [
            values[0].as_f64().expect("x"),
            values[1].as_f64().expect("y"),
            values[2].as_f64().expect("z"),
        ];
        builder = builder
            .site(Site::new(format!("C{}", index + 1), "C", lattice.to_frac(cart)).expect("site"));
    }
    let structure = builder.build().expect("structure");
    let grid = |value: &serde_json::Value| {
        xtalforge_ml::MbtrGrid::new(
            value["min"].as_f64().expect("min"),
            value["max"].as_f64().expect("max"),
            value["n"].as_u64().expect("n") as usize,
            value["sigma"].as_f64().expect("sigma"),
        )
    };
    for (name, key) in [("k2", "k2"), ("k3", "k3")] {
        let options = xtalforge_ml::MbtrOptions {
            k2: if name == "k2" {
                Some(grid(&fixture[key]["grid"]))
            } else {
                None
            },
            k3: if name == "k3" {
                Some(grid(&fixture[key]["grid"]))
            } else {
                None
            },
            cutoff: 5.0,
            elements: vec!["C".to_owned()],
        };
        let descriptor = mbtr(&structure, &options).expect("mbtr");
        let expected: Vec<f64> = fixture[key]["features"]
            .as_array()
            .expect("features")
            .iter()
            .map(|value| value.as_f64().expect("number"))
            .collect();
        assert_eq!(descriptor.features.len(), expected.len(), "{name} length");
        for (value, expected) in descriptor.features.iter().zip(&expected) {
            assert!(
                (value - expected).abs() < 5e-6 * expected.abs().max(1.0),
                "MBTR {name} feature {value} differs from DScribe {expected}"
            );
        }
    }
}

#[test]
fn mbtr_shape_and_errors() {
    let structure = nacl();
    let options = MbtrOptions {
        k2: Some(xtalforge_ml::MbtrGrid::new(0.0, 6.0, 16, 0.1)),
        k3: Some(xtalforge_ml::MbtrGrid::new(0.0, 180.0, 16, 0.1)),
        cutoff: 6.0,
        elements: Vec::new(),
    };
    let descriptor = mbtr(&structure, &options).expect("mbtr");
    let elements = descriptor.elements.len();
    let pairs = elements * (elements + 1) / 2;
    let triples = elements * elements * (elements + 1) / 2;
    assert_eq!(descriptor.features.len(), (pairs + triples) * 16);
    assert!(descriptor.features.iter().all(|value| value.is_finite()));
    assert!(
        mbtr(
            &structure,
            &MbtrOptions {
                k2: None,
                k3: None,
                ..options
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
