//! Exact SOAP tests against the recorded DScribe fixture (Task 21.1.2).

#![allow(
    clippy::unwrap_used,
    clippy::expect_used,
    clippy::panic,
    clippy::indexing_slicing
)]

use xtalforge_core::{Lattice, Site, StructureBuilder};
use xtalforge_ml::soap::{SoapOptions, soap_power_spectrum};

const FIXTURE: &str = include_str!("../data/soap_dscribe.json");

fn structure() -> xtalforge_core::Structure {
    let fixture: serde_json::Value = serde_json::from_str(FIXTURE).expect("fixture");
    let lattice = fixture["structure"]["lattice"]
        .as_array()
        .expect("lattice")
        .iter()
        .map(|row| {
            let values = row.as_array().expect("row");
            [
                values[0].as_f64().expect("a"),
                values[1].as_f64().expect("b"),
                values[2].as_f64().expect("c"),
            ]
        })
        .collect::<Vec<_>>();
    let lattice = Lattice::new([lattice[0], lattice[1], lattice[2]]).expect("cell");
    let mut builder = StructureBuilder::new("soap", lattice);
    for (index, (element, frac)) in fixture["structure"]["elements"]
        .as_array()
        .expect("elements")
        .iter()
        .zip(fixture["structure"]["frac"].as_array().expect("frac"))
        .enumerate()
    {
        let values = frac.as_array().expect("frac row");
        builder = builder.site(
            Site::new(
                format!("{}{}", element.as_str().expect("element"), index + 1),
                element.as_str().expect("element"),
                [
                    values[0].as_f64().expect("x"),
                    values[1].as_f64().expect("y"),
                    values[2].as_f64().expect("z"),
                ],
            )
            .expect("site"),
        );
    }
    builder.build().expect("structure")
}

#[test]
fn feature_layout_matches_the_dscribe_fixture() {
    let fixture: serde_json::Value = serde_json::from_str(FIXTURE).expect("fixture");
    let rows = fixture["rows"].as_array().expect("rows");
    assert_eq!(rows.len(), 2);
    assert_eq!(rows[0].as_array().expect("row").len(), 18);
    let descriptor = soap_power_spectrum(
        &structure(),
        &SoapOptions {
            cutoff: 4.0,
            n_max: 3,
            l_max: 2,
            sigma: 0.5,
            average: false,
        },
    )
    .expect("soap");
    assert_eq!(descriptor.n_features(), 18);
    assert_eq!(descriptor.layout[0], [0, 0, 0]);
    assert_eq!(descriptor.layout[17], [2, 2, 2]);
    let expected: Vec<Vec<f64>> = rows
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
            let tolerance = 5e-4 * expected.abs().max(1.0);
            assert!(
                (value - expected).abs() < tolerance,
                "SOAP feature {value} differs from DScribe {expected}"
            );
        }
    }
}

#[test]
fn invariants_and_errors() {
    let structure = structure();
    let options = SoapOptions {
        cutoff: 5.0,
        n_max: 2,
        l_max: 1,
        sigma: 0.7,
        average: true,
    };
    let first = soap_power_spectrum(&structure, &options).expect("first");
    assert_eq!(first.rows.len(), 1);
    assert_eq!(first.n_features(), 2 * 3 / 2 * 2);
    assert!(first.rows[0].iter().all(|value| value.is_finite()));
    let mut bad = options.clone();
    bad.cutoff = 0.5;
    assert!(soap_power_spectrum(&structure, &bad).is_err());
    bad = options;
    bad.l_max = 6;
    assert!(soap_power_spectrum(&structure, &bad).is_err());
}
