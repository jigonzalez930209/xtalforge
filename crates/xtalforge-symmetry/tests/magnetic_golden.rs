//! Magnetic space groups against spglib's recorded axial (is_axial=true) data.
//!
//! spglib's magnetic test files are not redistributed; the JSON in
//! `crates/xtalforge-symmetry/tests/fixtures/` stores only the numeric reference generated from
//! spglib 2.8.0, with the same cell, moments and tolerances.

#![allow(
    clippy::unwrap_used,
    clippy::expect_used,
    clippy::panic,
    clippy::indexing_slicing
)]

use serde::Deserialize;
use serde_json::json;
use xtalforge_core::elements::element_by_z;
use xtalforge_core::{Lattice, Site, StructureBuilder};
use xtalforge_symmetry::{MagneticSymmetryOptions, detect_magnetic_symmetry_with_options};

#[derive(Debug, Deserialize)]
struct Golden {
    cases: Vec<Case>,
}

#[derive(Debug, Deserialize)]
struct Case {
    name: String,
    lattice: [[f64; 3]; 3],
    positions: Vec<[f64; 3]>,
    numbers: Vec<u8>,
    moments: Vec<f64>,
    symprec: f64,
    mag_symprec: f64,
    uni_number: u32,
    n_operations: usize,
    n_time_reversal: usize,
}

#[test]
fn magnetic_space_groups_match_spglib_axial_data() {
    let golden: Golden =
        serde_json::from_str(include_str!("fixtures/magnetic-spglib.json")).expect("golden JSON");
    let mut compared = 0;
    for case in golden.cases {
        let lattice = Lattice::new(case.lattice).expect("lattice");
        let mut builder = StructureBuilder::new(case.name.clone(), lattice);
        for (index, (position, number)) in case.positions.iter().zip(&case.numbers).enumerate() {
            let element = element_by_z(*number).expect("element");
            let site = Site::new(
                format!("{}{}", element.symbol, index + 1),
                &element.symbol,
                *position,
            )
            .expect("site")
            .with_prop("magmom", json!(case.moments[index]));
            builder = builder.site(site);
        }
        let structure = builder.build().expect("structure");
        let options = MagneticSymmetryOptions {
            symprec: case.symprec,
            mag_symprec: Some(case.mag_symprec),
            rotate_basis: false,
        };
        let dataset = detect_magnetic_symmetry_with_options(&structure, &options)
            .unwrap_or_else(|error| panic!("{}: {error}", case.name));
        assert_eq!(
            dataset.uni_number, case.uni_number,
            "{} UNI number",
            case.name
        );
        assert_eq!(
            dataset.operations.len(),
            case.n_operations,
            "{} operation count",
            case.name
        );
        let time_reversal = dataset
            .operations
            .iter()
            .filter(|operation| operation.time_reversal)
            .count();
        assert_eq!(
            time_reversal, case.n_time_reversal,
            "{} time-reversal count",
            case.name
        );
        compared += 1;
    }
    assert_eq!(compared, 3, "all golden magnetic cases must run");
}
