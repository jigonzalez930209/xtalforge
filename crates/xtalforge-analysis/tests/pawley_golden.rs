//! Pawley fit against the Phase 7 golden pattern (Task 7.6.2).

#![allow(
    clippy::unwrap_used,
    clippy::expect_used,
    clippy::panic,
    clippy::indexing_slicing
)]

use std::collections::BTreeMap;

use serde_json::Value;
use xtalforge_analysis::{PawleyOptions, ProfileKind, pawley_refine};
use xtalforge_core::{Lattice, Site, StructureBuilder};
use xtalforge_io::xy::PatternPoint;

const GOLDEN: &str = include_str!("../../../tests/golden/phase-7/pawley/one-atom-cubic.json");

#[test]
fn pawley_recovers_the_golden_pattern() {
    let golden: Value = serde_json::from_str(GOLDEN).unwrap();
    let cell = golden["cell"].as_array().unwrap();
    let value = |index: usize| cell[index].as_f64().unwrap();
    let structure = StructureBuilder::new(
        "one-atom",
        Lattice::from_parameters(value(0), value(1), value(2), value(3), value(4), value(5))
            .unwrap(),
    )
    .site(Site::new("C1", "C", [0.0, 0.0, 0.0]).unwrap())
    .build()
    .unwrap();
    let peaks: Vec<PatternPoint> = golden["peaks"]
        .as_array()
        .unwrap()
        .iter()
        .map(|peak| PatternPoint {
            two_theta: peak[0].as_f64().unwrap(),
            intensity: peak[1].as_f64().unwrap(),
            uncertainty: None,
        })
        .collect();
    let options = PawleyOptions {
        kind: ProfileKind::Gaussian,
        fwhm: golden["fwhm"].as_f64().unwrap(),
        background_order: 1,
        refine_cell: true,
        refine_zero: true,
        max_iterations: 30,
        ..PawleyOptions::default()
    };
    let result = pawley_refine(&structure, &peaks, &options).unwrap();
    println!(
        "reflections {} rwp {} rp {} zero {} a {} conv {}",
        result.reflections.len(),
        result.rwp,
        result.rp,
        result.zero_shift,
        result.cell[0],
        result.converged
    );
    println!("cell {:?}", result.cell);
    for reflection in result.reflections.iter().take(6) {
        println!(
            "  {:?} {} {}",
            reflection.hkl, reflection.two_theta, reflection.intensity
        );
    }
    let expected_zero = golden["zero_shift"].as_f64().unwrap();
    assert!(
        (result.zero_shift - expected_zero).abs() < 5e-3,
        "zero shift {} != {expected_zero}",
        result.zero_shift
    );
    assert!(
        (result.cell[0] - value(0)).abs() < 2e-3,
        "a {} != {}",
        result.cell[0],
        value(0)
    );
    assert!(result.rwp < 0.02, "Rwp {}", result.rwp);
    assert!(result.rp < 0.02, "Rp {}", result.rp);
    assert!(result.chi_squared < 0.02, "chi2 {}", result.chi_squared);
    assert!(result.converged);

    let found: BTreeMap<[i32; 3], (f64, f64)> = result
        .reflections
        .iter()
        .map(|reflection| (reflection.hkl, (reflection.intensity, reflection.esd)))
        .collect();
    let mut checked = 0;
    for reflection in golden["reflections"].as_array().unwrap() {
        let hkl = [
            reflection["hkl"][0].as_i64().unwrap() as i32,
            reflection["hkl"][1].as_i64().unwrap() as i32,
            reflection["hkl"][2].as_i64().unwrap() as i32,
        ];
        let expected = reflection["intensity"].as_f64().unwrap();
        let (intensity, esd) = found
            .get(&hkl)
            .copied()
            .unwrap_or_else(|| panic!("missing reflection {hkl:?}"));
        assert!(
            (intensity - expected).abs() <= 3.0 * esd.max(1e-6),
            "{hkl:?}: {intensity} != {expected} (esd {esd})"
        );
        checked += 1;
    }
    assert_eq!(checked, golden["reflections"].as_array().unwrap().len());
}
