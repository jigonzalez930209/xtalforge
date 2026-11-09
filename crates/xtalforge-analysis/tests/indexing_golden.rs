//! Indexing and figures of merit against the Phase 7 golden data (Task 7.5).

#![allow(
    clippy::unwrap_used,
    clippy::expect_used,
    clippy::panic,
    clippy::indexing_slicing
)]

use serde_json::Value;
use xtalforge_analysis::{
    CellMetric, IndexingOptions, LatticeSystem, Peak, figures_of_merit, index_pattern,
};

const GOLDEN: &str = include_str!("../../../tests/golden/phase-7/indexing/cells.json");

/// Two cells describe the same lattice when their sorted `Q` spectra agree.
fn same_lattice(first: &CellMetric, second: &CellMetric) -> bool {
    let spectrum = |cell: &CellMetric| {
        let mut values = Vec::new();
        for h in 0..=5 {
            for k in 0..=5 {
                for l in 0..=5 {
                    if let Some(q) = cell.reciprocal_squared([h, k, l]) {
                        values.push(q);
                    }
                }
            }
        }
        values.sort_by(f64::total_cmp);
        values.dedup_by(|a, b| (*a - *b).abs() < 1e-9);
        values
    };
    let first = spectrum(first);
    let second = spectrum(second);
    first.len() == second.len()
        && first
            .iter()
            .zip(second.iter())
            .all(|(a, b)| (a - b).abs() <= 5e-4 * a.abs().max(1e-6))
}

fn system(name: &str) -> LatticeSystem {
    match name {
        "cubic" => LatticeSystem::Cubic,
        "tetragonal" => LatticeSystem::Tetragonal,
        "hexagonal" => LatticeSystem::Hexagonal,
        "orthorhombic" => LatticeSystem::Orthorhombic,
        "monoclinic" => LatticeSystem::Monoclinic,
        other => panic!("unknown golden system {other}"),
    }
}

#[test]
fn figures_of_merit_match_the_python_reference() {
    let golden: Value = serde_json::from_str(GOLDEN).unwrap();
    let mut checked = 0;
    for entry in golden["cells"].as_array().unwrap() {
        let cell = &entry["cell"];
        let value = |key: &str| cell[key].as_f64().unwrap();
        let cell = CellMetric::new(
            value("a"),
            value("b"),
            value("c"),
            value("alpha"),
            value("beta"),
            value("gamma"),
        );
        let wavelength = entry["wavelength"].as_f64().unwrap();
        let tolerance = entry["tolerance"].as_f64().unwrap();
        let peaks: Vec<Peak> = entry["peaks"]
            .as_array()
            .unwrap()
            .iter()
            .map(|value| Peak {
                two_theta: value.as_f64().unwrap(),
                intensity: 1.0,
                fwhm: 0.1,
            })
            .collect();
        let found = figures_of_merit(&cell, &peaks, wavelength, tolerance).unwrap();
        let expected = &entry["figures"];
        let name = entry["name"].as_str().unwrap();
        for key in ["m20", "f_n", "n20", "n_possible"] {
            let want = expected[key].as_f64().unwrap();
            let got = match key {
                "m20" => found.m20,
                "f_n" => found.f_n,
                "n20" => found.n20 as f64,
                _ => found.n_possible as f64,
            };
            let tolerance = if key == "m20" || key == "f_n" {
                1e-4 * want.abs().max(1.0)
            } else {
                1e-6
            };
            assert!(
                (got - want).abs() <= tolerance,
                "{name}: {key} {got} != {want}"
            );
        }
        checked += 1;
    }
    assert_eq!(checked, 5);
}

#[test]
fn index_pattern_recovers_every_golden_cell() {
    let golden: Value = serde_json::from_str(GOLDEN).unwrap();
    let mut checked = 0;
    for entry in golden["cells"].as_array().unwrap() {
        let cell = &entry["cell"];
        let value = |key: &str| cell[key].as_f64().unwrap();
        let true_cell = CellMetric::new(
            value("a"),
            value("b"),
            value("c"),
            value("alpha"),
            value("beta"),
            value("gamma"),
        );
        let name = entry["name"].as_str().unwrap();
        let wavelength = entry["wavelength"].as_f64().unwrap();
        let tolerance = entry["tolerance"].as_f64().unwrap();
        let peaks: Vec<Peak> = entry["peaks"]
            .as_array()
            .unwrap()
            .iter()
            .map(|value| Peak {
                two_theta: value.as_f64().unwrap(),
                intensity: 1.0,
                fwhm: 0.1,
            })
            .collect();
        let options = IndexingOptions {
            max_index: 3,
            tolerance,
            ..IndexingOptions::default()
        };
        let result = index_pattern(&peaks, wavelength, &options).unwrap();
        let expected_system = system(entry["system"].as_str().unwrap());
        assert_eq!(
            result.best().map(|candidate| candidate.system),
            Some(expected_system),
            "{name}: wrong best family; candidates {:?}",
            result
                .candidates
                .iter()
                .take(5)
                .map(|candidate| (candidate.system, candidate.cell))
                .collect::<Vec<_>>()
        );
        // A supercell of the same family can index the same lines; the
        // expected lattice must appear among the candidates.
        let best = result
            .candidates
            .iter()
            .find(|candidate| same_lattice(&true_cell, &candidate.cell))
            .unwrap_or_else(|| {
                panic!(
                    "{name}: expected lattice not among {:?}",
                    result
                        .candidates
                        .iter()
                        .take(5)
                        .map(|candidate| candidate.cell)
                        .collect::<Vec<_>>()
                )
            });
        assert!(best.completeness >= 0.95, "{name}: coverage");
        // The true-cell figures are checked exactly by the test above; here a
        // recovered cell only has to be a good fit (a slightly different cell
        // setting can move its M20 by a factor).
        let want = entry["figures"]["m20"].as_f64().unwrap();
        assert!(
            best.m20 >= 0.1 * want,
            "{name}: M20 {} is far below {want}",
            best.m20
        );
        checked += 1;
    }
    assert_eq!(checked, 5);
}

#[test]
fn a_wrong_cell_scores_below_the_golden_cell() {
    let golden: Value = serde_json::from_str(GOLDEN).unwrap();
    let entry = &golden["cells"][0];
    let cell = &entry["cell"];
    let value = |key: &str| cell[key].as_f64().unwrap();
    let true_cell = CellMetric::new(
        value("a"),
        value("b"),
        value("c"),
        value("alpha"),
        value("beta"),
        value("gamma"),
    );
    let peaks: Vec<Peak> = entry["peaks"]
        .as_array()
        .unwrap()
        .iter()
        .map(|value| Peak {
            two_theta: value.as_f64().unwrap(),
            intensity: 1.0,
            fwhm: 0.1,
        })
        .collect();
    let good = figures_of_merit(&true_cell, &peaks, 1.540_598, 0.03).unwrap();
    let wrong = figures_of_merit(
        &CellMetric::new(4.0, 4.0, 4.0, 90.0, 90.0, 90.0),
        &peaks,
        1.540_598,
        0.03,
    )
    .map(|figures| figures.m20)
    .unwrap_or(0.0);
    assert!(good.m20 > wrong);
}
