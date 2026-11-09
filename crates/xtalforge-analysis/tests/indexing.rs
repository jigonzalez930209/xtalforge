//! Integration tests for powder indexing (Task 7.5).

#![allow(
    clippy::unwrap_used,
    clippy::expect_used,
    clippy::panic,
    clippy::indexing_slicing
)]

use std::time::Instant;

use xtalforge_analysis::{
    CellMetric, LatticeSystem, Peak, figures_of_merit, index_cubic, index_pattern,
};

const CU_KA1: f64 = 1.540_598;

fn peaks_from_cell(cell: &CellMetric, max_index: i32, max_two_theta: f64) -> Vec<Peak> {
    let mut positions = Vec::new();
    let mut seed = 0x1234_5678_u64;
    let mut jitter = || {
        seed = seed
            .wrapping_mul(6_364_136_223_846_793_005)
            .wrapping_add(1_442_695_040_888_963_407);
        (((seed >> 11) as f64) / ((1_u64 << 53) as f64) - 0.5) * 0.006
    };
    for h in 0..=max_index {
        for k in 0..=max_index {
            for l in 0..=max_index {
                let Some(q) = cell.reciprocal_squared([h, k, l]) else {
                    continue;
                };
                let sine = (q.sqrt() * CU_KA1 / 2.0).clamp(0.0, 1.0);
                let two_theta = 2.0 * sine.asin().to_degrees();
                if two_theta > 5.0 && two_theta <= max_two_theta {
                    positions.push(two_theta);
                }
            }
        }
    }
    positions.sort_by(f64::total_cmp);
    positions.dedup_by(|first, second| (*first - *second).abs() < 0.01);
    positions.truncate(60);
    positions
        .into_iter()
        .map(|two_theta| Peak {
            two_theta: two_theta + jitter(),
            intensity: 1.0,
            fwhm: 0.1,
        })
        .collect()
}

fn check(cell: CellMetric, system: LatticeSystem, max_index: i32) {
    let peaks = peaks_from_cell(&cell, 10, 100.0);
    let start = Instant::now();
    let options = xtalforge_analysis::IndexingOptions {
        max_index,
        tolerance: 0.03,
        ..Default::default()
    };
    let result = index_pattern(&peaks, CU_KA1, &options).unwrap();
    let elapsed = start.elapsed();
    println!(
        "{system:?} peaks={} elapsed={elapsed:?} candidates={}",
        peaks.len(),
        result.candidates.len()
    );
    for candidate in result.candidates.iter().take(5) {
        println!(
            "  {:?} a={:.4} b={:.4} c={:.4} beta={:.3} m20={:.1} fn={:.1} completeness={:.2}",
            candidate.system,
            candidate.cell.a,
            candidate.cell.b,
            candidate.cell.c,
            candidate.cell.beta,
            candidate.m20,
            candidate.f_n,
            candidate.completeness
        );
    }
    let best = result.best().expect("a candidate");
    assert_eq!(best.system, system);
    assert!(best.completeness >= 0.95);
    // The same lattice can be written in several settings; any of them that
    // matches the expected lengths (and angle, for the monoclinic case) is the
    // right answer.
    let matches = |cell_parameters: &[f64; 3]| {
        let mut want = [cell.a, cell.b, cell.c];
        let mut got = *cell_parameters;
        want.sort_by(f64::total_cmp);
        got.sort_by(f64::total_cmp);
        got.iter()
            .zip(want.iter())
            .all(|(got, want)| (got - want).abs() < 0.05)
    };
    let candidate = result
        .candidates
        .iter()
        .find(|candidate| {
            candidate.system == system
                && matches(&[candidate.cell.a, candidate.cell.b, candidate.cell.c])
                && (system != LatticeSystem::Monoclinic || {
                    let mut got = [
                        candidate.cell.alpha,
                        candidate.cell.beta,
                        candidate.cell.gamma,
                    ];
                    let mut want = [cell.alpha, cell.beta, cell.gamma];
                    got.sort_by(f64::total_cmp);
                    want.sort_by(f64::total_cmp);
                    got.iter()
                        .zip(want.iter())
                        .all(|(got, want)| (got - want).abs() < 0.1)
                })
        })
        .unwrap_or_else(|| {
            panic!(
                "no {system:?} candidate matches {cell:?}; got {:?}",
                result
                    .candidates
                    .iter()
                    .take(5)
                    .map(|candidate| candidate.cell)
                    .collect::<Vec<_>>()
            )
        });
    assert!(candidate.completeness >= 0.95);
}

#[test]
fn indexes_a_cubic_cell_generally() {
    check(
        CellMetric::new(5.430941, 5.430941, 5.430941, 90.0, 90.0, 90.0),
        LatticeSystem::Cubic,
        3,
    );
}

#[test]
fn indexes_a_tetragonal_cell() {
    check(
        CellMetric::new(4.5937, 4.5937, 2.9587, 90.0, 90.0, 90.0),
        LatticeSystem::Tetragonal,
        3,
    );
}

#[test]
fn indexes_a_hexagonal_cell() {
    check(
        CellMetric::new(3.2417, 3.2417, 5.1876, 90.0, 90.0, 120.0),
        LatticeSystem::Hexagonal,
        3,
    );
}

#[test]
fn indexes_an_orthorhombic_cell() {
    check(
        CellMetric::new(20.0511, 19.8757, 13.36823, 90.0, 90.0, 90.0),
        LatticeSystem::Orthorhombic,
        3,
    );
}

#[test]
fn indexes_a_monoclinic_cell() {
    check(
        CellMetric::new(5.0, 6.0, 7.0, 90.0, 100.0, 90.0),
        LatticeSystem::Monoclinic,
        3,
    );
}

#[test]
fn figures_of_merit_are_finite_and_prefer_the_true_cell() {
    let cell = CellMetric::new(4.5937, 4.5937, 2.9587, 90.0, 90.0, 90.0);
    let peaks = peaks_from_cell(&cell, 8, 90.0);
    let good = figures_of_merit(&cell, &peaks, CU_KA1, 0.05).unwrap();
    let wrong = figures_of_merit(
        &CellMetric::new(4.7, 4.7, 3.0, 90.0, 90.0, 90.0),
        &peaks,
        CU_KA1,
        0.05,
    );
    assert!(good.m20 > 0.0 && good.f_n > 0.0);
    if let Ok(wrong) = wrong {
        assert!(good.m20 > wrong.m20);
    }
}

#[test]
fn cubic_compatibility_path_still_works() {
    let a = 5.430941;
    let hkls = [[1, 1, 1], [2, 2, 0], [3, 1, 1], [4, 0, 0], [3, 3, 1]];
    let peaks: Vec<Peak> = hkls
        .iter()
        .map(|hkl| {
            let sum = (hkl[0] * hkl[0] + hkl[1] * hkl[1] + hkl[2] * hkl[2]) as f64;
            let d = a / sum.sqrt();
            let two_theta = 2.0 * (CU_KA1 / (2.0 * d)).asin().to_degrees();
            Peak {
                two_theta,
                intensity: 1.0,
                fwhm: 0.1,
            }
        })
        .collect();
    let index = index_cubic(&peaks, CU_KA1, 5, 0.05).unwrap().unwrap();
    assert!((index.a - a).abs() < 0.05);
    assert!(index.completeness >= 0.8);
    assert!(index.mean_error < 0.05);
}
