//! Integration tests for EOS fits and convex hulls (Task 9.6).

#![allow(
    clippy::unwrap_used,
    clippy::expect_used,
    clippy::panic,
    clippy::indexing_slicing
)]

use std::collections::BTreeMap;

use xtalforge_core::{EosKind, EosPoint, PhaseEntry, convex_hull, eos_energy, fit_eos};

#[test]
fn birch_murnaghan_fit_recovers_the_synthetic_parameters() {
    let truth = [20.0, -10.0, 0.5, 4.0];
    let mut points = Vec::new();
    for step in 0..11 {
        let volume = 18.0 + 0.4 * f64::from(step);
        points.push(EosPoint {
            volume,
            energy: eos_energy(EosKind::BirchMurnaghan, &truth, volume),
        });
    }
    let fit = fit_eos(&points, EosKind::BirchMurnaghan).unwrap();
    assert!((fit.v0 - truth[0]).abs() < 1e-3, "v0 {}", fit.v0);
    assert!((fit.e0 - truth[1]).abs() < 1e-4, "e0 {}", fit.e0);
    assert!((fit.b0 - truth[2] * 160.217_662_08).abs() < 0.01);
    assert!((fit.b0_prime - truth[3]).abs() < 1e-3);
    assert!(fit.rms < 1e-8);
}

#[test]
fn all_three_eos_models_fit() {
    let truth = [20.0, -10.0, 0.5, 4.0];
    for kind in [EosKind::BirchMurnaghan, EosKind::Vinet, EosKind::Murnaghan] {
        let points: Vec<EosPoint> = (0..11)
            .map(|step| {
                let volume = 18.0 + 0.4 * f64::from(step);
                EosPoint {
                    volume,
                    energy: eos_energy(kind, &truth, volume),
                }
            })
            .collect();
        let fit = fit_eos(&points, kind).unwrap();
        assert!((fit.v0 - 20.0).abs() < 0.01, "{kind:?} v0 {}", fit.v0);
        assert!(fit.rms < 1e-7, "{kind:?} rms {}", fit.rms);
    }
}

#[test]
fn eos_fit_requires_enough_points() {
    let points = [
        EosPoint {
            volume: 1.0,
            energy: 0.0,
        },
        EosPoint {
            volume: 2.0,
            energy: 1.0,
        },
    ];
    assert!(fit_eos(&points, EosKind::BirchMurnaghan).is_err());
}

fn entry(name: &str, a: f64, b: f64, energy: f64) -> PhaseEntry {
    let mut composition = BTreeMap::new();
    composition.insert("A".to_owned(), a);
    composition.insert("B".to_owned(), b);
    PhaseEntry {
        name: name.to_owned(),
        composition,
        energy,
    }
}

#[test]
fn convex_hull_marks_the_stable_and_unstable_phases() {
    let entries = [
        entry("A", 1.0, 0.0, 0.0),
        entry("B", 0.0, 1.0, 0.0),
        entry("AB", 1.0, 1.0, -1.0),
        entry("A2B", 2.0, 1.0, -1.2),
        entry("A3B", 3.0, 1.0, -0.8),
    ];
    let hull = convex_hull(&entries, None).unwrap();
    assert!(hull.points[0].stable);
    assert!(hull.points[1].stable);
    assert!(hull.points[2].stable);
    // A2B (x_B = 1/3) is below the A-AB tie line.
    assert!(hull.points[3].stable);
    // A3B (x_B = 0.25) lies above the A-A2B tie line.
    assert!(!hull.points[4].stable);
    assert!(hull.points[4].distance > 0.0);
    assert!(!hull.points[4].decomposition.is_empty());
}

#[test]
fn convex_hull_uses_supplied_references() {
    let entries = [
        entry("A", 1.0, 0.0, 0.0),
        entry("B", 0.0, 1.0, 0.0),
        entry("AB", 1.0, 1.0, -2.0),
    ];
    let mut references = BTreeMap::new();
    references.insert("A".to_owned(), -1.0);
    references.insert("B".to_owned(), -1.0);
    let hull = convex_hull(&entries, Some(&references)).unwrap();
    let ab = &hull.points[2];
    assert!((ab.formation_energy - 0.0).abs() < 1e-9);
    assert!(ab.stable);
}
