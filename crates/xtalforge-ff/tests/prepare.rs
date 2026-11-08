//! Integration tests for the preparation pipeline (Task 4.7).
//!
//! The water case exercises the default pipeline: FIRE, QEq charges and the
//! report fields.  Ewald and vibrations are checked as separate opt-in steps so
//! a failure names the step that broke rather than the whole workflow.

#![allow(
    clippy::unwrap_used,
    clippy::expect_used,
    clippy::panic,
    clippy::indexing_slicing
)]

use xtalforge_core::Structure;
use xtalforge_ff::{OptimizeOptions, PrepareOptions, PrepareStep, UffCalculator, prepare};

fn molecule(name: &str, atoms: &[(&str, [f64; 3])]) -> Structure {
    let points: Vec<(String, [f64; 3])> = atoms
        .iter()
        .map(|(element, point)| ((*element).to_owned(), *point))
        .collect();
    Structure::molecule(name, &points, 8.0).expect("molecule")
}

fn distorted_water() -> Structure {
    molecule(
        "distorted-water",
        &[
            ("O", [0.0, 0.0, 0.0]),
            ("H", [1.08, 0.0, 0.0]),
            ("H", [-0.27, 1.04, 0.0]),
        ],
    )
}

fn options() -> PrepareOptions {
    PrepareOptions {
        minimize: true,
        charges: true,
        ewald: false,
        vibrations: false,
        optimize: OptimizeOptions {
            max_steps: 300,
            force_tolerance: 0.05,
            ..OptimizeOptions::default()
        },
        ..PrepareOptions::default()
    }
}

fn names(steps: &[PrepareStep]) -> Vec<&str> {
    steps.iter().map(|step| step.name.as_str()).collect()
}

#[test]
fn prepare_relaxes_water_and_records_every_step() {
    let calculator = UffCalculator::with_defaults().unwrap();
    let report = prepare(&calculator, &distorted_water(), &options()).unwrap();
    assert_eq!(
        names(&report.steps),
        ["minimize", "charges", "ewald", "vibrations"]
    );
    assert!(report.steps[0].done && report.steps[1].done);
    assert!(!report.steps[2].done && !report.steps[3].done);
    assert!(report.energy.is_finite(), "energy {}", report.energy);
    assert!(report.total_charge.abs() < 1e-8, "{}", report.total_charge);
    assert!(report.ewald_energy.is_none());
    assert!(report.frequencies.is_empty());
    let optimized = report.optimized.expect("minimizer result");
    assert!(optimized.steps > 0);
}

#[test]
fn prepare_can_evaluate_ewald_on_the_assigned_charges() {
    let calculator = UffCalculator::with_defaults().unwrap();
    let mut settings = options();
    settings.minimize = false;
    settings.ewald = true;
    let report = prepare(&calculator, &distorted_water(), &settings).unwrap();
    assert!(report.steps[2].done);
    let energy = report.ewald_energy.expect("Ewald energy");
    assert!(energy.is_finite(), "Ewald energy {energy}");
    assert!(report.total_charge.abs() < 1e-8);
    assert!(report.optimized.is_none());
}

#[test]
fn prepare_can_skip_optional_steps() {
    let calculator = UffCalculator::with_defaults().unwrap();
    let settings = PrepareOptions {
        minimize: false,
        charges: false,
        ..PrepareOptions::default()
    };
    let report = prepare(&calculator, &distorted_water(), &settings).unwrap();
    assert!(report.steps.iter().all(|step| !step.done));
    assert!(report.energy.is_finite());
    assert!(report.optimized.is_none());
    assert!(report.ewald_energy.is_none());
    assert!(report.frequencies.is_empty());
}

#[test]
fn prepare_can_build_gamma_frequencies() {
    let calculator = UffCalculator::with_defaults().unwrap();
    let settings = PrepareOptions {
        minimize: false,
        charges: false,
        vibrations: true,
        ..PrepareOptions::default()
    };
    let report = prepare(&calculator, &distorted_water(), &settings).unwrap();
    assert!(report.steps[3].done);
    assert_eq!(report.frequencies.len(), 9);
    assert!(report.frequencies.iter().all(|value| value.is_finite()));
}

#[test]
fn prepare_rejects_invalid_options() {
    let calculator = UffCalculator::with_defaults().unwrap();
    let invalid_minimizer = PrepareOptions {
        optimize: OptimizeOptions {
            max_steps: 0,
            ..OptimizeOptions::default()
        },
        ..PrepareOptions::default()
    };
    assert!(prepare(&calculator, &distorted_water(), &invalid_minimizer).is_err());
    let invalid_charges = PrepareOptions {
        charge_options: xtalforge_ff::ChargeOptions {
            total_charge: f64::INFINITY,
            ..xtalforge_ff::ChargeOptions::default()
        },
        ..PrepareOptions::default()
    };
    assert!(prepare(&calculator, &distorted_water(), &invalid_charges).is_err());
}
