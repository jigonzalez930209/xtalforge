//! Integration tests for molecular dynamics (Task 14.2.5).
//!
//! The systems are small on purpose: methane at its UFF minimum for the
//! thermostats and energy conservation, and fcc argon for the barostat.

#![allow(
    clippy::unwrap_used,
    clippy::expect_used,
    clippy::panic,
    clippy::indexing_slicing
)]

use xtalforge_core::{Lattice, Site, Structure, StructureBuilder};
use xtalforge_ff::uff::natural_bond_length;
use xtalforge_ff::{
    CancelToken, Ensemble, MdOptions, MdResult, ProgressEvent, RunControl, UffCalculator,
    molecular_dynamics, molecular_dynamics_with,
};

fn methane() -> Structure {
    let length = natural_bond_length("C_3", "H_", 1.0).unwrap();
    let a = length / 3.0_f64.sqrt();
    let points: Vec<(String, [f64; 3])> = [
        ("C", [0.0, 0.0, 0.0]),
        ("H", [a, a, a]),
        ("H", [a, -a, -a]),
        ("H", [-a, a, -a]),
        ("H", [-a, -a, a]),
    ]
    .iter()
    .map(|(element, point)| ((*element).to_owned(), *point))
    .collect();
    Structure::molecule("methane", &points, 12.0).expect("methane")
}

fn argon_fcc(edge: f64) -> Structure {
    let lattice = Lattice::cubic(edge).unwrap();
    let mut builder = StructureBuilder::new("argon-fcc", lattice);
    for (index, point) in [
        [0.0, 0.0, 0.0],
        [0.0, 0.5, 0.5],
        [0.5, 0.0, 0.5],
        [0.5, 0.5, 0.0],
    ]
    .iter()
    .enumerate()
    {
        builder = builder.site(Site::new(format!("Ar{index}"), "Ar", *point).unwrap());
    }
    builder.build().unwrap()
}

fn nve_options() -> MdOptions {
    MdOptions {
        ensemble: Ensemble::Nve,
        time_step_fs: 0.25,
        steps: 400,
        seed: 42,
        initial_temperature: 300.0,
        sample_every: 10,
        ..MdOptions::default()
    }
}

fn mean_volume_tail(result: &MdResult) -> f64 {
    let start = result.samples.len() / 2;
    let tail = &result.samples[start..];
    tail.iter().map(|sample| sample.volume).sum::<f64>() / tail.len() as f64
}

#[test]
fn nve_conserves_energy_and_is_deterministic() {
    let calculator = UffCalculator::with_defaults().unwrap();
    let structure = methane();
    let first = molecular_dynamics(&calculator, &structure, &nve_options()).unwrap();
    assert!(!first.cancelled);
    assert_eq!(first.steps, 400);
    assert!(first.samples.len() > 20);
    let drift = (first.final_energy - first.initial_energy).abs();
    assert!(
        drift < 1e-3,
        "NVE drift {drift} eV from {}",
        first.initial_energy
    );
    assert!(first.initial_energy > 0.0);

    let second = molecular_dynamics(&calculator, &structure, &nve_options()).unwrap();
    assert_eq!(first.final_energy, second.final_energy);
    for (left, right) in first.structure.sites().iter().zip(second.structure.sites()) {
        assert_eq!(left.frac(), right.frac());
    }

    let other_seed = molecular_dynamics(
        &calculator,
        &structure,
        &MdOptions {
            seed: 43,
            ..nve_options()
        },
    )
    .unwrap();
    assert_ne!(first.final_energy, other_seed.final_energy);
}

#[test]
fn langevin_and_nose_hoover_hold_the_target_temperature() {
    let calculator = UffCalculator::with_defaults().unwrap();
    let structure = methane();
    let langevin = molecular_dynamics(
        &calculator,
        &structure,
        &MdOptions {
            ensemble: Ensemble::NvtLangevin {
                temperature: 300.0,
                friction: 0.01,
            },
            time_step_fs: 0.5,
            steps: 600,
            seed: 7,
            ..MdOptions::default()
        },
    )
    .unwrap();
    assert!(!langevin.cancelled);
    assert!(
        (250.0..350.0).contains(&langevin.mean_temperature),
        "Langevin mean temperature {}",
        langevin.mean_temperature
    );
    assert!(
        langevin
            .samples
            .iter()
            .all(|sample| sample.temperature.is_finite() && sample.temperature >= 0.0)
    );

    let nose_hoover = molecular_dynamics(
        &calculator,
        &structure,
        &MdOptions {
            ensemble: Ensemble::NvtNoseHoover {
                temperature: 300.0,
                relaxation_fs: 50.0,
            },
            time_step_fs: 0.5,
            steps: 600,
            seed: 7,
            ..MdOptions::default()
        },
    )
    .unwrap();
    assert!(!nose_hoover.cancelled);
    assert!(
        (250.0..350.0).contains(&nose_hoover.mean_temperature),
        "Nose-Hoover mean temperature {}",
        nose_hoover.mean_temperature
    );
    assert_eq!(langevin.ensemble.label(), "nvt-langevin");
    assert_eq!(nose_hoover.ensemble.label(), "nvt-nose-hoover");
}

#[test]
fn npt_moves_the_volume_with_the_external_pressure() {
    let calculator = UffCalculator::with_defaults().unwrap();
    let structure = argon_fcc(5.0);
    let run = |pressure| {
        molecular_dynamics(
            &calculator,
            &structure,
            &MdOptions {
                ensemble: Ensemble::Npt {
                    temperature: 300.0,
                    pressure_gpa: pressure,
                    friction: 0.02,
                    barostat_interval: 5,
                    volume_move: 0.02,
                },
                time_step_fs: 0.5,
                steps: 600,
                seed: 3,
                initial_temperature: 300.0,
                sample_every: 5,
                ..MdOptions::default()
            },
        )
        .unwrap()
    };
    let ambient = run(0.0);
    let compressed = run(5.0);
    assert!(!ambient.cancelled && !compressed.cancelled);
    let ambient_volume = mean_volume_tail(&ambient);
    let compressed_volume = mean_volume_tail(&compressed);
    assert!(
        compressed_volume < ambient_volume,
        "5 GPa gave {compressed_volume} A^3 and 0 GPa gave {ambient_volume} A^3"
    );
    assert!(
        compressed.mean_pressure_gpa > ambient.mean_pressure_gpa,
        "mean pressures {} and {}",
        compressed.mean_pressure_gpa,
        ambient.mean_pressure_gpa
    );
    assert!(ambient.samples.iter().any(|sample| sample.volume > 0.0));
    assert!(ambient.ensemble.pressure_gpa() == Some(0.0));
}

#[test]
fn dynamics_can_be_cancelled_and_reports_progress() {
    let calculator = UffCalculator::with_defaults().unwrap();
    let structure = methane();
    let token = CancelToken::new();
    let stop_after = 20_usize;
    let mut seen = 0_usize;
    let result = {
        let mut callback = |event: &ProgressEvent| {
            if let ProgressEvent::Dynamics { step, .. } = event {
                seen = seen.max(*step);
                if *step >= stop_after {
                    token.cancel();
                }
            }
        };
        let mut control = RunControl::new()
            .with_token(&token)
            .with_progress(&mut callback);
        molecular_dynamics_with(
            &calculator,
            &structure,
            &MdOptions {
                steps: 500,
                sample_every: 1,
                ..nve_options()
            },
            &mut control,
        )
        .unwrap()
    };
    assert!(result.cancelled);
    assert!(result.steps < 500);
    assert!(result.steps >= stop_after);
    assert!(seen >= stop_after);
    assert!(!result.samples.is_empty());
}

#[test]
fn invalid_dynamics_inputs_are_refused() {
    let calculator = UffCalculator::with_defaults().unwrap();
    let structure = methane();
    assert!(
        molecular_dynamics(
            &calculator,
            &structure,
            &MdOptions {
                time_step_fs: 0.0,
                ..nve_options()
            }
        )
        .is_err()
    );
    assert!(
        molecular_dynamics(
            &calculator,
            &structure,
            &MdOptions {
                steps: 0,
                ..nve_options()
            }
        )
        .is_err()
    );
    assert!(
        molecular_dynamics(
            &calculator,
            &structure,
            &MdOptions {
                ensemble: Ensemble::NvtLangevin {
                    temperature: -1.0,
                    friction: 0.01,
                },
                ..nve_options()
            }
        )
        .is_err()
    );
    let lattice = Lattice::cubic(5.0).unwrap();
    let symmetric = StructureBuilder::new("symmetric", lattice)
        .hall("P 4/m m m")
        .site(Site::new("C1", "C", [0.0, 0.0, 0.0]).unwrap())
        .build()
        .unwrap();
    let error = molecular_dynamics(&calculator, &symmetric, &nve_options()).unwrap_err();
    assert_eq!(error.code(), "invalid_input");
    assert!(error.message().contains("P1"));
}

#[test]
fn dynamics_options_and_samples_round_trip_as_json() {
    let options = MdOptions {
        ensemble: Ensemble::Npt {
            temperature: 320.0,
            pressure_gpa: 1.5,
            friction: 0.02,
            barostat_interval: 4,
            volume_move: 0.03,
        },
        steps: 25,
        seed: 99,
        ..MdOptions::default()
    };
    let text = serde_json::to_string(&options).unwrap();
    let parsed: MdOptions = serde_json::from_str(&text).unwrap();
    assert_eq!(parsed, options);
    assert_eq!(options.ensemble.temperature(), Some(320.0));
    assert_eq!(options.ensemble.pressure_gpa(), Some(1.5));
    assert_eq!(Ensemble::Nve.label(), "nve");
}
