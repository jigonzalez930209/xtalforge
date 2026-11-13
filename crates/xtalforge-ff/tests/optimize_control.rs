//! Integration tests for cell freedom, external pressure, cancellation,
//! progress reporting and trajectories (Tasks 14.2.1 and 14.2.2).

#![allow(
    clippy::unwrap_used,
    clippy::expect_used,
    clippy::panic,
    clippy::indexing_slicing
)]

use std::sync::Arc;
use std::sync::atomic::AtomicBool;

use xtalforge_core::{Lattice, Site, Structure, StructureBuilder};
use xtalforge_ff::{
    CancelToken, CellFreedom, Minimizer, OptimizeOptions, ProgressEvent, RelaxMode, RunControl,
    Trajectory, TrajectoryOptions, UffCalculator, fire, fire_with, lbfgs, lbfgs_with, minimize,
    minimize_with,
};

fn molecule(name: &str, atoms: &[(&str, [f64; 3])]) -> Structure {
    let points: Vec<(String, [f64; 3])> = atoms
        .iter()
        .map(|(element, point)| ((*element).to_owned(), *point))
        .collect();
    Structure::molecule(name, &points, 12.0).expect("molecule")
}

fn distorted_methane() -> Structure {
    let length = 1.15;
    let a = length / 3.0_f64.sqrt();
    molecule(
        "distorted-methane",
        &[
            ("C", [0.0, 0.0, 0.0]),
            ("H", [a, a, a]),
            ("H", [a, -a, -a]),
            ("H", [-a, a, -a]),
            ("H", [-a, -a, a]),
        ],
    )
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

fn cell_options() -> OptimizeOptions {
    OptimizeOptions {
        max_steps: 300,
        force_tolerance: 0.02,
        energy_tolerance: 1e-7,
        cell: true,
        mode: RelaxMode::P1,
        max_step: 0.2,
        ..OptimizeOptions::default()
    }
}

#[test]
fn held_angles_and_lengths_stay_exact() {
    let calculator = UffCalculator::with_defaults().unwrap();
    let structure = argon_fcc(5.5);
    let initial = structure.lattice().parameters();

    let lengths = OptimizeOptions {
        cell_freedom: CellFreedom::lengths_only(),
        ..cell_options()
    };
    let result = lbfgs(&calculator, &structure, &lengths).unwrap();
    let parameters = result.structure.lattice().parameters();
    for index in 0..3 {
        assert!(
            (parameters[3 + index] - initial[3 + index]).abs() < 1e-9,
            "held angle {} moved from {} to {}",
            index,
            initial[3 + index],
            parameters[3 + index]
        );
    }
    assert!(parameters[0] > 0.0 && parameters[1] > 0.0 && parameters[2] > 0.0);

    let fixed_volume = OptimizeOptions {
        cell_freedom: CellFreedom::fixed_volume(),
        ..cell_options()
    };
    let held = lbfgs(&calculator, &structure, &fixed_volume).unwrap();
    let volume = held.structure.lattice().volume().abs();
    assert!(
        (volume - structure.lattice().volume().abs()).abs() < 1e-6,
        "held volume moved from {} to {volume}",
        structure.lattice().volume()
    );
    assert!(held.energy >= result.energy - 1e-6);

    // Held angles and a held volume can coexist: the volume rescale is
    // absorbed by the free lengths only.
    let shape = OptimizeOptions {
        cell_freedom: CellFreedom::lengths_only().hold_volume(),
        ..cell_options()
    };
    let shaped = lbfgs(&calculator, &structure, &shape).unwrap();
    let parameters = shaped.structure.lattice().parameters();
    for index in 0..3 {
        assert!((parameters[3 + index] - initial[3 + index]).abs() < 1e-9);
    }
    assert!(
        (shaped.structure.lattice().volume().abs() - structure.lattice().volume().abs()).abs()
            < 1e-6
    );
}

#[test]
fn external_pressure_contracts_the_cell() {
    let calculator = UffCalculator::with_defaults().unwrap();
    let structure = argon_fcc(5.5);
    let ambient = lbfgs(&calculator, &structure, &cell_options()).unwrap();
    let compressed = lbfgs(
        &calculator,
        &structure,
        &OptimizeOptions {
            pressure_gpa: 5.0,
            ..cell_options()
        },
    )
    .unwrap();
    let ambient_volume = ambient.structure.lattice().volume().abs();
    let compressed_volume = compressed.structure.lattice().volume().abs();
    assert!(
        compressed_volume < ambient_volume,
        "5 GPa gave {compressed_volume} and 0 GPa gave {ambient_volume}"
    );
    assert!((ambient.enthalpy - ambient.energy).abs() < 1e-12);
    assert!(compressed.enthalpy > compressed.energy);
}

#[test]
fn cancelled_runs_return_the_partial_geometry() {
    let calculator = UffCalculator::with_defaults().unwrap();
    let structure = distorted_methane();
    let options = OptimizeOptions {
        max_steps: 200,
        force_tolerance: 0.001,
        mode: RelaxMode::P1,
        ..OptimizeOptions::default()
    };
    for minimizer in [Minimizer::Fire, Minimizer::Lbfgs] {
        let token = CancelToken::new();
        token.cancel();
        let mut control = RunControl::new().with_token(&token);
        let result =
            minimize_with(&calculator, &structure, &options, minimizer, &mut control).unwrap();
        assert!(result.cancelled);
        assert!(!result.converged);
        assert_eq!(result.steps, 0);
        assert!(result.trajectory.is_empty());
        token.reset();
        assert!(!token.is_cancelled());
        assert!(fire(&calculator, &structure, &options).is_ok());
        assert!(lbfgs(&calculator, &structure, &options).is_ok());
    }
    let flag = Arc::new(AtomicBool::new(true));
    let result = fire_with(
        &calculator,
        &structure,
        &options,
        &mut RunControl::new().with_cancel(flag.as_ref()),
    )
    .unwrap();
    assert!(result.cancelled);

    // A run cancelled in the middle keeps every finished step.
    let slow = OptimizeOptions {
        max_step: 0.005,
        force_tolerance: 1e-4,
        ..options.clone()
    };
    let token = CancelToken::new();
    let mid_run = {
        let mut callback = |event: &ProgressEvent| {
            if let ProgressEvent::Optimizer { step, .. } = event
                && *step >= 2
            {
                token.cancel();
            }
        };
        let mut control = RunControl::new()
            .with_token(&token)
            .with_progress(&mut callback);
        lbfgs_with(&calculator, &structure, &slow, &mut control).unwrap()
    };
    assert!(mid_run.cancelled);
    assert!(mid_run.steps >= 2 && mid_run.steps < slow.max_steps);
    assert!(mid_run.energy.is_finite());
}

#[test]
fn progress_and_trajectory_cover_the_run() {
    let calculator = UffCalculator::with_defaults().unwrap();
    let structure = distorted_methane();
    let options = OptimizeOptions {
        max_steps: 100,
        force_tolerance: 0.002,
        mode: RelaxMode::P1,
        trajectory: TrajectoryOptions::full(),
        ..OptimizeOptions::default()
    };
    let mut events = Vec::new();
    let result = {
        let mut callback = |event: &ProgressEvent| events.push(event.clone());
        let mut control = RunControl::new().with_progress(&mut callback);
        minimize_with(
            &calculator,
            &structure,
            &options,
            Minimizer::Lbfgs,
            &mut control,
        )
        .unwrap()
    };
    assert!(!events.is_empty());
    let last = events.last().unwrap();
    match last {
        ProgressEvent::Optimizer { step, .. } => assert_eq!(*step, result.steps),
        other => panic!("unexpected event {other:?}"),
    }
    assert!(result.trajectory.len() >= 2);
    let first = &result.trajectory.frames[0];
    assert_eq!(first.step, 0);
    assert_eq!(
        first.positions.as_ref().map(Vec::len),
        Some(structure.sites().len())
    );
    assert!(first.volume > 0.0);
    let last_frame = result.trajectory.frames.last().unwrap();
    assert_eq!(last_frame.step, result.steps);

    let json = result.trajectory.to_json().unwrap();
    let parsed = Trajectory::from_json(&json).unwrap();
    // serde_json's default float parser is not bit-exact for every decimal
    // spelling, so compare what a caller reads rather than raw bits.
    assert_eq!(parsed.len(), result.trajectory.len());
    for (first, second) in parsed.frames.iter().zip(&result.trajectory.frames) {
        assert_eq!(first.step, second.step);
        assert!((first.energy - second.energy).abs() < 1e-12);
        assert!((first.volume - second.volume).abs() < 1e-9);
    }
    let jsonl = result.trajectory.to_jsonl().unwrap();
    assert_eq!(Trajectory::from_jsonl(&jsonl).unwrap().len(), parsed.len());

    let base = std::env::temp_dir().join(format!("xtalforge-trajectory-{}", std::process::id()));
    let json_path = base.with_extension("json");
    let jsonl_path = base.with_extension("jsonl");
    result.trajectory.write(&json_path).unwrap();
    result.trajectory.write(&jsonl_path).unwrap();
    assert_eq!(
        Trajectory::read(&json_path).unwrap().len(),
        result.trajectory.len()
    );
    assert_eq!(
        Trajectory::read(&jsonl_path).unwrap().len(),
        result.trajectory.len()
    );
    let _ = std::fs::remove_file(json_path);
    let _ = std::fs::remove_file(jsonl_path);
}

#[test]
fn invalid_cell_options_are_refused() {
    let calculator = UffCalculator::with_defaults().unwrap();
    let structure = distorted_methane();
    let invalid_pressure = OptimizeOptions {
        pressure_gpa: f64::NAN,
        ..OptimizeOptions::default()
    };
    assert!(minimize(&calculator, &structure, &invalid_pressure, Minimizer::Lbfgs).is_err());

    let held_cell = OptimizeOptions {
        cell: true,
        cell_freedom: CellFreedom::none(),
        fixed_atoms: vec![0, 1, 2, 3, 4],
        ..OptimizeOptions::default()
    };
    let error = minimize(&calculator, &structure, &held_cell, Minimizer::Lbfgs).unwrap_err();
    assert_eq!(error.code(), "invalid_input");
    assert!(error.message().contains("held"));

    let freedom = CellFreedom::fixed_volume().hold_length(0).hold_angle(2);
    assert!(!freedom.lengths[0]);
    assert!(freedom.lengths[1]);
    assert!(!freedom.angles[2]);
    assert!(!freedom.volume);
    assert!(!freedom.is_full());
    assert!(CellFreedom::full().is_full());
}

#[test]
fn older_options_json_still_loads() {
    // A record written before Task 14.2 has none of the new fields.
    let text = r#"{"max_steps":10,"force_tolerance":0.1,"energy_tolerance":1e-6,
        "cell":false,"fixed_atoms":[],"mode":"Asymmetric","max_step":0.2}"#;
    let options: OptimizeOptions = serde_json::from_str(text).unwrap();
    assert_eq!(options.cell_freedom, CellFreedom::full());
    assert_eq!(options.pressure_gpa, 0.0);
    assert!(!options.trajectory.enabled());

    let round_trip = serde_json::to_string(&options).unwrap();
    assert_eq!(
        serde_json::from_str::<OptimizeOptions>(&round_trip).unwrap(),
        options
    );

    let scan: xtalforge_ff::ScanOptions = serde_json::from_str(
        r#"{"optimizer": {"max_steps":10,"force_tolerance":0.1,"energy_tolerance":1e-6,
        "cell":false,"fixed_atoms":[],"mode":"Asymmetric","max_step":0.2},
        "penalty": 50.0, "minimizer": "Lbfgs"}"#,
    )
    .unwrap();
    assert!(scan.checkpoint.is_none());
    assert!(!scan.trajectory.enabled());
}
