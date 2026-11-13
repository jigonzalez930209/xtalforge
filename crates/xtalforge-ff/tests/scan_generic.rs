//! Integration tests for generic scan coordinates, checkpoints, cancellation
//! and scan trajectories (Tasks 14.2.2 and 14.2.3).

#![allow(
    clippy::unwrap_used,
    clippy::expect_used,
    clippy::panic,
    clippy::indexing_slicing
)]

use std::path::PathBuf;

use xtalforge_core::{Lattice, Site, Structure, StructureBuilder};
use xtalforge_ff::{
    CancelToken, CellFreedom, CellParameter, CollectiveVariable, Constraint, Minimizer,
    OptimizeOptions, ProgressEvent, RelaxMode, RunControl, ScanAxis, ScanCheckpoint,
    ScanCoordinate, ScanGeometry, ScanOptions, TrajectoryOptions, UffCalculator, relaxed_scan,
    relaxed_scan_with,
};

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

fn water_dimer() -> Structure {
    let points: Vec<(String, [f64; 3])> = [
        ("O", [0.0, 0.0, 0.0]),
        ("H", [0.96, 0.0, 0.0]),
        ("H", [-0.24, 0.93, 0.0]),
        ("O", [4.0, 0.0, 0.0]),
        ("H", [4.96, 0.0, 0.0]),
        ("H", [3.76, 0.93, 0.0]),
    ]
    .iter()
    .map(|(element, point)| ((*element).to_owned(), *point))
    .collect();
    Structure::molecule("water-dimer", &points, 15.0).expect("dimer")
}

fn scan_options() -> ScanOptions {
    ScanOptions {
        optimizer: OptimizeOptions {
            max_steps: 200,
            force_tolerance: 0.02,
            energy_tolerance: 1e-8,
            mode: RelaxMode::P1,
            ..OptimizeOptions::default()
        },
        penalty: 500.0,
        minimizer: Minimizer::Lbfgs,
        ..ScanOptions::default()
    }
}

fn temporary(name: &str) -> PathBuf {
    std::env::temp_dir().join(format!("xtalforge-scan-{}-{name}", std::process::id()))
}

#[test]
fn cell_parameter_scan_imposes_the_parameter_and_finds_a_minimum() {
    let calculator = UffCalculator::with_defaults().unwrap();
    let structure = argon_fcc(5.0);
    let coordinate = ScanCoordinate::cell_parameter(CellParameter::A);
    assert_eq!(coordinate.label(), "cell a");
    assert_eq!(coordinate.units(), "A");
    assert!(coordinate.is_cell_parameter());
    assert_eq!(coordinate.parameter(), Some(CellParameter::A));
    assert!(coordinate.indices().is_empty());
    assert_eq!(CellParameter::Alpha.units(), "deg");

    let axis = ScanAxis::new(coordinate.clone(), vec![4.8, 5.4, 6.0, 6.6]).unwrap();
    let result = relaxed_scan(&calculator, &structure, &axis, &scan_options()).unwrap();
    assert_eq!(result.points.len(), 4);
    for point in &result.points {
        assert!(
            (point.achieved - point.target).abs() < 1e-9,
            "target {} reached {}",
            point.target,
            point.achieved
        );
        assert!((point.structure.lattice().parameters()[0] - point.target).abs() < 1e-9);
    }
    // Imposing a cell parameter keeps every fractional coordinate and scales
    // the Cartesian components along the changed cell vector.
    let mut seeded_geometry = ScanGeometry::from_structure(&structure);
    coordinate.set_in(&mut seeded_geometry, 4.8).unwrap();
    let seeded = seeded_geometry.into_structure(&structure).unwrap();
    assert!((seeded.lattice().parameters()[0] - 4.8).abs() < 1e-12);
    assert!((seeded.sites()[2].frac()[0] - 0.5).abs() < 1e-12);
    assert!((seeded.cart_sites()[2][0] - 2.4).abs() < 1e-12);
    assert!((seeded.cart_sites()[1][1] - 2.5).abs() < 1e-12);
    let energies = result.energies();
    eprintln!("cell scan targets {:?} energies {energies:?}", axis.values);
    assert!(
        energies[1] < energies[0] || energies[2] < energies[0],
        "energies {energies:?}"
    );
    let minimum = result.minimum().unwrap();
    assert!(minimum.target > 4.8 && minimum.target < 6.6);

    let mut positions = structure.cart_sites();
    let lattice = *structure.lattice();
    assert!(coordinate.set(&mut positions, &lattice, 5.4).is_err());
    let mut geometry = ScanGeometry::from_structure(&structure);
    coordinate.set_in(&mut geometry, 5.4).unwrap();
    assert!((geometry.lattice.parameters()[0] - 5.4).abs() < 1e-12);
    assert_eq!(
        coordinate.value(&positions, &geometry.lattice).unwrap(),
        5.4
    );
    assert!(coordinate.gradient(&positions, &lattice).unwrap()[0][0] == 0.0);
}

#[test]
fn cell_parameter_scan_holds_the_axis_under_cell_relaxation() {
    let calculator = UffCalculator::with_defaults().unwrap();
    let structure = argon_fcc(5.0);
    let options = ScanOptions {
        optimizer: OptimizeOptions {
            cell: true,
            max_steps: 120,
            force_tolerance: 0.05,
            ..scan_options().optimizer
        },
        ..scan_options()
    };
    let axis = ScanAxis::new(
        ScanCoordinate::cell_parameter(CellParameter::A),
        vec![5.0, 5.2],
    )
    .unwrap();
    let result = relaxed_scan(&calculator, &structure, &axis, &options).unwrap();
    let first = &result.points[0];
    assert!((first.achieved - 5.0).abs() < 1e-9);
    let parameters = first.structure.lattice().parameters();
    assert!((parameters[0] - 5.0).abs() < 1e-9);
    assert!(
        (first.structure.lattice().volume().abs() - structure.lattice().volume().abs()).abs()
            > 1e-6,
        "the free lengths and angles should have moved"
    );
}

#[test]
fn cell_axis_scan_can_hold_the_volume() {
    let calculator = UffCalculator::with_defaults().unwrap();
    let structure = argon_fcc(5.0);
    let options = ScanOptions {
        optimizer: OptimizeOptions {
            cell: true,
            cell_freedom: CellFreedom::fixed_volume(),
            max_steps: 80,
            force_tolerance: 0.1,
            ..scan_options().optimizer
        },
        ..scan_options()
    };
    let axis = ScanAxis::new(
        ScanCoordinate::cell_parameter(CellParameter::A),
        vec![5.0, 5.1],
    )
    .unwrap();
    let result = relaxed_scan(&calculator, &structure, &axis, &options).unwrap();
    let initial_volume = structure.lattice().volume().abs();
    for point in &result.points {
        assert!((point.achieved - point.target).abs() < 1e-9);
        assert!(
            (point.structure.lattice().volume().abs() - initial_volume).abs() < 1e-6,
            "volume {} left {initial_volume}",
            point.structure.lattice().volume()
        );
    }
}

#[test]
fn cell_angle_scan_imposes_gamma() {
    let calculator = UffCalculator::with_defaults().unwrap();
    let structure = argon_fcc(5.0);
    let options = ScanOptions {
        optimizer: OptimizeOptions {
            cell: true,
            max_steps: 60,
            force_tolerance: 0.2,
            ..scan_options().optimizer
        },
        ..scan_options()
    };
    let axis = ScanAxis::new(
        ScanCoordinate::cell_parameter(CellParameter::Gamma),
        vec![92.0, 96.0],
    )
    .unwrap();
    let result = relaxed_scan(&calculator, &structure, &axis, &options).unwrap();
    for point in &result.points {
        assert!(
            (point.achieved - point.target).abs() < 1e-9,
            "gamma target {} reached {}",
            point.target,
            point.achieved
        );
    }
}

#[test]
fn collective_variables_have_values_gradients_and_seeds() {
    let structure = water_dimer();
    let positions = structure.cart_sites();
    let lattice = *structure.lattice();

    let distance = CollectiveVariable::centroid_distance(vec![0, 1, 2], vec![3, 4, 5]).unwrap();
    let value = distance.value(&positions, &lattice).unwrap();
    assert!((value - 4.0).abs() < 1e-9, "centroid distance {value}");
    let gradient = distance.gradient(&positions, &lattice).unwrap();
    assert_eq!(gradient.len(), 6);
    let total: [f64; 3] = gradient.iter().fold([0.0; 3], |sum, value| {
        [sum[0] + value[0], sum[1] + value[1], sum[2] + value[2]]
    });
    assert!(total[0].abs() < 1e-9 && total[1].abs() < 1e-9 && total[2].abs() < 1e-9);
    let mut moved = positions.clone();
    distance.set(&mut moved, &lattice, 3.25).unwrap();
    eprintln!(
        "centroid after set {} lattice {:?}",
        distance.value(&moved, &lattice).unwrap(),
        structure.lattice().parameters()
    );
    assert!((distance.value(&moved, &lattice).unwrap() - 3.25).abs() < 1e-9);
    assert_eq!(distance.units(), "A");
    assert_eq!(distance.indices().len(), 6);
    assert!(CollectiveVariable::centroid_distance(Vec::new(), vec![1]).is_err());
    assert!(CollectiveVariable::centroid_distance(vec![0], vec![0]).is_err());

    let linear = CollectiveVariable::linear(vec![0, 1], vec![-1.0, 1.0], [1.0, 0.0, 0.0]).unwrap();
    let value = linear.value(&positions, &lattice).unwrap();
    assert!((value - 0.96).abs() < 1e-12, "linear value {value}");
    let mut moved = positions.clone();
    linear.set(&mut moved, &lattice, 1.5).unwrap();
    assert!((linear.value(&moved, &lattice).unwrap() - 1.5).abs() < 1e-12);
    assert!(CollectiveVariable::linear(vec![0], vec![1.0, 2.0], [1.0, 0.0, 0.0]).is_err());
    assert!(CollectiveVariable::linear(vec![0], vec![1.0], [0.0, 0.0, 0.0]).is_err());

    let coordinate = ScanCoordinate::collective(distance);
    assert_eq!(coordinate.units(), "A");
    assert_eq!(coordinate.indices().len(), 6);
    assert!(coordinate.validate(6).is_ok());
    assert!(coordinate.validate(2).is_err());
}

#[test]
fn collective_scan_is_monotone_and_reaches_its_targets() {
    let calculator = UffCalculator::with_defaults().unwrap();
    let structure = water_dimer();
    let variable = CollectiveVariable::centroid_distance(vec![0, 1, 2], vec![3, 4, 5]).unwrap();
    let axis = ScanAxis::linspace(ScanCoordinate::collective(variable), 3.4, 4.4, 3).unwrap();
    let result = relaxed_scan(&calculator, &structure, &axis, &scan_options()).unwrap();
    assert_eq!(result.points.len(), 3);
    let mut previous = f64::NEG_INFINITY;
    for point in &result.points {
        assert!(
            (point.achieved - point.target).abs() < 0.05,
            "target {} reached {}",
            point.target,
            point.achieved
        );
        assert!(point.achieved > previous);
        previous = point.achieved;
    }
}

#[test]
fn checkpoints_resume_at_the_first_unfinished_point() {
    let calculator = UffCalculator::with_defaults().unwrap();
    let structure = argon_fcc(5.0);
    let path = temporary("resume.json");
    let _ = std::fs::remove_file(&path);
    let axis = ScanAxis::new(
        ScanCoordinate::cell_parameter(CellParameter::A),
        vec![4.7, 4.9, 5.1, 5.3],
    )
    .unwrap();
    let options = ScanOptions {
        checkpoint: Some(path.clone()),
        ..scan_options()
    };
    let full = relaxed_scan(&calculator, &structure, &axis, &options).unwrap();
    assert_eq!(full.points.len(), 4);
    assert_eq!(full.resumed_from, 0);

    let _ = std::fs::remove_file(&path);
    let token = CancelToken::new();
    let interrupted = {
        let mut callback = |event: &ProgressEvent| {
            if let ProgressEvent::Scan { point, .. } = event
                && *point == 1
            {
                token.cancel();
            }
        };
        let mut control = RunControl::new()
            .with_token(&token)
            .with_progress(&mut callback);
        relaxed_scan_with(&calculator, &structure, &axis, &options, &mut control).unwrap()
    };
    assert!(interrupted.cancelled);
    assert_eq!(interrupted.points.len(), 2);
    let checkpoint = ScanCheckpoint::read(&path).unwrap();
    assert_eq!(checkpoint.points.len(), 2);
    assert_eq!(checkpoint.resumed_len(), 2);
    assert!(!checkpoint.is_complete());
    checkpoint
        .matches(&axis.coordinate, &axis.values, 4)
        .unwrap();

    let resumed = relaxed_scan(&calculator, &structure, &axis, &options).unwrap();
    assert_eq!(resumed.resumed_from, 2);
    assert_eq!(resumed.points.len(), 4);
    assert!(!resumed.cancelled);
    for (resumed, reference) in resumed.points.iter().zip(&full.points) {
        assert!((resumed.energy - reference.energy).abs() < 1e-6);
        assert!((resumed.achieved - reference.achieved).abs() < 1e-9);
    }

    let jsonl_path = temporary("resume.jsonl");
    let _ = std::fs::remove_file(&jsonl_path);
    let jsonl_options = ScanOptions {
        checkpoint: Some(jsonl_path.clone()),
        ..scan_options()
    };
    let shorter = ScanAxis::new(
        ScanCoordinate::cell_parameter(CellParameter::A),
        vec![4.8, 5.0],
    )
    .unwrap();
    let _ = relaxed_scan(&calculator, &structure, &shorter, &jsonl_options).unwrap();
    let text = std::fs::read_to_string(&jsonl_path).unwrap();
    assert_eq!(text.lines().count(), 3);
    let parsed = ScanCheckpoint::from_jsonl(&text).unwrap();
    assert_eq!(parsed.points.len(), 2);
    assert!(parsed.is_complete());
    assert_eq!(
        ScanCheckpoint::from_json(&parsed.to_json().unwrap())
            .unwrap()
            .points
            .len(),
        2
    );

    let mismatched = ScanAxis::new(
        ScanCoordinate::cell_parameter(CellParameter::A),
        vec![4.7, 4.9, 5.2, 5.3],
    )
    .unwrap();
    let error = relaxed_scan(&calculator, &structure, &mismatched, &options).unwrap_err();
    assert_eq!(error.code(), "invalid_input");
    let _ = std::fs::remove_file(path);
    let _ = std::fs::remove_file(jsonl_path);
}

#[test]
fn scan_trajectories_and_invalid_coordinates() {
    let calculator = UffCalculator::with_defaults().unwrap();
    let structure = argon_fcc(5.0);
    let options = ScanOptions {
        trajectory: TrajectoryOptions::full(),
        ..scan_options()
    };
    let axis = ScanAxis::new(
        ScanCoordinate::cell_parameter(CellParameter::A),
        vec![4.9, 5.1],
    )
    .unwrap();
    let result = relaxed_scan(&calculator, &structure, &axis, &options).unwrap();
    assert!(!result.trajectory.is_empty());
    assert_eq!(result.trajectory.frames[0].point, Some(0));
    assert_eq!(
        result.trajectory.frames[0].positions.as_ref().map(Vec::len),
        Some(structure.sites().len())
    );

    assert!(Constraint::new(ScanCoordinate::cell_parameter(CellParameter::A), 5.0).is_err());
    let cell_axis = ScanAxis::new(
        ScanCoordinate::cell_parameter(CellParameter::A),
        vec![5.0, 5.1],
    )
    .unwrap();
    assert!((cell_axis.label()).contains("cell a"));
    assert!(ScanAxis::new(ScanCoordinate::Distance { i: 0, j: 1 }, vec![]).is_err());
    let bad_checkpoint = ScanOptions {
        checkpoint: Some(PathBuf::from("/")),
        ..ScanOptions::default()
    };
    assert!(bad_checkpoint.validate().is_err());
    let bad_collective = ScanCoordinate::Collective {
        variable: CollectiveVariable::CentroidDistance {
            first: vec![0],
            second: vec![99],
        },
    };
    assert!(bad_collective.validate(4).is_err());
    let axis = ScanAxis::new(bad_collective, vec![1.0, 2.0]).unwrap();
    assert!(relaxed_scan(&calculator, &structure, &axis, &scan_options()).is_err());
}
