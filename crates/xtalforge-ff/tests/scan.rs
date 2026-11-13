//! Integration tests for scan coordinates, penalties and the relaxed scan
//! driver (Task 4.5).

#![allow(
    clippy::unwrap_used,
    clippy::expect_used,
    clippy::panic,
    clippy::indexing_slicing
)]

use xtalforge_core::Structure;
use xtalforge_ff::{
    Constraint, Minimizer, OptimizeOptions, RelaxMode, ScanAxis, ScanCoordinate, ScanOptions,
    UffCalculator, relaxed_scan,
};

fn molecule(name: &str, atoms: &[(&str, [f64; 3])]) -> Structure {
    let points: Vec<(String, [f64; 3])> = atoms
        .iter()
        .map(|(element, point)| ((*element).to_owned(), *point))
        .collect();
    Structure::molecule(name, &points, 8.0).expect("molecule")
}

fn water() -> Structure {
    molecule(
        "water",
        &[
            ("O", [0.0, 0.0, 0.0]),
            ("H", [0.96, 0.0, 0.0]),
            ("H", [-0.24, 0.93, 0.0]),
        ],
    )
}

fn argon_dimer() -> Structure {
    molecule(
        "argon-dimer",
        &[("Ar", [0.0, 0.0, 0.0]), ("Ar", [3.4, 0.0, 0.0])],
    )
}

fn options() -> ScanOptions {
    ScanOptions {
        optimizer: OptimizeOptions {
            max_steps: 200,
            force_tolerance: 0.02,
            energy_tolerance: 1e-8,
            cell: false,
            fixed_atoms: Vec::new(),
            mode: RelaxMode::P1,
            max_step: 0.1,
            ..OptimizeOptions::default()
        },
        penalty: 100.0,
        minimizer: Minimizer::Lbfgs,
        ..ScanOptions::default()
    }
}

#[test]
fn distance_scan_is_monotone_and_holds_its_targets() {
    let calculator = UffCalculator::with_defaults().unwrap();
    let coordinate = ScanCoordinate::distance(0, 1).unwrap();
    let axis = ScanAxis::linspace(coordinate, 3.4, 4.4, 3).unwrap();
    let result = relaxed_scan(&calculator, &argon_dimer(), &axis, &options()).unwrap();
    assert_eq!(result.points.len(), 3);
    let mut previous = f64::NEG_INFINITY;
    for point in &result.points {
        assert!(point.energy.is_finite());
        assert!(
            (point.achieved - point.target).abs() < 0.02,
            "target {} reached {}",
            point.target,
            point.achieved
        );
        assert!(point.achieved > previous, "scan went backwards");
        previous = point.achieved;
    }
    assert!(result.minimum().is_some());
}

#[test]
fn angle_scan_is_monotone() {
    let calculator = UffCalculator::with_defaults().unwrap();
    let coordinate = ScanCoordinate::angle(1, 0, 2).unwrap();
    let axis = ScanAxis::linspace(coordinate, 96.0, 116.0, 3).unwrap();
    let result = relaxed_scan(&calculator, &water(), &axis, &options()).unwrap();
    assert_eq!(result.points.len(), 3);
    let mut previous = f64::NEG_INFINITY;
    for point in &result.points {
        assert!(
            (point.achieved - point.target).abs() < 1.0,
            "target {} reached {}",
            point.target,
            point.achieved
        );
        assert!(point.achieved > previous);
        previous = point.achieved;
    }
}

#[test]
fn torsion_value_and_gradient_are_finite() {
    let structure = molecule(
        "chain",
        &[
            ("C", [0.0, 0.0, 0.0]),
            ("C", [1.5, 0.0, 0.0]),
            ("C", [2.0, 1.4, 0.0]),
            ("H", [2.5, 1.8, 0.9]),
        ],
    );
    let coordinate = ScanCoordinate::torsion(0, 1, 2, 3).unwrap();
    let positions = structure.cart_sites();
    let value = coordinate.value(&positions, structure.lattice()).unwrap();
    assert!(value.is_finite() && value.abs() <= 180.0);
    let gradient = coordinate
        .gradient(&positions, structure.lattice())
        .unwrap();
    assert_eq!(gradient.len(), 4);
    let total: [f64; 3] = gradient.iter().fold([0.0; 3], |sum, force| {
        [sum[0] + force[0], sum[1] + force[1], sum[2] + force[2]]
    });
    assert!(total[0].abs() < 1e-9 && total[1].abs() < 1e-9 && total[2].abs() < 1e-9);
}

#[test]
fn constraint_penalty_is_quadratic() {
    let structure = water();
    let positions = structure.cart_sites();
    let constraint = Constraint::angle(1, 0, 2, 104.5).unwrap();
    let (energy, gradient) = constraint
        .penalty(&positions, structure.lattice(), 2.0)
        .unwrap();
    let residual = constraint
        .residual(&positions, structure.lattice())
        .unwrap();
    assert!((energy - residual * residual).abs() < 1e-12);
    assert_eq!(gradient.len(), 3);
}

#[test]
fn invalid_indices_are_refused() {
    let coordinate = ScanCoordinate::Distance { i: 0, j: 99 };
    assert!(coordinate.validate(3).is_err());
    assert!(ScanCoordinate::distance(1, 1).is_err());
    assert!(ScanCoordinate::angle(0, 0, 1).is_err());
    assert!(ScanCoordinate::torsion(0, 1, 1, 2).is_err());

    let calculator = UffCalculator::with_defaults().unwrap();
    let axis = ScanAxis::new(coordinate, vec![2.0, 2.5]).unwrap();
    let error = relaxed_scan(&calculator, &water(), &axis, &options()).unwrap_err();
    assert_eq!(error.code(), "invalid_input");
}

#[test]
fn scan_set_moves_the_coordinate_to_its_target() {
    let structure = water();
    let mut positions = structure.cart_sites();
    let coordinate = ScanCoordinate::distance(1, 2).unwrap();
    coordinate
        .set(&mut positions, structure.lattice(), 2.0)
        .unwrap();
    let achieved = coordinate.value(&positions, structure.lattice()).unwrap();
    assert!((achieved - 2.0).abs() < 1e-9);

    let coordinate = ScanCoordinate::angle(1, 0, 2).unwrap();
    coordinate
        .set(&mut positions, structure.lattice(), 120.0)
        .unwrap();
    let achieved = coordinate.value(&positions, structure.lattice()).unwrap();
    assert!((achieved - 120.0).abs() < 1e-9);
}
