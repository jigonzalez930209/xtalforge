//! Integration tests for harmonic thermochemistry and the quasi-harmonic
//! free-energy evaluation (Task 4.8.4 and 4.8.5).
//!
//! The checks are the identities a reader can reproduce: the zero-point energy
//! is half the frequency sum, entropy grows with temperature, imaginary and
//! zero modes are excluded, pressure lowers translational entropy, and the
//! quasi-harmonic minimum lies at the expected volume.

#![allow(
    clippy::unwrap_used,
    clippy::expect_used,
    clippy::panic,
    clippy::indexing_slicing
)]

use xtalforge_core::Structure;
use xtalforge_ff::{
    CM_TO_EV, QhaPoint, ThermoOptions, VibrationalModes, ideal_gas_thermochemistry,
    quasi_harmonic_free_energy, thermochemistry,
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

fn modes_with(frequencies: Vec<f64>, n_atoms: usize) -> VibrationalModes {
    let n_modes = frequencies.len();
    VibrationalModes {
        frequencies,
        modes: vec![vec![[0.0; 3]; n_atoms]; n_modes],
        zero_modes: Vec::new(),
        imaginary: Vec::new(),
    }
}

#[test]
fn zero_point_energy_is_half_the_frequency_sum() {
    let frequencies = vec![1000.0, 2000.0, 3000.0];
    let result = thermochemistry(
        &frequencies,
        &ThermoOptions {
            temperature: 298.15,
            pressure: 1.0,
        },
    )
    .unwrap();
    let expected = 0.5 * frequencies.iter().sum::<f64>() * CM_TO_EV;
    assert!(
        (result.zpe - expected).abs() < 1.0e-12 * expected.max(1.0),
        "ZPE {} vs {expected}",
        result.zpe
    );
}

#[test]
fn entropy_increases_with_temperature() {
    let frequencies = vec![500.0, 1500.0, 3000.0];
    let mut previous = f64::NEG_INFINITY;
    for temperature in [50.0, 100.0, 298.15, 1000.0, 3000.0] {
        let result = thermochemistry(
            &frequencies,
            &ThermoOptions {
                temperature,
                pressure: 1.0,
            },
        )
        .unwrap();
        assert!(result.entropy.is_finite() && result.entropy > 0.0);
        assert!(
            result.entropy > previous,
            "entropy did not increase from the previous temperature at {temperature} K"
        );
        previous = result.entropy;
    }
}

#[test]
fn imaginary_and_zero_modes_are_excluded() {
    let frequencies = vec![-500.0, -1.0, 0.0, 4.0, 1000.0];
    let result = thermochemistry(
        &frequencies,
        &ThermoOptions {
            temperature: 300.0,
            pressure: 1.0,
        },
    )
    .unwrap();
    let expected = 0.5 * 1000.0 * CM_TO_EV;
    assert!((result.zpe - expected).abs() < 1.0e-12);
    assert!(result.thermal_u > 0.0);
    assert!(result.entropy > 0.0);
    assert!(result.gibbs.is_finite());
}

#[test]
fn invalid_thermo_options_are_rejected() {
    for temperature in [0.0, -1.0, f64::NAN, f64::INFINITY] {
        let error = thermochemistry(
            &[1000.0],
            &ThermoOptions {
                temperature,
                pressure: 1.0,
            },
        )
        .unwrap_err();
        assert_eq!(error.code(), "invalid_input");
    }
    for pressure in [0.0, -1.0, f64::NAN, f64::INFINITY] {
        let error = thermochemistry(
            &[1000.0],
            &ThermoOptions {
                temperature: 300.0,
                pressure,
            },
        )
        .unwrap_err();
        assert_eq!(error.code(), "invalid_input");
    }
    let error = thermochemistry(
        &[f64::NAN],
        &ThermoOptions {
            temperature: 300.0,
            pressure: 1.0,
        },
    )
    .unwrap_err();
    assert_eq!(error.code(), "invalid_input");
}

#[test]
fn ideal_gas_entropy_falls_with_pressure() {
    let structure = water();
    let modes = modes_with(vec![500.0, 1000.0, 2000.0], 3);
    let base = ThermoOptions {
        temperature: 500.0,
        pressure: 1.0,
    };
    let low = ideal_gas_thermochemistry(&structure, &modes, &base).unwrap();
    let high = ideal_gas_thermochemistry(
        &structure,
        &modes,
        &ThermoOptions {
            pressure: 10.0,
            ..base
        },
    )
    .unwrap();
    assert!(low.entropy > high.entropy);
    assert!(low.gibbs.is_finite() && high.gibbs.is_finite());
    assert!(low.thermal_u <= low.thermal_h);
}

#[test]
fn quasi_harmonic_minimum_is_refined_between_volume_points() {
    let points = vec![
        QhaPoint {
            volume: 10.0,
            static_energy: 0.1,
            frequencies: vec![1000.0, 1200.0, 1400.0],
        },
        QhaPoint {
            volume: 20.0,
            static_energy: -0.1,
            frequencies: vec![950.0, 1150.0, 1350.0],
        },
        QhaPoint {
            volume: 30.0,
            static_energy: 0.1,
            frequencies: vec![900.0, 1100.0, 1300.0],
        },
    ];
    let result = quasi_harmonic_free_energy(
        &points,
        &ThermoOptions {
            temperature: 300.0,
            pressure: 1.0,
        },
    )
    .unwrap();
    assert_eq!(result.points.len(), 3);
    // The interior sample is refined by a parabola through its neighbours, so
    // the minimum need not sit exactly on the sampled 20 Angstrom^3 point.
    assert!(
        (result.minimum_volume - 20.0).abs() < 1.0,
        "minimum volume {} is not near 20",
        result.minimum_volume
    );
    assert!(result.minimum_free_energy < result.points[0].free_energy);
    assert!(result.minimum_free_energy < result.points[2].free_energy);
}

#[test]
fn quasi_harmonic_rejects_non_increasing_volumes() {
    let points = vec![
        QhaPoint {
            volume: 20.0,
            static_energy: 0.0,
            frequencies: vec![1000.0],
        },
        QhaPoint {
            volume: 20.0,
            static_energy: 0.0,
            frequencies: vec![1000.0],
        },
    ];
    let error = quasi_harmonic_free_energy(
        &points,
        &ThermoOptions {
            temperature: 300.0,
            pressure: 1.0,
        },
    )
    .unwrap_err();
    assert_eq!(error.code(), "invalid_input");
}
