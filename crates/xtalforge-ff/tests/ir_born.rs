//! Integration tests for Born-effective-charge IR intensities (Task 19.4.1).
//!
//! The diatomic case has a closed-form intensity:
//! `I = C (Z_1 - Z_2)^2 / (m_1 + m_2)` for opposite unit displacements, so the
//! test checks the tensor path against that expression rather than a golden
//! number from this implementation.

#![allow(
    clippy::unwrap_used,
    clippy::expect_used,
    clippy::panic,
    clippy::indexing_slicing
)]

use xtalforge_core::Structure;
use xtalforge_core::elements::mass_of;
use xtalforge_ff::charges::Charges;
use xtalforge_ff::vibrations::E2_PER_AMU_TO_KM_PER_MOL;
use xtalforge_ff::{BornCharges, IrOptions, VibrationalModes, ir_spectrum_with_born, qeq_charges};

fn molecule(name: &str, atoms: &[(&str, [f64; 3])]) -> Structure {
    let points: Vec<(String, [f64; 3])> = atoms
        .iter()
        .map(|(element, point)| ((*element).to_owned(), *point))
        .collect();
    Structure::molecule(name, &points, 8.0).expect("molecule")
}

fn modes(frequencies: Vec<f64>, vectors: Vec<Vec<[f64; 3]>>) -> VibrationalModes {
    let zero_modes = frequencies
        .iter()
        .enumerate()
        .filter(|(_, frequency)| frequency.abs() <= 5.0)
        .map(|(index, _)| index)
        .collect();
    let imaginary = frequencies
        .iter()
        .enumerate()
        .filter(|(_, frequency)| **frequency < -5.0)
        .map(|(index, _)| index)
        .collect();
    VibrationalModes {
        frequencies,
        modes: vectors,
        zero_modes,
        imaginary,
    }
}

fn diatomic() -> (Structure, VibrationalModes) {
    let structure = molecule(
        "carbon-hydrogen",
        &[("C", [0.0, 0.0, 0.0]), ("H", [1.09, 0.0, 0.0])],
    );
    let mut vectors = vec![vec![[0.0; 3]; 2]; 5];
    vectors.push(vec![[1.0, 0.0, 0.0], [-1.0, 0.0, 0.0]]);
    let modes = modes(vec![0.0, 0.0, 0.0, 0.0, 0.0, 3000.0], vectors);
    (structure, modes)
}

fn uniform(value: f64) -> [[f64; 3]; 3] {
    [[value, 0.0, 0.0], [0.0, value, 0.0], [0.0, 0.0, value]]
}

#[test]
fn diatomic_born_intensity_matches_the_analytic_dipole_derivative() {
    let (structure, modes) = diatomic();
    let born = BornCharges::from_tensors(vec![uniform(2.0), uniform(-2.0)]);
    let spectrum = ir_spectrum_with_born(&structure, &modes, &born, &IrOptions::default()).unwrap();

    assert_eq!(spectrum.peaks.len(), 1);
    let mass_c = mass_of("C").unwrap();
    let mass_h = mass_of("H").unwrap();
    let expected = E2_PER_AMU_TO_KM_PER_MOL * 16.0 / (mass_c + mass_h);
    assert!(
        (spectrum.peaks[0].intensity - expected).abs() <= 1e-9 * expected,
        "{} != {expected}",
        spectrum.peaks[0].intensity
    );
    assert_eq!(spectrum.peaks[0].mode, 5);
    assert!(spectrum.max_intensity() > 0.0);
}

#[test]
fn doubling_the_born_tensors_quadruples_the_intensity() {
    let (structure, modes) = diatomic();
    let base = BornCharges::from_tensors(vec![uniform(2.0), uniform(-2.0)]);
    let doubled = BornCharges::from_tensors(vec![uniform(4.0), uniform(-4.0)]);
    let first = ir_spectrum_with_born(&structure, &modes, &base, &IrOptions::default()).unwrap();
    let second =
        ir_spectrum_with_born(&structure, &modes, &doubled, &IrOptions::default()).unwrap();
    assert!(
        (second.peaks[0].intensity - 4.0 * first.peaks[0].intensity).abs()
            <= 1e-9 * first.peaks[0].intensity
    );
}

#[test]
fn a_uniform_translation_has_no_intensity_for_neutral_born_tensors() {
    let structure = molecule(
        "carbon-hydrogen",
        &[("C", [0.0, 0.0, 0.0]), ("H", [1.09, 0.0, 0.0])],
    );
    let mut vectors = vec![vec![[0.0; 3]; 2]; 5];
    vectors.push(vec![[1.0, 0.0, 0.0], [1.0, 0.0, 0.0]]);
    let modes = modes(vec![0.0, 0.0, 0.0, 0.0, 0.0, 120.0], vectors);
    let born = BornCharges::from_tensors(vec![uniform(2.0), uniform(-2.0)]);
    let spectrum = ir_spectrum_with_born(&structure, &modes, &born, &IrOptions::default()).unwrap();
    assert_eq!(spectrum.peaks.len(), 1);
    assert!(spectrum.peaks[0].intensity.abs() < 1e-12);
}

#[test]
fn anisotropic_tensor_only_sees_the_displacement_direction() {
    let (structure, modes) = diatomic();
    // The mode moves along x, so the y and z diagonal elements are invisible.
    let carbon = [[2.0, 0.0, 0.0], [0.0, 9.0, 0.0], [0.0, 0.0, 9.0]];
    let hydrogen = [[-2.0, 0.0, 0.0], [0.0, -9.0, 0.0], [0.0, 0.0, -9.0]];
    let born = BornCharges::from_tensors(vec![carbon, hydrogen]);
    let anisotropic =
        ir_spectrum_with_born(&structure, &modes, &born, &IrOptions::default()).unwrap();
    let isotropic = ir_spectrum_with_born(
        &structure,
        &modes,
        &BornCharges::from_tensors(vec![uniform(2.0), uniform(-2.0)]),
        &IrOptions::default(),
    )
    .unwrap();
    assert!(
        (anisotropic.peaks[0].intensity - isotropic.peaks[0].intensity).abs()
            <= 1e-9 * isotropic.peaks[0].intensity
    );
}

#[test]
fn born_charges_from_scalar_charges_are_diagonal() {
    let charges = Charges {
        per_atom: vec![0.5, -0.25, -0.25],
        total: 0.0,
        iterations: 1,
    };
    let born = BornCharges::from_charges(&charges);
    assert_eq!(born.len(), 3);
    for (tensor, charge) in born.per_atom.iter().zip(&charges.per_atom) {
        for (row, values) in tensor.iter().enumerate() {
            for (column, value) in values.iter().enumerate() {
                let expected = if row == column { *charge } else { 0.0 };
                assert!((value - expected).abs() < 1e-15);
            }
        }
    }
}

#[test]
fn qeq_and_eqeq_construct_diagonal_born_tensors_of_the_right_size() {
    let water = molecule(
        "water",
        &[
            ("O", [0.0, 0.0, 0.0]),
            ("H", [0.96, 0.0, 0.0]),
            ("H", [-0.24, 0.93, 0.0]),
        ],
    );
    let qeq = BornCharges::qeq(&water).unwrap();
    let eqeq = BornCharges::eqeq(&water).unwrap();
    assert_eq!(qeq.len(), 3);
    assert_eq!(eqeq.len(), 3);
    let charges = qeq_charges(&water, &xtalforge_ff::ChargeOptions::default()).unwrap();
    for (tensor, charge) in qeq.per_atom.iter().zip(&charges.per_atom) {
        assert!((tensor[0][0] - charge).abs() < 1e-15);
        assert!((tensor[1][1] - charge).abs() < 1e-15);
        assert!((tensor[2][2] - charge).abs() < 1e-15);
    }
    for tensor in &eqeq.per_atom {
        assert!(tensor.iter().flatten().all(|value| value.is_finite()));
    }
    assert!(qeq.validate(3).is_ok());
}

#[test]
fn born_tensor_text_file_matches_explicit_tensors() {
    let text = "\
        # isotropic charges\n\
        1  2 0 0  0 2 0  0 0 2\n\
        2 -2 0 0  0 -2 0  0 0 -2\n";
    let parsed = BornCharges::from_text(text).unwrap();
    let explicit = BornCharges::from_tensors(vec![uniform(2.0), uniform(-2.0)]);
    assert_eq!(parsed, explicit);
    let (structure, modes) = diatomic();
    let spectrum =
        ir_spectrum_with_born(&structure, &modes, &parsed, &IrOptions::default()).unwrap();
    assert_eq!(spectrum.peaks.len(), 1);
}

#[test]
fn neutralised_born_charges_remove_the_mean_tensor() {
    let born = BornCharges::from_tensors(vec![uniform(3.0), uniform(-1.0)]);
    let neutral = born.neutralised();
    let total_x = neutral
        .per_atom
        .iter()
        .map(|tensor| tensor[0][0])
        .sum::<f64>();
    assert!(total_x.abs() < 1e-15);
    assert!((neutral.per_atom[0][0][0] - 2.0).abs() < 1e-15);
    assert!((neutral.per_atom[1][0][0] + 2.0).abs() < 1e-15);
    let empty = BornCharges::from_tensors(Vec::new());
    assert!(empty.is_empty());
    assert_eq!(empty.neutralised(), empty);
}

#[test]
fn invalid_born_inputs_are_rejected() {
    let (structure, modes) = diatomic();
    let options = IrOptions::default();

    let wrong_count = BornCharges::from_tensors(vec![uniform(2.0)]);
    let error = ir_spectrum_with_born(&structure, &modes, &wrong_count, &options).unwrap_err();
    assert_eq!(error.code(), "invalid_input");

    let mut not_finite = uniform(2.0);
    not_finite[1][2] = f64::NAN;
    let error = ir_spectrum_with_born(
        &structure,
        &modes,
        &BornCharges::from_tensors(vec![not_finite, uniform(-2.0)]),
        &options,
    )
    .unwrap_err();
    assert_eq!(error.code(), "invalid_input");

    let error = ir_spectrum_with_born(
        &structure,
        &modes,
        &BornCharges::from_tensors(vec![uniform(2.0), uniform(-2.0)]),
        &IrOptions {
            fwhm: 0.0,
            ..IrOptions::default()
        },
    )
    .unwrap_err();
    assert_eq!(error.code(), "invalid_input");

    let error = BornCharges::from_text("1 2 0 0 0 2 0 0").unwrap_err();
    assert_eq!(error.code(), "invalid_input");
}
