//! Integration tests for Gamma-point vibrations, IR spectra and Gamma phonons
//! (Task 4.8.1 to 4.8.3).
//!
//! The diatomic check compares the finite-difference frequency with the
//! analytic harmonic value from UFF's own bond force constant.  Methane checks
//! the 3N - 6 rule and the six projected rigid motions.  The IR test only
//! requires physical positivity and a non-empty broadened curve; the absolute
//! intensity is a documented fixed-charge approximation and has no reference
//! value in this repository.

#![allow(
    clippy::unwrap_used,
    clippy::expect_used,
    clippy::panic,
    clippy::indexing_slicing
)]

use xtalforge_core::elements::mass_of;
use xtalforge_core::{Lattice, Site, Structure, StructureBuilder};
use xtalforge_ff::uff::{bond_force_constant, natural_bond_length};
use xtalforge_ff::{
    Broadening, CM_PER_SQRT_EV_PER_A2_AMU, Calculator, CalculatorOptions, DisplacementMode, Error,
    HessianOptions, IrOptions, PhononOptions, UffCalculator, VibrationalModes, hessian,
    ir_spectrum, phonons, vibrational_modes,
};

const KCAL_PER_MOL_TO_EV: f64 = 0.043_364_104_241_800_934;

fn molecule(name: &str, atoms: &[(&str, [f64; 3])]) -> Structure {
    let points: Vec<(String, [f64; 3])> = atoms
        .iter()
        .map(|(element, point)| ((*element).to_owned(), *point))
        .collect();
    Structure::molecule(name, &points, 8.0).expect("molecule")
}

fn tetrahedral_methane() -> Structure {
    let length = natural_bond_length("C_3", "H_", 1.0).unwrap();
    let a = length / 3.0_f64.sqrt();
    molecule(
        "methane",
        &[
            ("C", [0.0, 0.0, 0.0]),
            ("H", [a, a, a]),
            ("H", [a, -a, -a]),
            ("H", [-a, a, -a]),
            ("H", [-a, -a, a]),
        ],
    )
}

fn carbon_dimer() -> (Structure, f64) {
    let length = natural_bond_length("C_3", "C_3", 1.0).unwrap();
    let structure = molecule(
        "carbon-dimer",
        &[("C", [0.0, 0.0, 0.0]), ("C", [length, 0.0, 0.0])],
    );
    (structure, length)
}

fn diamond() -> Structure {
    let lattice = Lattice::cubic(3.566_986).unwrap();
    let mut builder = StructureBuilder::new("diamond", lattice);
    for (index, point) in [
        [0.0, 0.0, 0.0],
        [0.0, 0.5, 0.5],
        [0.5, 0.0, 0.5],
        [0.5, 0.5, 0.0],
    ]
    .iter()
    .enumerate()
    {
        builder = builder.site(Site::new(format!("C{index}"), "C", *point).unwrap());
    }
    for (index, point) in [
        [0.25, 0.25, 0.25],
        [0.25, 0.75, 0.75],
        [0.75, 0.25, 0.75],
        [0.75, 0.75, 0.25],
    ]
    .iter()
    .enumerate()
    {
        builder = builder.site(Site::new(format!("D{index}"), "C", *point).unwrap());
    }
    builder.build().unwrap()
}

#[test]
fn diatomic_frequency_matches_the_analytic_harmonic_value() {
    let (structure, length) = carbon_dimer();
    let options = CalculatorOptions {
        vdw: false,
        ..CalculatorOptions::default()
    };
    let calculator = UffCalculator::new(options).unwrap();
    let hessian = hessian(&calculator, &structure, &HessianOptions::molecule()).unwrap();
    let modes = vibrational_modes(&hessian, &structure, &HessianOptions::molecule()).unwrap();

    let stiffness = bond_force_constant("C_3", "C_3", length).unwrap();
    let reduced_mass = mass_of("C").unwrap() / 2.0;
    let expected =
        (stiffness * KCAL_PER_MOL_TO_EV / reduced_mass).sqrt() * CM_PER_SQRT_EV_PER_A2_AMU;

    assert_eq!(
        modes.zero_modes.len(),
        5,
        "a linear dimer has 5 rigid modes"
    );
    assert_eq!(
        modes.zero_modes.len(),
        5,
        "a linear dimer has 5 rigid modes"
    );
    assert!(modes.imaginary.is_empty(), "unexpected imaginary mode");
    let real: Vec<(usize, f64)> = modes.real_modes().collect();
    assert_eq!(real.len(), 1, "a diatomic has one vibration");
    let (_, frequency) = real[0];
    assert!(
        (frequency - expected).abs() / expected < 2.0e-3,
        "finite-difference {frequency} cm^-1 vs analytic {expected} cm^-1"
    );
}

#[test]
fn methane_has_3n_minus_6_nonzero_modes_and_6_zero_modes() {
    let structure = tetrahedral_methane();
    let calculator = UffCalculator::with_defaults().unwrap();
    let options = HessianOptions::molecule();
    let hessian = hessian(&calculator, &structure, &options).unwrap();
    let modes = vibrational_modes(&hessian, &structure, &options).unwrap();

    assert_eq!(modes.frequencies.len(), 15, "methane has 3N = 15 modes");
    assert_eq!(modes.zero_modes.len(), 6, "3 translations and 3 rotations");
    assert!(modes.imaginary.is_empty(), "methane should be a minimum");
    assert_eq!(modes.real_modes().count(), 9, "3N - 6 = 9 vibrations");
    for (_, frequency) in modes.real_modes() {
        assert!(
            frequency > 100.0,
            "a non-rigid methane mode at {frequency} cm^-1 is suspiciously soft"
        );
    }
}

#[test]
fn ir_spectrum_has_positive_intensities_and_a_broadened_curve() {
    let structure = tetrahedral_methane();
    let calculator = UffCalculator::with_defaults().unwrap();
    let options = HessianOptions::molecule();
    let hessian = hessian(&calculator, &structure, &options).unwrap();
    let modes = vibrational_modes(&hessian, &structure, &options).unwrap();
    let spectrum = ir_spectrum(&structure, &modes, &IrOptions::default()).unwrap();

    assert_eq!(spectrum.peaks.len(), 9, "one stick per real mode");
    assert!(
        spectrum.peaks.iter().any(|peak| peak.intensity > 0.0),
        "methane's T2 modes must absorb"
    );
    assert!(!spectrum.points.is_empty(), "spectrum grid is empty");
    assert!(spectrum.max_intensity() > 0.0);
    for point in &spectrum.points {
        assert!(point[0].is_finite() && point[1].is_finite());
        assert!(point[1] >= 0.0);
    }
    for peak in &spectrum.peaks {
        assert!(peak.intensity.is_finite() && peak.intensity >= 0.0);
    }

    let gaussian = ir_spectrum(
        &structure,
        &modes,
        &IrOptions {
            broadening: Broadening::Gaussian,
            ..IrOptions::default()
        },
    )
    .unwrap();
    assert_eq!(gaussian.points.len(), spectrum.points.len());
    assert!(gaussian.max_intensity() > 0.0);
}

#[test]
fn invalid_ir_options_are_rejected() {
    let structure = carbon_dimer().0;
    let modes = VibrationalModes {
        frequencies: Vec::new(),
        modes: Vec::new(),
        zero_modes: Vec::new(),
        imaginary: Vec::new(),
    };
    for (fwhm, step) in [(0.0, 5.0), (-1.0, 5.0), (10.0, 0.0), (10.0, f64::NAN)] {
        let options = IrOptions {
            fwhm,
            step,
            ..IrOptions::default()
        };
        let error = ir_spectrum(&structure, &modes, &options).unwrap_err();
        assert_eq!(error.code(), "invalid_input");
    }
}

#[test]
fn gamma_phonons_expose_three_acoustic_modes() {
    let structure = diamond();
    let calculator = UffCalculator::with_defaults().unwrap();
    let result = phonons(&calculator, &structure, &PhononOptions::default()).unwrap();

    assert_eq!(result.q_points, vec![[0.0, 0.0, 0.0]]);
    assert_eq!(result.gamma_frequencies().len(), 24);
    assert_eq!(result.modes.zero_modes.len(), 3, "three acoustic modes");
    assert!(result.note.contains("Gamma"));

    let error = phonons(
        &calculator,
        &structure,
        &PhononOptions {
            q_grid: [2, 1, 1],
            ..PhononOptions::default()
        },
    )
    .unwrap_err();
    assert_eq!(error.code(), "invalid_input");
    assert!(error.message().contains("not implemented"));
}

#[test]
fn invalid_hessian_options_are_rejected() {
    let structure = carbon_dimer().0;
    let calculator = UffCalculator::with_defaults().unwrap();
    for step in [0.0, -0.01, f64::NAN, f64::INFINITY] {
        let options = HessianOptions {
            step,
            displacement: DisplacementMode::Molecule,
        };
        let error = hessian(&calculator, &structure, &options).unwrap_err();
        assert_eq!(error.code(), "invalid_input", "step {step} was accepted");
    }
    assert!(HessianOptions::with_step(0.0).is_err());
}

#[test]
fn non_p1_structures_are_rejected() {
    let structure = carbon_dimer().0;
    let mut builder = StructureBuilder::new("not-p1", *structure.lattice()).hall("P 2");
    for site in structure.sites() {
        builder = builder.site(site.clone());
    }
    builder = builder.operation(xtalforge_core::SymOp::identity());
    let structure = builder.build().unwrap();
    let calculator = UffCalculator::with_defaults().unwrap();
    let error = hessian(&calculator, &structure, &HessianOptions::molecule()).unwrap_err();
    assert_eq!(error.code(), "invalid_input");
}

struct NanForces;

impl Calculator for NanForces {
    fn name(&self) -> &'static str {
        "nan-forces"
    }

    fn energy(&self, _structure: &Structure) -> Result<f64, Error> {
        Ok(0.0)
    }

    fn energy_forces(&self, structure: &Structure) -> Result<(f64, Vec<[f64; 3]>), Error> {
        Ok((0.0, vec![[f64::NAN, 0.0, 0.0]; structure.sites().len()]))
    }
}

#[test]
fn non_finite_forces_are_rejected() {
    let structure = carbon_dimer().0;
    let error = hessian(&NanForces, &structure, &HessianOptions::molecule()).unwrap_err();
    assert_eq!(error.code(), "invalid_input");
    assert!(error.message().contains("non-finite"));
}
