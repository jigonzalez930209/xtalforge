//! Integration tests for FIRE, L-BFGS, applying geometry and the EOS fit
//! (Task 4.4).
//!
//! The molecular cases start from the reference UFF equilibrium and distort it
//! by a known amount, so the energy minimum is known to be zero by
//! construction.  The EOS case fits points generated from a known
//! Birch-Murnaghan curve.

#![allow(
    clippy::unwrap_used,
    clippy::expect_used,
    clippy::panic,
    clippy::indexing_slicing
)]

use xtalforge_core::{Lattice, Site, Structure, StructureBuilder};
use xtalforge_ff::uff::natural_bond_length;
use xtalforge_ff::{
    Calculator, CalculatorOptions, EosForm, OptimizeOptions, RelaxMode, UffCalculator,
    apply_optimized_geometry, birch_murnaghan, fire, fit_eos, lbfgs, murnaghan,
};

fn molecule(name: &str, atoms: &[(&str, [f64; 3])]) -> Structure {
    let points: Vec<(String, [f64; 3])> = atoms
        .iter()
        .map(|(element, point)| ((*element).to_owned(), *point))
        .collect();
    Structure::molecule(name, &points, 8.0).expect("molecule")
}

fn distorted_methane() -> Structure {
    let length = natural_bond_length("C_3", "H_", 1.0).unwrap() * 1.05;
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

fn carbon_dimer() -> Structure {
    molecule(
        "carbon-dimer",
        &[("C", [0.0, 0.0, 0.0]), ("C", [1.70, 0.0, 0.0])],
    )
}

fn argon_cell() -> Structure {
    let lattice = Lattice::cubic(6.0).unwrap();
    let mut builder = StructureBuilder::new("argon-cell", lattice);
    for (index, point) in [[0.0, 0.0, 0.0], [0.25, 0.0, 0.0]].iter().enumerate() {
        builder = builder.site(Site::new(format!("Ar{index}"), "Ar", *point).unwrap());
    }
    builder.build().unwrap()
}

fn options() -> OptimizeOptions {
    OptimizeOptions {
        max_steps: 500,
        force_tolerance: 0.01,
        energy_tolerance: 1e-8,
        cell: false,
        fixed_atoms: Vec::new(),
        mode: RelaxMode::P1,
        max_step: 0.2,
    }
}

#[test]
fn fire_relaxes_distorted_methane_to_the_uff_minimum() {
    let calculator = UffCalculator::with_defaults().unwrap();
    let structure = distorted_methane();
    let initial = calculator.energy(&structure).unwrap();
    let result = fire(&calculator, &structure, &options()).unwrap();
    assert!(result.converged, "FIRE did not converge");
    assert!(result.max_force <= 0.01, "max force {}", result.max_force);
    assert!(result.energy < initial - 0.1);
    assert!(result.energy.abs() < 1e-3, "final energy {}", result.energy);
    let lattice = result.structure.lattice();
    let cart = result.structure.cart_sites();
    for index in 1..5 {
        let distance = lattice.distance(
            result.structure.sites()[0].frac(),
            result.structure.sites()[index].frac(),
        );
        assert!(
            (distance - 1.109_400_794_877_744).abs() < 0.01,
            "C-H {index} is {distance}"
        );
    }
    assert!(cart[1][0].is_finite());
}

#[test]
fn lbfgs_is_not_slower_than_fire_on_a_dimer() {
    let calculator = UffCalculator::new(CalculatorOptions {
        vdw: false,
        ..CalculatorOptions::default()
    })
    .unwrap();
    let structure = carbon_dimer();
    let mut settings = options();
    settings.max_steps = 300;
    let fire_result = fire(&calculator, &structure, &settings).unwrap();
    let lbfgs_result = lbfgs(&calculator, &structure, &settings).unwrap();
    assert!(fire_result.converged);
    assert!(lbfgs_result.converged);
    assert!(
        lbfgs_result.steps <= fire_result.steps,
        "L-BFGS {} steps, FIRE {}",
        lbfgs_result.steps,
        fire_result.steps
    );
    assert!((lbfgs_result.energy - fire_result.energy).abs() < 1e-5);
}

#[test]
fn optimized_geometry_applies_back_to_the_input() {
    let calculator = UffCalculator::with_defaults().unwrap();
    let original = distorted_methane();
    let charged = Structure::molecule(
        "charged-methane",
        &original
            .sites()
            .iter()
            .map(|site| (site.element().to_owned(), site.frac()))
            .collect::<Vec<_>>(),
        8.0,
    )
    .unwrap();
    let first = charged.sites()[0].clone().with_charge(-0.25).unwrap();
    let mut builder = StructureBuilder::new(charged.name(), *charged.lattice());
    builder = builder.site(first);
    for site in charged.sites().iter().skip(1) {
        builder = builder.site(site.clone());
    }
    let charged = builder.build().unwrap();
    let result = lbfgs(&calculator, &charged, &options()).unwrap();
    let applied = apply_optimized_geometry(&result, &charged).unwrap();
    assert_eq!(applied.sites().len(), charged.sites().len());
    assert_eq!(applied.sites()[0].charge(), Some(-0.25));
    for (before, after) in charged.sites().iter().zip(applied.sites()) {
        assert_eq!(before.label(), after.label());
        assert_eq!(before.element(), after.element());
    }
    let moved = applied
        .sites()
        .iter()
        .zip(charged.sites())
        .any(|(after, before)| after.frac() != before.frac());
    assert!(moved);
}

#[test]
fn cell_relaxation_expands_a_repulsive_argon_cell() {
    let calculator = UffCalculator::with_defaults().unwrap();
    let structure = argon_cell();
    let initial_energy = calculator.energy(&structure).unwrap();
    let settings = OptimizeOptions {
        cell: true,
        fixed_atoms: vec![0, 1],
        max_steps: 200,
        force_tolerance: 0.02,
        energy_tolerance: 1e-7,
        ..options()
    };
    let result = lbfgs(&calculator, &structure, &settings).unwrap();
    assert!(result.energy < initial_energy);
    assert!(
        result.structure.lattice().volume() > structure.lattice().volume(),
        "volume {} did not grow from {}",
        result.structure.lattice().volume(),
        structure.lattice().volume()
    );
}

#[test]
fn eos_recovers_a_known_birch_murnaghan_minimum() {
    let v0 = 100.0;
    let e0 = -5.0;
    let b0_gpa = 50.0;
    let b0 = b0_gpa / 160.217_663_4;
    let b0_prime = 4.0;
    let volumes = [92.0, 96.0, 100.0, 104.0, 108.0];
    let energies: Vec<f64> = volumes
        .iter()
        .map(|volume| birch_murnaghan(*volume, e0, v0, b0, b0_prime))
        .collect();
    let fit = fit_eos(&volumes, &energies, EosForm::BirchMurnaghan).unwrap();
    assert!((fit.v0 - v0).abs() < 1e-3, "V0 {}", fit.v0);
    assert!((fit.e0 - e0).abs() < 1e-6, "E0 {}", fit.e0);
    assert!((fit.b0 - b0_gpa).abs() < 0.05, "B0 {}", fit.b0);
    assert!((fit.b0_prime - b0_prime).abs() < 0.05);
    assert!(fit.rms < 1e-6);
    assert_eq!(fit.n_points, 5);
    assert!((fit.energy(v0) - e0).abs() < 1e-6);
}

#[test]
fn eos_murnaghan_also_recovers_its_own_points() {
    let v0 = 80.0;
    let e0 = -2.0;
    let b0 = 30.0 / 160.217_663_4;
    let b0_prime = 5.5;
    let volumes = [72.0, 76.0, 80.0, 84.0, 88.0];
    let energies: Vec<f64> = volumes
        .iter()
        .map(|volume| murnaghan(*volume, e0, v0, b0, b0_prime))
        .collect();
    let fit = fit_eos(&volumes, &energies, EosForm::Murnaghan).unwrap();
    assert!((fit.v0 - v0).abs() < 1e-2);
    assert!((fit.b0 - 30.0).abs() < 0.2);
    assert!((fit.b0_prime - b0_prime).abs() < 0.05);
}

#[test]
fn invalid_optimization_inputs_error() {
    let calculator = UffCalculator::with_defaults().unwrap();
    let structure = carbon_dimer();
    let mut settings = options();
    settings.fixed_atoms = vec![99];
    assert!(fire(&calculator, &structure, &settings).is_err());
    settings.fixed_atoms = Vec::new();
    settings.force_tolerance = 0.0;
    assert!(lbfgs(&calculator, &structure, &settings).is_err());

    let lattice = Lattice::cubic(4.0).unwrap();
    let builder = StructureBuilder::new("tetragonal", lattice)
        .hall("P 4/m m m")
        .site(Site::new("C1", "C", [0.0, 0.0, 0.0]).unwrap());
    let symmetric = builder.build().unwrap();
    let mut settings = options();
    settings.mode = RelaxMode::P1;
    assert!(lbfgs(&calculator, &symmetric, &settings).is_err());
}
