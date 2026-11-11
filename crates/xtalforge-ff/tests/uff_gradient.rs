//! Analytic UFF forces, the Gamma-point Hessian and the silicon modes
//! against central differences of the energy (functional audit 2026-11-11).
//! Silicon is read through a Cartesian round trip, which used to leave a
//! -1e-32 coordinate that the P1 expansion turned into 1.0.

#![allow(
    clippy::unwrap_used,
    clippy::expect_used,
    clippy::panic,
    clippy::indexing_slicing
)]

use xtalforge_core::Structure;
use xtalforge_ff::UffModel;

const SILICON: &str = include_str!("../../../docs/public/samples/Si.cif");

fn check(structure: &Structure, label: &str) {
    let model = UffModel::from_structure(structure).unwrap();
    let mut base: Vec<[f64; 3]> = structure
        .p1()
        .unwrap()
        .atoms
        .iter()
        .map(|a| structure.lattice().to_cart(a.frac))
        .collect();
    // A deterministic, non-symmetric perturbation away from the reference.
    for (index, point) in base.iter_mut().enumerate() {
        for (axis, value) in point.iter_mut().enumerate() {
            *value += 0.03 * (((index * 3 + axis) as f64 * 1.7).sin());
        }
    }
    let (_, forces) = model.energy_forces_at(&base).unwrap();
    let h = 1e-5;
    let mut worst: f64 = 0.0;
    for atom in 0..base.len() {
        for axis in 0..3 {
            let mut plus = base.clone();
            plus[atom][axis] += h;
            let mut minus = base.clone();
            minus[atom][axis] -= h;
            let numeric =
                -(model.energy_at(&plus).unwrap() - model.energy_at(&minus).unwrap()) / (2.0 * h);
            worst = worst.max((numeric - forces[atom][axis]).abs());
        }
    }
    assert!(worst < 1e-3, "{label}: max force error {worst} eV/A");
}

#[test]
fn periodic_silicon_forces_match_the_energy() {
    let structure = xtalforge_io::cif::read_cif(SILICON).unwrap();
    check(&structure, "Si");
}

#[test]
fn rebuilt_calculator_energy_is_smooth_in_silicon() {
    use xtalforge_ff::{Calculator, UffCalculator};
    let structure = xtalforge_io::cif::read_cif(SILICON)
        .unwrap()
        .as_p1()
        .unwrap();
    let calculator = UffCalculator::default();
    let lattice = *structure.lattice();
    let energy = |shift: f64| {
        let mut sites = structure.sites().to_vec();
        let mut cart = lattice.to_cart(sites[0].frac());
        cart[0] += shift;
        sites[0] = sites[0].clone().with_frac(lattice.to_frac(cart)).unwrap();
        let mut builder = xtalforge_core::StructureBuilder::new("si", lattice);
        for site in sites {
            builder = builder.site(site);
        }
        let moved = builder.build().unwrap();
        calculator.energy_forces(&moved).unwrap().0
    };
    let values: Vec<f64> = [-0.02, -0.01, 0.0, 0.01, 0.02]
        .iter()
        .map(|s| energy(*s))
        .collect();
    let curvature = (values[1] - 2.0 * values[2] + values[3]) / 1e-4;
    let wide = (values[0] - 2.0 * values[2] + values[4]) / 4e-4;
    assert!(
        (curvature - wide).abs() < 0.05 * curvature.abs().max(1.0),
        "{values:?} {curvature} {wide}"
    );
}

#[test]
fn silicon_hessian_is_symmetric_and_matches_the_energy() {
    use xtalforge_ff::vibrations::{HessianOptions, hessian};
    use xtalforge_ff::{Calculator, UffCalculator};
    let structure = xtalforge_io::cif::read_cif(SILICON)
        .unwrap()
        .as_p1()
        .unwrap();
    let calculator = UffCalculator::default();
    let matrix = hessian(&calculator, &structure, &HessianOptions::crystal()).unwrap();
    let size = matrix.size();
    let mut asymmetry: f64 = 0.0;
    for row in 0..size {
        for column in 0..size {
            asymmetry =
                asymmetry.max((matrix.matrix[row][column] - matrix.matrix[column][row]).abs());
        }
    }
    let lattice = *structure.lattice();
    let energy = |shift: f64| {
        let mut builder = xtalforge_core::StructureBuilder::new("si", lattice);
        for (index, site) in structure.sites().iter().enumerate() {
            let mut cart = lattice.to_cart(site.frac());
            if index == 0 {
                cart[0] += shift;
            }
            builder = builder.site(site.clone().with_frac(lattice.to_frac(cart)).unwrap());
        }
        calculator
            .energy_forces(&builder.build().unwrap())
            .unwrap()
            .0
    };
    let curvature = (energy(0.01) - 2.0 * energy(0.0) + energy(-0.01)) / 1e-4;
    assert!(
        (matrix.matrix[0][0] - curvature).abs() < 0.05 * curvature.abs().max(1.0)
            && asymmetry < 0.5,
        "H00 {} curvature {curvature} asymmetry {asymmetry}",
        matrix.matrix[0][0]
    );
}

#[test]
fn silicon_gamma_modes_are_physical() {
    use xtalforge_ff::UffCalculator;
    use xtalforge_ff::vibrations::{HessianOptions, hessian, vibrational_modes};
    let structure = xtalforge_io::cif::read_cif(SILICON)
        .unwrap()
        .as_p1()
        .unwrap();
    let calculator = UffCalculator::with_defaults().unwrap();
    let options = HessianOptions::crystal();
    let matrix = hessian(&calculator, &structure, &options).unwrap();
    let modes = vibrational_modes(&matrix, &structure, &options).unwrap();
    let low = modes
        .frequencies
        .iter()
        .copied()
        .fold(f64::INFINITY, f64::min);
    let high = modes
        .frequencies
        .iter()
        .copied()
        .fold(f64::NEG_INFINITY, f64::max);
    assert!(low > -50.0 && high < 1500.0, "{:?}", modes.frequencies);
}

#[test]
fn calculator_forces_match_energy_differences_off_symmetry() {
    use xtalforge_ff::{Calculator, UffCalculator};
    let structure = xtalforge_io::cif::read_cif(SILICON)
        .unwrap()
        .as_p1()
        .unwrap();
    let calculator = UffCalculator::with_defaults().unwrap();
    let lattice = *structure.lattice();
    let build = |shift: [f64; 3], atom: usize, axis: usize, h: f64| {
        let mut builder = xtalforge_core::StructureBuilder::new("si", lattice);
        for (index, site) in structure.sites().iter().enumerate() {
            let mut cart = lattice.to_cart(site.frac());
            if index == 0 {
                for (value, delta) in cart.iter_mut().zip(shift) {
                    *value += delta;
                }
            }
            if index == atom {
                cart[axis] += h;
            }
            builder = builder.site(site.clone().with_frac(lattice.to_frac(cart)).unwrap());
        }
        builder.build().unwrap()
    };
    let shift = [0.01, 0.0, 0.0];
    let (_, forces) = calculator.energy_forces(&build(shift, 0, 0, 0.0)).unwrap();
    let mut worst: f64 = 0.0;
    for (atom, force) in forces.iter().enumerate() {
        for (axis, component) in force.iter().enumerate() {
            let plus = calculator.energy(&build(shift, atom, axis, 1e-5)).unwrap();
            let minus = calculator.energy(&build(shift, atom, axis, -1e-5)).unwrap();
            let numeric = -(plus - minus) / 2e-5;
            let error = (numeric - component).abs();
            worst = worst.max(error);
        }
    }
    assert!(worst < 1e-3, "worst {worst}");
}
