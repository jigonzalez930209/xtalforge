//! Phonons of a Lennard-Jones solid: fcc argon with UFF (physical sanity).

#![allow(
    clippy::unwrap_used,
    clippy::expect_used,
    clippy::panic,
    clippy::indexing_slicing
)]

use xtalforge_core::{Lattice, Site, Structure, StructureBuilder};
use xtalforge_ff::{Calculator, CalculatorOptions, PhononOptions, UffCalculator, phonons};

fn fcc(nearest: f64) -> Structure {
    let a = nearest * 2f64.sqrt() / 2.0;
    let matrix = [[0.0, a, a], [a, 0.0, a], [a, a, 0.0]];
    StructureBuilder::new("Ar", Lattice::new(matrix).unwrap())
        .site(Site::new("Ar1", "Ar", [0.0, 0.0, 0.0]).unwrap())
        .build()
        .unwrap()
}

fn branches(calculator: &UffCalculator, structure: &Structure, q: [f64; 3]) -> Vec<f64> {
    let options = PhononOptions {
        supercell: [3, 3, 3],
        q_path: vec![q],
        ..PhononOptions::default()
    };
    phonons(calculator, structure, &options)
        .unwrap()
        .frequencies[0]
        .clone()
}

#[test]
fn fcc_argon_has_a_physical_dispersion() {
    let calculator = UffCalculator::new(CalculatorOptions::default()).unwrap();
    let mut best = (f64::MAX, 0.0);
    for step in 0..60 {
        let nearest = 3.5 + 0.01 * f64::from(step);
        let energy = calculator.energy(&fcc(nearest)).unwrap();
        if energy < best.0 {
            best = (energy, nearest);
        }
    }
    assert!(
        (3.7..3.85).contains(&best.1),
        "nearest-neighbour {}",
        best.1
    );
    let structure = fcc(best.1);

    // Gamma: three acoustic modes at zero (acoustic sum rule).
    let gamma = branches(&calculator, &structure, [0.0; 3]);
    assert!(gamma.iter().all(|value| value.abs() < 1.0), "{gamma:?}");

    // L point: doubly degenerate transverse pair below one longitudinal mode.
    let l = branches(&calculator, &structure, [0.5, 0.5, 0.5]);
    assert!(l.iter().all(|value| *value > 0.0), "{l:?}");
    assert!((l[0] - l[1]).abs() < 0.01 * l[2], "{l:?}");
    assert!(l[2] > 1.5 * l[0], "{l:?}");

    // Small q: all branches positive and increasing with q (stable crystal).
    let small = branches(&calculator, &structure, [0.05, 0.0, 0.05]);
    let larger = branches(&calculator, &structure, [0.1, 0.0, 0.1]);
    for (low, high) in small.iter().zip(larger.iter()) {
        assert!(*low > 0.0 && high > low, "{small:?} {larger:?}");
    }
}
