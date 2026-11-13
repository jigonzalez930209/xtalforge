//! Bit-identical parallel UFF energy and forces (Task 15.2.3).

#![allow(
    clippy::unwrap_used,
    clippy::expect_used,
    clippy::panic,
    clippy::indexing_slicing
)]

use xtalforge_core::ParallelOptions;
use xtalforge_ff::{CalculatorOptions, UffModel};
use xtalforge_io::read_cif;

fn mof5() -> xtalforge_core::Structure {
    read_cif(include_str!("../../../tests/corpus/cod/MOF-5.cif")).unwrap()
}

fn assert_identical(parallel: (f64, Vec<[f64; 3]>), sequential: (f64, Vec<[f64; 3]>)) -> usize {
    assert_eq!(parallel.0.to_bits(), sequential.0.to_bits(), "energy");
    assert_eq!(parallel.1.len(), sequential.1.len());
    for (index, (left, right)) in parallel.1.iter().zip(sequential.1.iter()).enumerate() {
        for axis in 0..3 {
            assert_eq!(
                left[axis].to_bits(),
                right[axis].to_bits(),
                "force {index} axis {axis}"
            );
        }
    }
    parallel.1.len()
}

#[test]
fn parallel_forces_are_bit_identical_on_a_mof() {
    let structure = mof5();
    let model = UffModel::build(&structure, CalculatorOptions::default()).unwrap();
    assert!(model.n_vdw_pairs() > 256, "the test must cross a block");
    let parallel = model.energy_forces().unwrap();
    let sequential = model.energy_forces_sequential().unwrap();
    let atoms = assert_identical(parallel, sequential);
    assert_eq!(atoms, model.n_p1_atoms());
}

#[test]
fn thread_count_does_not_change_a_bit() {
    let structure = mof5();
    let model = UffModel::build(&structure, CalculatorOptions::default()).unwrap();
    let reference = model.energy_forces_sequential().unwrap();
    for threads in [1_usize, 2, 4] {
        let result = ParallelOptions::new(Some(threads))
            .install(|| model.energy_forces().unwrap())
            .unwrap();
        assert_identical(result, reference.clone());
    }
}

#[test]
fn the_sequential_reference_matches_the_default_at_other_positions() {
    let structure = mof5();
    let model = UffModel::build(&structure, CalculatorOptions::default()).unwrap();
    let mut moved = model.positions().to_vec();
    for (index, point) in moved.iter_mut().enumerate() {
        point[0] += (index % 7) as f64 * 1e-4;
        point[1] -= (index % 5) as f64 * 1e-4;
    }
    let parallel = model.energy_forces_at(&moved).unwrap();
    let sequential = model.energy_forces_at_sequential(&moved).unwrap();
    assert_identical(parallel, sequential);
}

#[test]
fn a_corpus_framework_also_matches() {
    let structure = read_cif(include_str!("../../../tests/corpus/cod/MIL-53.cif")).unwrap();
    let model = UffModel::build(&structure, CalculatorOptions::default()).unwrap();
    let parallel = model.energy_forces().unwrap();
    let sequential = model.energy_forces_sequential().unwrap();
    assert_identical(parallel, sequential);
}
