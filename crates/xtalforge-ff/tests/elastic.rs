//! Integration tests for the elastic tensor (Task 9.6.1).

#![allow(
    clippy::unwrap_used,
    clippy::expect_used,
    clippy::panic,
    clippy::indexing_slicing
)]

use xtalforge_core::{Lattice, Site, StructureBuilder};
use xtalforge_ff::{
    CalculatorOptions, ElasticOptions, UffCalculator, deformation_set, elastic_tensor,
};

fn silicon() -> xtalforge_core::Structure {
    let a = 5.430941;
    StructureBuilder::new("Si", Lattice::cubic(a).unwrap())
        .site(Site::new("Si1", "Si", [0.0, 0.0, 0.0]).unwrap())
        .site(Site::new("Si2", "Si", [0.25, 0.25, 0.25]).unwrap())
        .build()
        .unwrap()
}

#[test]
fn deformation_set_has_six_tensors() {
    let set = deformation_set(0.01);
    assert_eq!(set.len(), 6);
    assert!((set[0][0][0] - 0.01).abs() < 1e-12);
    assert!((set[3][1][2] - 0.005).abs() < 1e-12);
}

#[test]
fn elastic_tensor_is_symmetric_and_stable_for_silicon() {
    let calculator = UffCalculator::new(CalculatorOptions::default()).unwrap();
    let tensor =
        elastic_tensor(&silicon(), &calculator, &ElasticOptions { strain: 0.002 }).unwrap();
    assert_eq!(tensor.evaluations, 12);
    for row in 0..6 {
        for column in 0..6 {
            assert!((tensor.stiffness[row][column] - tensor.stiffness[column][row]).abs() < 1e-8);
        }
    }
    assert!(tensor.bulk_modulus_hill.is_finite());
    assert!(tensor.bulk_modulus_hill.is_finite());
    assert!(tensor.shear_modulus_hill.is_finite());
}

#[test]
fn elastic_tensor_rejects_a_bad_strain() {
    let calculator = UffCalculator::new(CalculatorOptions::default()).unwrap();
    assert!(elastic_tensor(&silicon(), &calculator, &ElasticOptions { strain: 0.0 }).is_err());
}
