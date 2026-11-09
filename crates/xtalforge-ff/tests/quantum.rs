//! Integration tests for the in-process quantum adapter (Task 4.10).

#![allow(
    clippy::unwrap_used,
    clippy::expect_used,
    clippy::panic,
    clippy::indexing_slicing
)]

#[cfg(feature = "sci-form")]
use xtalforge_core::Structure;
#[cfg(feature = "sci-form")]
use xtalforge_ff::Calculator;
use xtalforge_ff::{CalculatorRegistry, QuantumMethod, adapter};

#[cfg(feature = "sci-form")]
fn water() -> Structure {
    Structure::molecule(
        "water",
        &[
            ("O".to_owned(), [0.0, 0.0, 0.0]),
            ("H".to_owned(), [0.96, 0.0, 0.0]),
            ("H".to_owned(), [-0.24, 0.93, 0.0]),
        ],
        5.0,
    )
    .unwrap()
}

#[cfg(not(feature = "sci-form"))]
#[test]
fn adapter_requires_the_feature() {
    let error = adapter(QuantumMethod::Gfn2).unwrap_err();
    assert!(
        error.message().contains("--features sci-form"),
        "{}",
        error.message()
    );
    let registry = CalculatorRegistry::builtin();
    assert!(registry.names().contains(&"gfn2"));
    let info = registry.info("gfn2").expect("a registry entry");
    assert!(!info.available());
}

#[cfg(feature = "sci-form")]
#[test]
fn gfn_methods_return_finite_water_energies() {
    let water = water();
    for method in [
        QuantumMethod::Gfn0,
        QuantumMethod::Gfn1,
        QuantumMethod::Gfn2,
    ] {
        let calculator = adapter(method).unwrap();
        assert_eq!(calculator.name(), method.name());
        let energy = calculator.energy(&water).unwrap();
        assert!(energy.is_finite(), "{} energy {energy}", method.name());
    }
}

#[cfg(feature = "sci-form")]
#[test]
fn gfn2_returns_conservative_forces() {
    let water = water();
    let calculator = adapter(QuantumMethod::Gfn2).unwrap();
    assert!(calculator.provides_forces());
    let (energy, forces) = calculator.energy_forces(&water).unwrap();
    assert!(energy.is_finite());
    assert_eq!(forces.len(), 3);
    assert!(
        forces
            .iter()
            .flatten()
            .all(|value| value.is_finite() && value.abs() < 100.0)
    );
    let sum = forces.iter().fold([0.0; 3], |total, force| {
        [
            total[0] + force[0],
            total[1] + force[1],
            total[2] + force[2],
        ]
    });
    assert!(
        sum.iter().all(|value| value.abs() < 0.05),
        "net force {sum:?}"
    );
    for method in [QuantumMethod::Gfn0, QuantumMethod::Gfn1] {
        let error = adapter(method).unwrap().energy_forces(&water).unwrap_err();
        assert!(error.message().contains("GFN2"));
    }
}

#[cfg(feature = "sci-form")]
#[test]
fn a_space_group_structure_is_refused() {
    let cif = include_str!("../../../tests/corpus/simple/TiO2.cif");
    let rutile = xtalforge_io::cif::read_cif(cif).unwrap();
    let error = adapter(QuantumMethod::Gfn2)
        .unwrap()
        .energy(&rutile)
        .unwrap_err();
    assert!(
        error.message().contains("molecular only"),
        "{}",
        error.message()
    );
}

#[cfg(feature = "sci-form")]
#[test]
fn the_registry_builds_the_in_process_methods() {
    let registry = CalculatorRegistry::builtin();
    for name in ["gfn0", "gfn1", "gfn2"] {
        let info = registry.info(name).expect("a registry entry");
        assert!(info.available(), "{name} must be available");
        assert!(!info.capabilities.periodic);
        assert!(info.capabilities.in_process);
    }
    let calculator = registry
        .build("gfn2", &xtalforge_ff::CalculatorOptions::default())
        .unwrap();
    assert_eq!(calculator.name(), "gfn2");
}
