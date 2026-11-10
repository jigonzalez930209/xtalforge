//! Scalar registry (Task 8.3.1).

#![allow(
    clippy::unwrap_used,
    clippy::expect_used,
    clippy::panic,
    clippy::indexing_slicing
)]

use xtalforge_core::{
    Error, Lattice, ScalarInfo, ScalarRegistry, Site, Structure, StructureBuilder,
};

fn nacl() -> Structure {
    StructureBuilder::new("nacl", Lattice::cubic(2.5).unwrap())
        .site(Site::new("Na", "Na", [0.0, 0.0, 0.0]).unwrap())
        .site(Site::new("Cl", "Cl", [0.5, 0.5, 0.5]).unwrap())
        .build()
        .unwrap()
}

fn first_x(structure: &Structure) -> Result<Vec<f64>, Error> {
    Ok(structure
        .sites()
        .iter()
        .map(|site| site.frac()[0])
        .collect())
}

#[test]
fn builtin_scalars_compute() {
    let registry = ScalarRegistry::new();
    let structure = nacl();
    for name in [
        "bond_lengths",
        "coordination_numbers",
        "atom_charges",
        "atom_masses",
    ] {
        assert!(registry.info(name).is_some(), "{name}");
        let values = registry.compute(name, &structure).unwrap();
        assert!(!values.is_empty(), "{name}");
        assert!(
            values
                .iter()
                .all(|value| !value.is_nan() || name == "atom_charges")
        );
    }
    let coordination = registry
        .compute("coordination_numbers", &structure)
        .unwrap();
    assert_eq!(coordination, vec![8.0, 8.0]);
    assert!(registry.names().contains(&"bond_lengths"));
}

#[test]
fn a_plugin_scalar_registers_and_runs() {
    let mut registry = ScalarRegistry::new();
    registry
        .register(
            ScalarInfo {
                name: "first_x".to_owned(),
                description: "First fractional coordinate".to_owned(),
                unit: String::new(),
            },
            first_x,
        )
        .unwrap();
    assert_eq!(
        registry.compute("first_x", &nacl()).unwrap(),
        vec![0.0, 0.5]
    );
    assert_eq!(registry.info("first_x").unwrap().unit, "");
}

#[test]
fn duplicate_and_unknown_names_are_errors() {
    let mut registry = ScalarRegistry::empty();
    let info = ScalarInfo {
        name: "first_x".to_owned(),
        description: String::new(),
        unit: String::new(),
    };
    registry.register(info.clone(), first_x).unwrap();
    assert!(registry.register(info, first_x).is_err());
    assert!(registry.compute("missing", &nacl()).is_err());
    let empty = ScalarInfo {
        name: " ".to_owned(),
        description: String::new(),
        unit: String::new(),
    };
    assert!(registry.register(empty, first_x).is_err());
}
