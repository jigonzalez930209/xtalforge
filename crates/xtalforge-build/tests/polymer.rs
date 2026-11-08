//! Integration tests for polymer packing and build records (Task 6.9).

#![allow(
    clippy::unwrap_used,
    clippy::expect_used,
    clippy::panic,
    clippy::indexing_slicing
)]

use std::collections::BTreeSet;

use xtalforge_build::{PolymerOptions, build_polymer, pack_chains, polymer_protocol};
use xtalforge_core::{Lattice, Site, Structure, StructureBuilder};

fn atom(element: &str, x: f64) -> Structure {
    StructureBuilder::new("atom", Lattice::cubic(8.0).unwrap())
        .site(Site::new("X1", element, [x, 0.5, 0.5]).unwrap())
        .build()
        .unwrap()
}

#[test]
fn packing_keeps_all_chains_and_returns_p1() {
    let chains = vec![atom("C", 0.1), atom("C", 0.2), atom("C", 0.3)];
    let packed = pack_chains(&chains, 1.0).unwrap();
    assert_eq!(packed.sites().len(), chains.len());
    assert!(packed.is_p1());
    assert!(packed.density().unwrap() > 0.0);
    let labels: BTreeSet<&str> = packed.sites().iter().map(Site::label).collect();
    assert_eq!(labels.len(), packed.sites().len());
    assert!(
        packed
            .sites()
            .iter()
            .flat_map(|site| site.frac())
            .all(|value| (0.0..1.0).contains(&value))
    );
}

#[test]
fn build_polymer_records_its_steps() {
    let options = PolymerOptions { density: 1.0 };
    let build = build_polymer("CCO", 2, &options).unwrap();
    assert_eq!(build.structure.sites().len(), 6);
    assert!(build.structure.is_p1());
    assert!(!build.steps.is_empty());
    assert!(build.steps.iter().all(|step| !step.is_empty()));
    assert!(build_polymer("CCO", 0, &options).is_err());
    assert!(build_polymer("", 1, &options).is_err());
    assert!(build_polymer("CCO", 1, &PolymerOptions { density: 0.0 }).is_err());
}

#[test]
fn invalid_density_and_empty_inputs_error() {
    let chains = vec![atom("C", 0.5)];
    assert!(pack_chains(&chains, 0.0).is_err());
    assert!(pack_chains(&chains, -1.0).is_err());
    assert!(pack_chains(&chains, f64::NAN).is_err());
    assert!(pack_chains(&[], 1.0).is_err());
}

#[test]
fn protocol_loads_documented_stages() {
    let protocol = polymer_protocol().unwrap();
    assert!(protocol.len() >= 4 && protocol.len() <= 6);
    let joined = protocol.join("\n");
    assert!(joined.contains("Minimise"));
    assert!(joined.contains("NVT"));
    assert!(joined.contains("NPT"));
    assert!(joined.contains("Production"));
    assert_eq!(protocol, polymer_protocol().unwrap());
}
