//! Integration tests for sites and structures (Task 1.4).

#![allow(
    clippy::unwrap_used,
    clippy::expect_used,
    clippy::panic,
    clippy::indexing_slicing
)]

use serde_json::json;
use xtalforge_core::{Lattice, Site, Structure, StructureBuilder};

#[test]
fn site_normalizes_elements_and_keeps_props() {
    let site = Site::new("Fe1", "Fe2+", [0.0, 0.0, 0.0])
        .unwrap()
        .with_occupancy(0.5)
        .unwrap()
        .with_u_iso(0.01)
        .unwrap()
        .with_charge(2.0)
        .unwrap()
        .with_disorder_group(3)
        .with_prop("uff", json!("Fe3+"));
    assert_eq!(site.label(), "Fe1");
    assert_eq!(site.element(), "Fe");
    assert_eq!(site.occupancy(), 0.5);
    assert!(site.is_partial());
    assert_eq!(site.u_iso(), Some(0.01));
    assert_eq!(site.charge(), Some(2.0));
    assert_eq!(site.disorder_group(), Some(3));
    assert_eq!(site.props()["uff"], json!("Fe3+"));
}

#[test]
fn site_validation_rejects_bad_values() {
    assert!(Site::new("", "C", [0.0, 0.0, 0.0]).is_err());
    assert!(Site::new("C1", "Zz", [0.0, 0.0, 0.0]).is_err());
    assert!(Site::new("C1", "C", [f64::NAN, 0.0, 0.0]).is_err());
    let site = Site::new("C1", "C", [0.0, 0.0, 0.0]).unwrap();
    assert!(site.clone().with_occupancy(1.5).is_err());
    assert!(site.clone().with_u_iso(-1.0).is_err());
    assert!(
        site.with_u_aniso([f64::INFINITY, 0.0, 0.0, 0.0, 0.0, 0.0])
            .is_err()
    );
}

#[test]
fn structure_properties() {
    let lattice = Lattice::cubic(5.64).unwrap();
    let structure = StructureBuilder::new("NaCl", lattice)
        .site(Site::new("Na1", "Na", [0.0, 0.0, 0.0]).unwrap())
        .site(Site::new("Cl1", "Cl", [0.5, 0.5, 0.5]).unwrap())
        .build()
        .unwrap();
    assert_eq!(structure.formula(), "Na Cl");
    assert!(structure.is_p1());
    let mass = structure.cell_mass().unwrap();
    assert!((mass - 58.442_769_77).abs() < 1e-6);
    let density = structure.density().unwrap();
    assert!((density - 0.541).abs() < 0.001);
    assert_eq!(structure.z(), 1);
    assert_eq!(structure.total_charge(), 0.0);
    assert_eq!(structure.cell_contents().len(), 2);
}

#[test]
fn partial_occupancy_changes_the_formula() {
    let lattice = Lattice::cubic(4.0).unwrap();
    let structure = StructureBuilder::new("partial", lattice)
        .site(
            Site::new("C1", "C", [0.0, 0.0, 0.0])
                .unwrap()
                .with_occupancy(0.5)
                .unwrap(),
        )
        .site(Site::new("O1", "O", [0.5, 0.5, 0.5]).unwrap())
        .build()
        .unwrap();
    assert_eq!(structure.formula(), "C0.500 O");
    assert!((structure.total_charge()).abs() < 1e-12);
}

#[test]
fn molecule_is_a_padded_p1_box() {
    let structure = Structure::molecule(
        "H2",
        &[
            ("H".to_owned(), [0.0, 0.0, 0.0]),
            ("H".to_owned(), [1.0, 0.0, 0.0]),
        ],
        5.0,
    )
    .unwrap();
    assert!(structure.is_p1());
    assert_eq!(structure.sites().len(), 2);
    let [a, b, c] = structure.lattice().lengths();
    assert!((a - 11.0).abs() < 1e-12);
    assert!((b - 10.0).abs() < 1e-12);
    assert!((c - 10.0).abs() < 1e-12);
    assert_eq!(structure.formula(), "H2");
}

#[test]
fn summary_and_metadata() {
    let lattice = Lattice::cubic(4.0).unwrap();
    let structure = StructureBuilder::new("diamond", lattice)
        .site(Site::new("C1", "C", [0.0, 0.0, 0.0]).unwrap())
        .metadata("source", json!("test"))
        .build()
        .unwrap();
    let summary = structure.summary().unwrap();
    assert_eq!(summary.name, "diamond");
    assert_eq!(summary.n_sites, 1);
    assert!((summary.volume - 64.0).abs() < 1e-12);
    assert_eq!(structure.metadata()["source"], json!("test"));
    let renamed = structure.with_name("renamed");
    assert_eq!(renamed.name(), "renamed");
}

#[test]
fn builder_rejects_an_empty_hall_symbol() {
    let lattice = Lattice::cubic(4.0).unwrap();
    let error = StructureBuilder::new("x", lattice)
        .hall("  ")
        .build()
        .expect_err("empty Hall symbol must fail");
    assert_eq!(error.code(), "invalid_input");
}

#[test]
fn shared_structure_vector_summary() {
    let summary =
        xtalforge_core::structure_info(include_str!("../../../tests/vectors/structure.json"))
            .unwrap();
    assert_eq!(summary.name, "NaCl");
    assert_eq!(summary.formula, "Na Cl");
    assert_eq!(summary.n_sites, 2);
    assert!((summary.volume - 179.406_143_999_999_95).abs() < 1e-9);
    assert!((summary.density - 0.540_931_879_932_264_1).abs() < 1e-9);
    assert!(summary.is_p1);
}

#[test]
fn structure_info_accepts_the_versioned_envelope() {
    let lattice = Lattice::cubic(4.0).unwrap();
    let structure = StructureBuilder::new("C", lattice)
        .site(Site::new("C1", "C", [0.0, 0.0, 0.0]).unwrap())
        .build()
        .unwrap();
    let bare = serde_json::to_string(&structure).unwrap();
    let wrapped = structure.to_document().to_json().unwrap();
    let from_bare = xtalforge_core::structure_info(&bare).unwrap();
    let from_wrapped = xtalforge_core::structure_info(&wrapped).unwrap();
    assert_eq!(from_bare, from_wrapped);
    assert_eq!(from_wrapped.formula, "C");
    let wrong_version = wrapped.replace("\"schema_version\":1", "\"schema_version\":99");
    assert_eq!(
        xtalforge_core::structure_info(&wrong_version)
            .unwrap_err()
            .code(),
        "invalid_input"
    );
}

#[test]
fn structure_round_trips_through_json() {
    let lattice = Lattice::cubic(4.0).unwrap();
    let structure = StructureBuilder::new("C", lattice)
        .site(Site::new("C1", "C", [0.0, 0.0, 0.0]).unwrap())
        .build()
        .unwrap();
    let text = serde_json::to_string(&structure).unwrap();
    let back: Structure = serde_json::from_str(&text).unwrap();
    assert_eq!(back, structure);
}
