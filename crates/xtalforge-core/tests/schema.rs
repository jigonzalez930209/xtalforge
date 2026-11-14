//! Integration tests for the versioned structure schema (Task 1.5).

#![allow(
    clippy::unwrap_used,
    clippy::expect_used,
    clippy::panic,
    clippy::indexing_slicing
)]

use xtalforge_core::{MagneticAxis, SCHEMA_VERSION, Site, StructureDocument};

const DOCUMENT: &str = include_str!("../../../tests/schema/structure-v1.json");

#[test]
fn fixture_loads_and_round_trips() {
    let document = StructureDocument::from_json(DOCUMENT).unwrap();
    assert_eq!(document.schema_version, SCHEMA_VERSION);
    assert_eq!(document.structure.name(), "NaCl");
    assert_eq!(document.structure.formula(), "Na Cl");
    let text = document.to_json().unwrap();
    let back = StructureDocument::from_json(&text).unwrap();
    assert_eq!(back, document);
}

#[test]
fn unknown_top_level_fields_are_rejected() {
    let text = DOCUMENT.replace(
        "\"schema_version\": 1,",
        "\"schema_version\": 1, \"unexpected\": true,",
    );
    assert!(StructureDocument::from_json(&text).is_err());
}

#[test]
fn future_schema_versions_are_rejected() {
    let text = DOCUMENT.replace("\"schema_version\": 1", "\"schema_version\": 2");
    let error = StructureDocument::from_json(&text).expect_err("future schema must fail");
    assert!(error.message().contains("unsupported"));
}

#[test]
fn structure_creates_a_document() {
    let document = StructureDocument::from_json(DOCUMENT).unwrap();
    let wrapped = document.structure.to_document();
    assert_eq!(wrapped.schema_version, SCHEMA_VERSION);
    assert_eq!(wrapped.structure, document.structure);
}

#[test]
fn magnetic_moment_is_an_additive_schema_field() {
    let document = StructureDocument::from_json(DOCUMENT).unwrap();
    assert!(
        document
            .structure
            .sites()
            .iter()
            .all(|site| site.magnetic_moment().is_none())
    );
    let plain = document.structure.sites()[0].clone();
    let older: Site = serde_json::from_str(&serde_json::to_string(&plain).unwrap()).unwrap();
    assert!(older.magnetic_moment().is_none());

    let magnetic = plain.with_collinear_moment(MagneticAxis::C, 3.0).unwrap();
    let text = serde_json::to_string(&magnetic).unwrap();
    assert!(text.contains("\"magnetic_moment\""));
    assert!(text.contains("\"collinear\""));
    let back: Site = serde_json::from_str(&text).unwrap();
    assert_eq!(back, magnetic);
    assert_eq!(
        back.magnetic_moment().unwrap().components(),
        [0.0, 0.0, 3.0]
    );
    assert_eq!(SCHEMA_VERSION, 1);
}
