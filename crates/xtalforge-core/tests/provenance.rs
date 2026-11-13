//! Integration tests for provenance records (Task 16.2.1).

#![allow(
    clippy::unwrap_used,
    clippy::expect_used,
    clippy::panic,
    clippy::indexing_slicing
)]

use std::collections::BTreeMap;

use serde_json::Value;
use xtalforge_core::{
    Lattice, PROVENANCE_METADATA_KEY, PROVENANCE_SCHEMA_VERSION, Provenance, ProvenanceInput, Site,
    StructureBuilder, sha256_hex,
};

fn structure() -> xtalforge_core::Structure {
    StructureBuilder::new("Si", Lattice::cubic(5.43).unwrap())
        .site(Site::new("Si1", "Si", [0.0, 0.0, 0.0]).unwrap())
        .build()
        .unwrap()
}

fn record() -> Provenance {
    Provenance::new("pxrd")
        .with_option("wavelength", 1.5406)
        .with_option("source", "Cu")
        .with_input_bytes(Some("sample.cif".to_owned()), b"crystal")
        .with_bundle("elements", "data-2026.11.13")
        .with_seed(7)
        .with_timestamp("2026-11-13T12:00:00Z")
}

#[test]
fn sha256_matches_known_digests() {
    assert_eq!(
        sha256_hex(b""),
        "e3b0c44298fc1c149afbf4c8996fb92427ae41e4649b934ca495991b7852b855"
    );
    assert_eq!(
        sha256_hex(b"abc"),
        "ba7816bf8f01cfea414140de5dae2223b00361a396177a9cb410ff61f20015ad"
    );
}

#[test]
fn record_round_trips_with_a_stable_schema_version() {
    let record = record();
    assert_eq!(record.schema_version, PROVENANCE_SCHEMA_VERSION);
    assert_eq!(record.xtalforge_version, xtalforge_core::version());
    let text = record.to_json().unwrap();
    let back = Provenance::from_json(&text).unwrap();
    assert_eq!(back, record);
}

#[test]
fn json_is_deterministic_and_ordered() {
    let mut first = BTreeMap::new();
    first.insert("beta".to_owned(), Value::from(2));
    first.insert("alpha".to_owned(), Value::from(1));
    let mut second = BTreeMap::new();
    second.insert("alpha".to_owned(), Value::from(1));
    second.insert("beta".to_owned(), Value::from(2));
    let left = Provenance::new("op").with_options(first).to_json().unwrap();
    let right = Provenance::new("op")
        .with_options(second)
        .to_json()
        .unwrap();
    assert_eq!(left, right);
    let alpha = left.find("\"alpha\"").unwrap();
    let beta = left.find("\"beta\"").unwrap();
    assert!(alpha < beta);
    assert!(left.starts_with(&format!(
        "{{\"schema_version\":{PROVENANCE_SCHEMA_VERSION},\"xtalforge_version\":\"{}\",\"operation\":\"op\"",
        xtalforge_core::version()
    )));
}

#[test]
fn optional_fields_are_omitted_when_absent() {
    let text = Provenance::new("identity").to_json().unwrap();
    assert!(!text.contains("\"seed\""));
    assert!(!text.contains("\"timestamp\""));
    assert!(!text.contains("\"options\""));
    assert!(!text.contains("\"inputs\""));
    assert!(!text.contains("\"bundles\""));
}

#[test]
fn attach_and_read_structure_metadata() {
    let record = record();
    let structure = record.attach(&structure()).unwrap();
    assert!(structure.metadata().contains_key(PROVENANCE_METADATA_KEY));
    let back = Provenance::from_structure(&structure).unwrap().unwrap();
    assert_eq!(back, record);

    let replacement = Provenance::new("relaxed");
    let replaced = replacement.attach(&structure).unwrap();
    assert_eq!(
        Provenance::from_structure(&replaced).unwrap().unwrap(),
        replacement
    );
    assert!(
        Provenance::from_structure(&replaced.without_metadata(PROVENANCE_METADATA_KEY))
            .unwrap()
            .is_none()
    );
}

#[test]
fn provenance_metadata_accepts_a_json_object() {
    let record = record();
    let value: Value = serde_json::from_str(&record.to_json().unwrap()).unwrap();
    let structure = structure().with_metadata(PROVENANCE_METADATA_KEY, value);
    assert_eq!(
        Provenance::from_structure(&structure).unwrap().unwrap(),
        record
    );
}

#[test]
fn invalid_records_are_rejected() {
    assert!(Provenance::from_json("not json").is_err());
    assert!(Provenance::from_json("{}").is_err());

    let mut future = record();
    future.schema_version = PROVENANCE_SCHEMA_VERSION + 1;
    let error = future.to_json().expect_err("future schema must fail");
    assert!(error.message().contains("schema_version"));

    let mut tampered = record().to_json().unwrap();
    tampered = tampered.replace(
        "\"operation\":\"pxrd\"",
        "\"operation\":\"pxrd\",\"extra\":1",
    );
    assert!(Provenance::from_json(&tampered).is_err());

    let mut empty = Provenance::new("op");
    empty.operation = "  ".to_owned();
    assert!(empty.to_json().is_err());

    let mut blank = Provenance::new("op");
    blank.xtalforge_version.clear();
    assert!(blank.to_json().is_err());

    assert!(ProvenanceInput::new(None, "abc").is_err());
    assert!(
        Provenance::new("op")
            .try_with_input_sha256(None, "xyz")
            .is_err()
    );

    let malformed = structure().with_metadata(PROVENANCE_METADATA_KEY, Value::from(42));
    assert!(Provenance::from_structure(&malformed).is_err());
}

#[test]
fn tampered_inputs_are_detected() {
    let record = Provenance::new("read").with_input_bytes(Some("a.cif".to_owned()), b"original");
    record.verify_input(Some("a.cif"), b"original").unwrap();
    let error = record
        .verify_input(Some("a.cif"), b"modified")
        .expect_err("a modified input must fail");
    assert!(error.message().contains("does not match"));
    assert!(record.verify_input(Some("missing"), b"original").is_err());
    assert!(record.verify_input(None, b"original").is_err());
}

#[test]
fn bundle_entries_reject_empty_names_and_versions() {
    let mut record = Provenance::new("op");
    record.bundles.insert(String::new(), "v1".to_owned());
    assert!(record.to_json().is_err());
    let mut record = Provenance::new("op");
    record.bundles.insert("elements".to_owned(), String::new());
    assert!(record.to_json().is_err());
}
