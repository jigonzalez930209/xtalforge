//! Integration tests for the JSON Lines dataset container (Task 16.3.1).

#![allow(
    clippy::unwrap_used,
    clippy::expect_used,
    clippy::panic,
    clippy::indexing_slicing
)]

use std::collections::BTreeMap;

use xtalforge_analysis::{DATASET_VERSION, Dataset, DatasetRecord, PropertyValue};
use xtalforge_core::{Lattice, Site, Structure, StructureBuilder};

fn exact_cubic(a: f64) -> Lattice {
    Lattice::new([[a, 0.0, 0.0], [0.0, a, 0.0], [0.0, 0.0, a]]).unwrap()
}

fn nacl() -> Structure {
    let positions = [
        ("Na", [0.0, 0.0, 0.0]),
        ("Na", [0.0, 0.5, 0.5]),
        ("Na", [0.5, 0.0, 0.5]),
        ("Na", [0.5, 0.5, 0.0]),
        ("Cl", [0.5, 0.5, 0.5]),
        ("Cl", [0.5, 0.0, 0.0]),
        ("Cl", [0.0, 0.5, 0.0]),
        ("Cl", [0.0, 0.0, 0.5]),
    ];
    let mut builder = StructureBuilder::new("NaCl", exact_cubic(5.64));
    for (index, (element, position)) in positions.iter().enumerate() {
        builder = builder.site(Site::new(format!("{element}{index}"), element, *position).unwrap());
    }
    builder.build().unwrap()
}

fn cscl() -> Structure {
    StructureBuilder::new("CsCl", exact_cubic(4.2))
        .site(Site::new("Cs1", "Cs", [0.0, 0.0, 0.0]).unwrap())
        .site(Site::new("Cl1", "Cl", [0.5, 0.5, 0.5]).unwrap())
        .build()
        .unwrap()
}

fn sample() -> Dataset {
    let first = DatasetRecord::new("nacl", nacl())
        .unwrap()
        .with_property("energy", -3.25)
        .unwrap()
        .with_property("label", "rock salt")
        .unwrap()
        .with_property("converged", true)
        .unwrap();
    let second = DatasetRecord::new("cscl", cscl())
        .unwrap()
        .with_property("energy", -2.5)
        .unwrap()
        .with_property("steps", 12_i64)
        .unwrap()
        .with_property("forces", vec![0.0, 0.1, -0.1])
        .unwrap();
    let dataset = Dataset::from_records("halides", vec![first, second]).unwrap();
    dataset
        .with_description("two alkali halides")
        .with_source("unit test")
        .with_label("license", "MIT")
        .unwrap()
}

#[test]
fn json_lines_round_trip_is_deterministic() {
    let dataset = sample();
    let text = dataset.to_json_lines().unwrap();
    assert_eq!(text, dataset.to_json_lines().unwrap());
    assert_eq!(text.lines().count(), 3);
    assert!(text.ends_with('\n'));
    assert!(
        text.lines()
            .next()
            .unwrap()
            .contains("\"dataset_version\":1")
    );

    let parsed = Dataset::from_json_lines(&text).unwrap();
    assert_eq!(parsed, dataset);
    assert_eq!(parsed.len(), 2);
    assert_eq!(parsed.metadata().name, "halides");
    assert_eq!(parsed.metadata().description, "two alkali halides");
    assert_eq!(parsed.metadata().source, "unit test");
    assert_eq!(
        parsed
            .metadata()
            .labels
            .get("license")
            .and_then(|value| value.as_str()),
        Some("MIT")
    );
}

#[test]
fn records_keep_their_order_and_types() {
    let dataset = sample();
    let ids: Vec<&str> = dataset.iter().map(|record| record.id.as_str()).collect();
    assert_eq!(ids, ["nacl", "cscl"]);
    assert_eq!(dataset.by_id("cscl").unwrap().structure.name(), "CsCl");
    assert!(dataset.by_id("missing").is_none());
    assert_eq!(dataset.get(1).unwrap().id, "cscl");
    assert!(dataset.get(2).is_none());

    let energy = dataset.column("energy");
    assert_eq!(energy.len(), 2);
    assert_eq!(energy[0].unwrap().as_f64(), Some(-3.25));
    assert!(dataset.column("missing")[1].is_none());

    assert_eq!(
        dataset.property_names(),
        ["converged", "energy", "forces", "label", "steps"]
    );
    let record = dataset.get(1).unwrap();
    assert_eq!(record.number("energy"), Some(-2.5));
    assert_eq!(record.number("steps"), Some(12.0));
    assert_eq!(record.property("forces").unwrap().kind(), "vector");
    assert_eq!(
        record.property("forces").unwrap().as_vector(),
        Some(&[0.0, 0.1, -0.1][..])
    );
}

#[test]
fn property_values_round_trip() {
    let mut properties = BTreeMap::new();
    properties.insert("boolean".to_owned(), PropertyValue::from(true));
    properties.insert("integer".to_owned(), PropertyValue::from(7_i64));
    properties.insert("number".to_owned(), PropertyValue::from(1.5));
    properties.insert("text".to_owned(), PropertyValue::from("hello"));
    properties.insert("vector".to_owned(), PropertyValue::from(vec![1.0, 2.0]));
    let text = serde_json::to_string(&properties).unwrap();
    let parsed: BTreeMap<String, PropertyValue> = serde_json::from_str(&text).unwrap();
    assert_eq!(parsed, properties);
    assert_eq!(parsed["integer"].kind(), "integer");
    assert_eq!(parsed["boolean"].as_bool(), Some(true));
    assert_eq!(parsed["text"].as_str(), Some("hello"));
}

#[test]
fn duplicate_identifiers_are_rejected() {
    let record = DatasetRecord::new("same", nacl()).unwrap();
    let error = Dataset::from_records("duplicates", vec![record.clone(), record]).unwrap_err();
    assert!(error.to_string().contains("duplicate"), "{error}");
}

#[test]
fn invalid_records_are_rejected() {
    assert!(DatasetRecord::new("  ", nacl()).is_err());
    let record = DatasetRecord::new("id", nacl()).unwrap();
    assert!(record.clone().with_property(" ", 1.0).is_err());
    assert!(record.with_property("energy", f64::NAN).is_err());
    assert!(
        DatasetRecord::new("id", nacl())
            .unwrap()
            .with_property("forces", vec![0.0, f64::INFINITY])
            .is_err()
    );
}

#[test]
fn malformed_records_report_their_line() {
    let text = sample().to_json_lines().unwrap();
    let mut lines: Vec<&str> = text.lines().collect();
    lines[2] = "{\"id\":\"broken\"}";
    let broken = format!("{}\n", lines.join("\n"));
    let error = Dataset::from_json_lines(&broken).unwrap_err();
    assert!(error.to_string().contains("line 3"), "{error}");
    assert!(error.to_string().contains("invalid record"), "{error}");
}

#[test]
fn an_unsupported_version_is_rejected() {
    let text = sample().to_json_lines().unwrap().replacen(
        "\"dataset_version\":1",
        "\"dataset_version\":2",
        1,
    );
    let error = Dataset::from_json_lines(&text).unwrap_err();
    assert!(
        error.to_string().contains("unsupported dataset_version"),
        "{error}"
    );
}

#[test]
fn a_truncated_file_is_rejected() {
    let text = sample().to_json_lines().unwrap();
    let mut lines: Vec<&str> = text.lines().collect();
    lines.pop();
    let truncated = format!("{}\n", lines.join("\n"));
    let error = Dataset::from_json_lines(&truncated).unwrap_err();
    assert!(error.to_string().contains("declares 2 records"), "{error}");
}

#[test]
fn empty_and_blank_inputs_are_rejected() {
    assert!(Dataset::from_json_lines("").is_err());
    assert!(Dataset::from_json_lines("\n\n").is_err());
}

#[test]
fn blank_lines_between_records_are_ignored() {
    let text = sample().to_json_lines().unwrap();
    let spaced = text.replace('\n', "\n\n");
    let parsed = Dataset::from_json_lines(&spaced).unwrap();
    assert_eq!(parsed, sample());
}

#[test]
fn an_empty_dataset_round_trips() {
    let dataset = Dataset::new("empty").with_description("no records");
    let text = dataset.to_json_lines().unwrap();
    assert_eq!(text.lines().count(), 1);
    let parsed = Dataset::from_json_lines(&text).unwrap();
    assert!(parsed.is_empty());
    assert_eq!(parsed.metadata().name, "empty");
    assert_eq!(parsed.metadata().description, "no records");
}

#[test]
fn missing_optional_metadata_is_accepted() {
    let text = "{\"dataset_version\":1,\"n_records\":0,\"metadata\":{\"name\":\"old\"}}\n";
    let dataset = Dataset::from_json_lines(text).unwrap();
    assert_eq!(dataset.metadata().name, "old");
    assert!(dataset.metadata().description.is_empty());
    assert!(dataset.metadata().source.is_empty());
    assert!(dataset.metadata().labels.is_empty());
}

#[test]
fn save_and_load_round_trip() {
    let dataset = sample();
    let stamp = std::time::SystemTime::now()
        .duration_since(std::time::UNIX_EPOCH)
        .unwrap()
        .as_nanos();
    let path = std::env::temp_dir().join(format!(
        "xtalforge-dataset-{}-{stamp}.jsonl",
        std::process::id()
    ));
    dataset.save(&path).unwrap();
    let loaded = Dataset::load(&path).unwrap();
    assert_eq!(loaded, dataset);
    std::fs::remove_file(&path).unwrap();
}

#[test]
fn the_version_constant_is_written() {
    let text = sample().to_json_lines().unwrap();
    assert!(text.contains(&format!("\"dataset_version\":{DATASET_VERSION}")));
}
