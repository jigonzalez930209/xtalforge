//! Integration tests for dataset deduplication (Task 16.3.2).

#![allow(
    clippy::unwrap_used,
    clippy::expect_used,
    clippy::panic,
    clippy::indexing_slicing
)]

use xtalforge_analysis::{Dataset, DatasetRecord, DedupOptions, deduplicate};
use xtalforge_core::{Lattice, Site, Structure, StructureBuilder};
use xtalforge_symmetry::{MatcherOptions, SpeciesComparator};

fn halide(cation: &str) -> Structure {
    let positions = [
        (cation, [0.0, 0.0, 0.0]),
        (cation, [0.0, 0.5, 0.5]),
        (cation, [0.5, 0.0, 0.5]),
        (cation, [0.5, 0.5, 0.0]),
        ("Cl", [0.5, 0.5, 0.5]),
        ("Cl", [0.5, 0.0, 0.0]),
        ("Cl", [0.0, 0.5, 0.0]),
        ("Cl", [0.0, 0.0, 0.5]),
    ];
    let mut builder = StructureBuilder::new("halide", Lattice::cubic(5.64).unwrap());
    for (index, (element, position)) in positions.iter().enumerate() {
        builder = builder.site(Site::new(format!("{element}{index}"), element, *position).unwrap());
    }
    builder.build().unwrap()
}

fn reversed_nacl() -> Structure {
    let mut sites: Vec<Site> = halide("Na").sites().to_vec();
    sites.reverse();
    let mut builder = StructureBuilder::new("NaCl-reversed", Lattice::cubic(5.64).unwrap());
    for site in sites {
        builder = builder.site(site);
    }
    builder.build().unwrap()
}

fn cscl() -> Structure {
    StructureBuilder::new("CsCl", Lattice::cubic(4.2).unwrap())
        .site(Site::new("Cs1", "Cs", [0.0, 0.0, 0.0]).unwrap())
        .site(Site::new("Cl1", "Cl", [0.5, 0.5, 0.5]).unwrap())
        .build()
        .unwrap()
}

fn primitive_nacl() -> Structure {
    let a = 5.64;
    let lattice = Lattice::new([
        [0.0, a / 2.0, a / 2.0],
        [a / 2.0, 0.0, a / 2.0],
        [a / 2.0, a / 2.0, 0.0],
    ])
    .unwrap();
    StructureBuilder::new("NaCl-primitive", lattice)
        .site(Site::new("Na1", "Na", [0.0, 0.0, 0.0]).unwrap())
        .site(Site::new("Cl1", "Cl", [0.5, 0.5, 0.5]).unwrap())
        .build()
        .unwrap()
}

fn sample() -> Dataset {
    let records = vec![
        DatasetRecord::new("nacl-a", halide("Na"))
            .unwrap()
            .with_property("energy", -3.0)
            .unwrap(),
        DatasetRecord::new("cscl", cscl()).unwrap(),
        DatasetRecord::new("nacl-b", reversed_nacl())
            .unwrap()
            .with_property("energy", -3.0)
            .unwrap(),
        DatasetRecord::new("nacl-primitive", primitive_nacl()).unwrap(),
    ];
    Dataset::from_records("halides", records).unwrap()
}

#[test]
fn duplicates_are_grouped_in_input_order() {
    let dataset = sample();
    let result = deduplicate(&dataset, &DedupOptions::default()).unwrap();

    assert_eq!(result.report.n_input, 4);
    assert_eq!(result.report.n_kept(), 2);
    assert_eq!(result.report.n_removed(), 2);
    assert_eq!(result.report.kept, [0, 1]);
    assert_eq!(result.report.removed(), [2, 3]);
    assert!(!result.report.is_duplicate_free());

    assert_eq!(result.report.groups.len(), 2);
    assert_eq!(result.report.groups[0].kept, 0);
    assert_eq!(result.report.groups[0].removed, [2, 3]);
    assert_eq!(result.report.groups[1].kept, 1);
    assert!(result.report.groups[1].removed.is_empty());

    let ids: Vec<&str> = result
        .dataset
        .iter()
        .map(|record| record.id.as_str())
        .collect();
    assert_eq!(ids, ["nacl-a", "cscl"]);
    assert_eq!(result.dataset.metadata().name, "halides");
}

#[test]
fn deduplication_is_deterministic_and_apply_matches() {
    let dataset = sample();
    let first = deduplicate(&dataset, &DedupOptions::default()).unwrap();
    let second = deduplicate(&dataset, &DedupOptions::default()).unwrap();
    assert_eq!(first, second);
    assert_eq!(first.report.apply(&dataset).unwrap(), first.dataset);
}

#[test]
fn a_duplicate_free_dataset_is_reported_as_such() {
    let records = vec![
        DatasetRecord::new("nacl", halide("Na")).unwrap(),
        DatasetRecord::new("cscl", cscl()).unwrap(),
    ];
    let dataset = Dataset::from_records("distinct", records).unwrap();
    let result = deduplicate(&dataset, &DedupOptions::default()).unwrap();
    assert!(result.report.is_duplicate_free());
    assert_eq!(result.dataset, dataset);
}

#[test]
fn an_empty_dataset_yields_an_empty_report() {
    let dataset = Dataset::new("empty");
    let result = deduplicate(&dataset, &DedupOptions::default()).unwrap();
    assert!(result.dataset.is_empty());
    assert_eq!(result.report.n_input, 0);
    assert_eq!(result.report.n_kept(), 0);
    assert!(result.report.is_duplicate_free());
}

#[test]
fn disabling_supercells_keeps_the_expanded_record() {
    let dataset = sample();
    let options = DedupOptions {
        matcher: MatcherOptions {
            allow_supercell: false,
            ..MatcherOptions::default()
        },
        ..DedupOptions::default()
    };
    let result = deduplicate(&dataset, &options).unwrap();
    assert_eq!(result.report.n_kept(), 3);
    assert_eq!(result.report.removed(), [2]);
}

#[test]
fn ignoring_species_merges_geometric_duplicates() {
    let records = vec![
        DatasetRecord::new("nacl", halide("Na")).unwrap(),
        DatasetRecord::new("kcl", halide("K")).unwrap(),
    ];
    let dataset = Dataset::from_records("geometry", records).unwrap();
    let exact = deduplicate(&dataset, &DedupOptions::default()).unwrap();
    assert_eq!(exact.report.n_kept(), 2);

    let options = DedupOptions {
        matcher: MatcherOptions {
            comparator: SpeciesComparator::Ignore,
            ..MatcherOptions::default()
        },
        ..DedupOptions::default()
    };
    let ignored = deduplicate(&dataset, &options).unwrap();
    assert_eq!(ignored.report.n_kept(), 1);
    assert_eq!(ignored.report.removed(), [1]);
}

#[test]
fn the_optional_similarity_gate_is_applied() {
    let dataset = sample();
    let strict = DedupOptions {
        min_fingerprint_similarity: Some(0.999),
        ..DedupOptions::default()
    };
    let result = deduplicate(&dataset, &strict).unwrap();
    // The exact duplicate still matches; the primitive cell has a lower
    // fingerprint similarity and is filtered before the matcher.
    assert_eq!(result.report.removed(), [2]);

    let disabled = DedupOptions {
        min_fingerprint_similarity: Some(0.0),
        ..DedupOptions::default()
    };
    assert_eq!(
        deduplicate(&dataset, &disabled).unwrap().report.removed(),
        [2, 3]
    );
}

#[test]
fn invalid_options_are_rejected() {
    let dataset = sample();
    for threshold in [-0.1, 1.1, f64::NAN] {
        let options = DedupOptions {
            min_fingerprint_similarity: Some(threshold),
            ..DedupOptions::default()
        };
        assert!(deduplicate(&dataset, &options).is_err());
    }
    let options = DedupOptions {
        fingerprint: xtalforge_core::RdfOptions {
            r_max: 0.0,
            ..xtalforge_core::RdfOptions::default()
        },
        ..DedupOptions::default()
    };
    assert!(deduplicate(&dataset, &options).is_err());
}

#[test]
fn an_empty_structure_is_reported_with_its_index() {
    let empty = StructureBuilder::new("empty", Lattice::cubic(5.0).unwrap())
        .build()
        .unwrap();
    let records = vec![
        DatasetRecord::new("nacl", halide("Na")).unwrap(),
        DatasetRecord::new("empty", empty).unwrap(),
    ];
    let dataset = Dataset::from_records("broken", records).unwrap();
    let error = deduplicate(&dataset, &DedupOptions::default()).unwrap_err();
    let text = error.to_string();
    assert!(text.contains("record 1"), "{text}");
    assert!(text.contains("empty"), "{text}");
}

#[test]
fn the_report_serializes_and_rejects_bad_indices() {
    let dataset = sample();
    let result = deduplicate(&dataset, &DedupOptions::default()).unwrap();
    let json = result.report.to_json().unwrap();
    assert!(json.contains("\"n_input\": 4"));
    let parsed: xtalforge_analysis::DedupReport = serde_json::from_str(&json).unwrap();
    assert_eq!(parsed, result.report);

    let mut broken = parsed;
    broken.kept.push(99);
    assert!(broken.apply(&dataset).is_err());

    let mut repeated = result.report;
    repeated.kept.push(0);
    assert!(repeated.apply(&dataset).is_err());
}
