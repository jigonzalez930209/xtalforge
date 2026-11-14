//! Checkpoint and workflow tests (Task 20.2.3).

#![allow(
    clippy::unwrap_used,
    clippy::expect_used,
    clippy::panic,
    clippy::indexing_slicing
)]

use serde_json::json;
use xtalforge_generation::{
    CHECKPOINT_VERSION, SearchCheckpoint, SearchKind, run_random_chunk, search_workflow_document,
    search_workflow_json,
};

#[test]
fn checkpoint_round_trips() {
    let checkpoint = SearchCheckpoint::new(SearchKind::Random, json!({ "seed": 1 }), 4);
    let text = checkpoint.to_json().expect("json");
    let parsed = SearchCheckpoint::from_json(&text).expect("parse");
    assert_eq!(parsed, checkpoint);
    assert_eq!(parsed.checkpoint_version, CHECKPOINT_VERSION);
    assert!(!parsed.is_complete());
}

#[test]
fn unknown_checkpoint_version_is_rejected() {
    let text = r#"{"checkpoint_version":99,"kind":"random","options":{},"completed":0,"total":0,"hits":[]}"#;
    assert!(SearchCheckpoint::from_json(text).is_err());
    assert!(SearchCheckpoint::from_json("not json").is_err());
}

#[test]
fn random_chunk_resumes_from_the_checkpoint() {
    let options = json!({
        "composition": [["Na", 4], ["Cl", 4]],
        "space_groups": [225],
        "n_structures": 2,
        "seed": 40,
        "volume_per_atom": 22.0,
        "min_distance": 1.2,
        "max_attempts": 32,
        "relax": false,
        "relax_cell": false,
        "max_steps": 10,
        "deduplicate": false,
        "matcher": {
            "ltol": 0.2,
            "stol": 0.3,
            "angle_tol": 5.0,
            "comparator": "exact",
            "allow_supercell": true
        },
        "minimizer": "fire"
    });
    let mut checkpoint = SearchCheckpoint::new(SearchKind::Random, options, 2);
    let first = run_random_chunk(&mut checkpoint, 1).expect("chunk 1");
    assert!(!first.complete);
    assert_eq!(first.checkpoint.completed, 1);
    assert_eq!(first.hits.len(), 1);
    let second = run_random_chunk(&mut checkpoint, 1).expect("chunk 2");
    assert!(second.complete);
    assert_eq!(second.checkpoint.completed, 2);
    assert_eq!(second.hits.len(), 2);
}

#[test]
fn workflow_document_has_read_search_write() {
    let document = search_workflow_document(
        "search",
        SearchKind::Random,
        "seed.cif",
        json!({ "seed": 1 }),
        Some("search.checkpoint.json"),
        "results.json",
    );
    assert_eq!(document["workflow_version"], 1);
    assert_eq!(document["steps"].as_array().expect("steps").len(), 3);
    assert_eq!(document["steps"][0]["kind"], "read");
    assert_eq!(document["steps"][1]["verb"], "search");
    assert_eq!(document["steps"][1]["arguments"]["kind"], "random");
    assert_eq!(document["steps"][2]["kind"], "write");
    let text = search_workflow_json(
        "search",
        SearchKind::Evolutionary,
        "seed.cif",
        json!({}),
        None,
        "results.json",
    )
    .expect("json");
    assert!(text.contains("\"verb\": \"search\""));
}
