//! Integration tests for RCSR identification (Tasks 5.6.1, 5.6.4, 5.6.5).

#![allow(
    clippy::unwrap_used,
    clippy::expect_used,
    clippy::panic,
    clippy::indexing_slicing
)]

use xtalforge_analysis::QuotientGraph;
use xtalforge_analysis::rcsr::{RCSR_LICENSE_NOTE, RcsrIndex};

fn key(key: &str) -> QuotientGraph {
    QuotientGraph::from_canonical_key(key).unwrap()
}

#[test]
fn common_nets_are_identified() {
    let index = RcsrIndex::from_embedded().unwrap();
    for (name, canonical) in [
        ("pcu", "3:1:0,0,-1,0,0;0,0,0,-1,0;0,0,0,0,-1"),
        ("dia", "3:2:0,1,0,0,0;0,1,0,0,1;0,1,0,1,0;0,1,1,0,0"),
        ("sql", "2:1:0,0,-1,0;0,0,0,-1"),
        (
            "tbo",
            "3:14:0,1,0,0,0;0,2,0,0,0;0,3,0,0,0;1,4,0,0,0;1,5,0,0,0;1,6,0,0,0;2,4,1,0,0;2,7,0,0,0;2,8,0,0,0;3,4,0,1,0;3,9,0,0,0;3,10,0,0,0;5,11,0,0,0;5,12,0,0,0;6,11,0,0,1;6,12,1,-1,1;7,12,1,0,0;7,13,0,0,0;8,12,1,-1,1;8,13,0,0,1;9,11,0,0,1;9,13,0,0,1;10,11,0,1,0;10,13,-1,1,0",
        ),
        (
            "fcu",
            "3:1:0,0,-1,0,0;0,0,-1,0,1;0,0,-1,1,0;0,0,-1,1,1;0,0,0,-1,0;0,0,0,0,-1",
        ),
        (
            "sod",
            "3:6:0,1,0,0,0;0,2,0,0,0;0,3,0,0,0;0,4,0,0,0;1,2,1,0,0;1,3,0,1,0;1,5,0,0,0;2,4,0,0,1;2,5,0,-1,1;3,4,1,-1,1;3,5,0,-1,1;4,5,0,0,0",
        ),
        (
            "msf",
            "3:24:0,1,0,0,0;0,2,0,0,0;0,3,0,0,0;1,4,0,0,0;1,5,0,0,0;2,5,0,0,0;2,6,0,0,0;3,6,0,0,0;3,6,1,0,0;3,7,0,0,0;3,8,0,0,0;4,9,0,0,0;4,10,0,0,0;4,11,-1,0,0;4,11,0,0,0;4,12,0,0,0;4,13,0,0,0;5,11,-1,0,0;6,7,0,0,0;6,8,-1,0,0;6,14,0,0,0;6,15,0,0,0;7,8,0,0,0;7,14,0,0,0;7,14,1,0,0;7,15,0,0,0;7,16,0,0,0;8,14,1,0,0;8,15,0,0,0;8,15,1,0,0;8,17,0,0,0;9,10,-1,0,0;9,10,0,0,0;9,12,0,0,0;9,18,0,0,0;10,11,0,0,0;10,12,1,0,0;10,13,0,0,0;10,19,0,0,0;11,12,1,0,0;12,13,-1,0,0;12,13,0,0,0;12,20,0,0,0;13,21,0,0,0;14,22,0,0,0;15,23,0,0,0;16,21,0,1,0;16,22,1,0,0;17,18,0,0,1;17,23,1,0,0;18,19,0,0,0;19,23,1,0,-1;20,21,-1,0,0;20,22,0,-1,0",
        ),
        (
            "srs",
            "3:4:0,1,0,0,0;0,2,0,0,0;0,3,0,0,0;1,2,1,0,0;1,3,0,1,0;2,3,0,0,1",
        ),
    ] {
        let found = index.identify(&key(canonical)).unwrap();
        assert_eq!(found.name(), Some(name), "{name}: {found:?}");
        assert_eq!(found.verdict(), "named");
    }
}

#[test]
fn an_unknown_graph_reports_its_fingerprint() {
    let index = RcsrIndex::from_embedded().unwrap();
    let unknown = key("3:1:0,0,1,0,0");
    let found = index.identify(&unknown).unwrap();
    assert!(found.names.is_empty());
    assert_eq!(found.verdict(), "unknown");
    assert_eq!(found.fingerprint.periodicity, 1);
    assert!(!found.fingerprint.sequences.is_empty());
}

#[test]
fn the_full_index_validates() {
    let index = RcsrIndex::from_embedded().unwrap();
    let count = index.validate().unwrap();
    assert_eq!(count, index.len());
    assert!(count > 2_900, "only {count} catalogue entries");
    assert!(index.names().contains(&"pcu"));
    assert!(index.names().contains(&"msf"));
}

#[test]
fn the_licence_note_names_the_bundled_snapshot() {
    assert!(!RCSR_LICENSE_NOTE.trim().is_empty());
    assert!(RCSR_LICENSE_NOTE.contains("RCSR"));
    assert!(RCSR_LICENSE_NOTE.contains("2019-06-01"));
}
