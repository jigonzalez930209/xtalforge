//! Integration tests for the CGD net format (Task 3.4.1).

#![allow(
    clippy::unwrap_used,
    clippy::expect_used,
    clippy::panic,
    clippy::indexing_slicing
)]

use xtalforge_io::cgd::{CgdEntry, CgdNode, read_cgd, write_cgd};

fn net() -> CgdEntry {
    CgdEntry {
        name: "test-net".to_owned(),
        group: "P1".to_owned(),
        nodes: vec![
            CgdNode {
                label: "A".to_owned(),
                coordination: 4,
                frac: vec![0.0, 0.0, 0.0],
            },
            CgdNode {
                label: "B".to_owned(),
                coordination: 4,
                frac: vec![0.5, 0.5, 0.5],
            },
        ],
        edges: vec![
            (vec![0.0, 0.0, 0.0], vec![0.5, 0.5, 0.5]),
            (vec![0.0, 0.0, 0.0], vec![0.0, 0.0, 1.0]),
        ],
        cell: Some(vec![4.0, 4.0, 4.0, 90.0, 90.0, 90.0]),
        line: 1,
    }
}

#[test]
fn cgd_round_trip() {
    let text = write_cgd(&[net()]).unwrap();
    assert!(text.contains("CRYSTAL"));
    let file = read_cgd(&text);
    assert!(file.problems.is_empty());
    assert_eq!(file.entries.len(), 1);
    assert_eq!(file.entries[0].name, "test-net");
    assert_eq!(file.entries[0].nodes.len(), 2);
    assert_eq!(file.entries[0].edges.len(), 2);
}

#[test]
fn cgd_reports_broken_blocks() {
    let file = read_cgd(
        "CRYSTAL\nNAME ok\nNODE A 2 0 0 0\nEDGE 0 0 0 0 0 1\nEND\nCRYSTAL\nNAME broken\nEND\n",
    );
    assert_eq!(file.entries.len(), 1);
    assert_eq!(file.problems.len(), 1);
    assert!(file.problems[0].contains("broken"));
    assert!(
        write_cgd(&[CgdEntry {
            name: String::new(),
            ..net()
        }])
        .is_err()
    );
}

#[test]
fn both_edge_dialects_and_bare_keywords_are_read() {
    let text = "CRYSTAL\nNAME pair\nGROUP P1\nNODE 1 2 0 0 0\nNODE 2 2 0.5 0 0\nEDGE 0 0 0 0.5 0 0\nEDGE 0.5 0 0 0 0 1\nEND\nCRYSTAL\nNAME neighbour\nGROUP P1\nNODE 1 2 0 0 0\nEDGE 1 0.5 0 0\nEND\n";
    let file = read_cgd(text);
    assert!(file.problems.is_empty(), "{:?}", file.problems);
    assert_eq!(file.entries.len(), 2);
    assert_eq!(
        file.entries[0].edges[1],
        (vec![0.5, 0.0, 0.0], vec![0.0, 0.0, 1.0])
    );
    assert_eq!(
        file.entries[1].edges[0],
        (vec![0.0, 0.0, 0.0], vec![0.5, 0.0, 0.0])
    );
}

#[test]
fn the_bundled_rcsr_snapshot_reads() {
    use std::io::Read;

    let path = concat!(
        env!("CARGO_MANIFEST_DIR"),
        "/../../data/rcsr/RCSRnets-2019-06-01.cgd.gz"
    );
    let mut text = String::new();
    flate2::read::GzDecoder::new(std::fs::File::open(path).unwrap())
        .read_to_string(&mut text)
        .unwrap();
    let file = read_cgd(&text);
    assert_eq!(file.entries.len(), 2_931);
    assert_eq!(file.problems.len(), 1);
    assert!(file.problems[0].contains("moo-a"));
}
