//! Integration tests for the optional full PORMAKE bundle (Task 17.1.3).

#![cfg(feature = "mof-data")]
#![allow(clippy::unwrap_used, clippy::expect_used, clippy::panic)]

use xtalforge_mof::{
    MofRequest, blocks, build_mof, curated_blocks, curated_topologies, find_block, find_topology,
    topologies,
};

#[test]
fn the_full_bundle_loads_lazily_and_keeps_the_curated_subset() {
    let full = topologies().unwrap();
    assert!(full.len() > 2000, "only {} topologies", full.len());
    assert!(find_topology("tbo").unwrap().is_some());
    assert!(find_topology("rht").unwrap().is_some());
    assert_eq!(curated_topologies().unwrap().len(), 6);

    let all_blocks = blocks().unwrap();
    assert!(all_blocks.len() >= 867, "only {} blocks", all_blocks.len());
    assert_eq!(curated_blocks().unwrap().len(), 8);
    let node = find_block("N59").unwrap().expect("N59 node");
    let linker = find_block("E32").unwrap().expect("E32 linker");
    assert_ne!(node.connection_count(), 0);
    assert!(linker.connection_count() >= 2);
}

#[test]
fn a_full_bundle_net_builds_a_finite_framework() {
    let built = build_mof(&MofRequest::new("tbo", "N59", "E32")).unwrap();
    assert!(!built.structure.sites().is_empty());
    assert!(built.structure.lattice().volume() > 0.0);
}

#[test]
fn unknown_full_bundle_names_still_error() {
    assert!(find_topology("not-a-real-net").unwrap().is_none());
    assert!(find_block("N000000").unwrap().is_none());
}
