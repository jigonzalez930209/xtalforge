//! Golden bond-count comparison (Task 1.11.5).

#![allow(
    clippy::unwrap_used,
    clippy::expect_used,
    clippy::panic,
    clippy::indexing_slicing
)]

use xtalforge_io::read_cif;

#[test]
fn mof5_bond_count_matches_the_golden_data() {
    let structure = read_cif(include_str!("../../../tests/corpus/cod/MOF-5.cif")).unwrap();
    let graph = structure.bond_graph().unwrap();
    let golden: serde_json::Value =
        serde_json::from_str(include_str!("../../../tests/golden/phase-1/cod/MOF-5.json")).unwrap();
    let expected = golden["bonds"].as_array().unwrap().len();
    assert_eq!(graph.bonds.len(), expected);
}

#[test]
fn silicon_bond_count_matches_the_golden_data() {
    let structure = read_cif(include_str!("../../../tests/corpus/simple/Si.cif")).unwrap();
    let graph = structure.bond_graph().unwrap();
    let golden: serde_json::Value =
        serde_json::from_str(include_str!("../../../tests/golden/phase-1/simple/Si.json")).unwrap();
    let expected = golden["bonds"].as_array().unwrap().len();
    assert_eq!(graph.bonds.len(), expected);
}
