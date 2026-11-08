//! Integration tests for rods, layers, merging and invented linkers (Tasks 6.14-6.16).

#![allow(
    clippy::unwrap_used,
    clippy::expect_used,
    clippy::panic,
    clippy::indexing_slicing
)]

use xtalforge_mof::{
    build_layered_framework, build_rod_framework, decorate_super_tetrahedron, invent_linker,
    linkers_by_span, merge_frameworks, partition_pore_space,
};

#[test]
fn linkers_have_spans_and_can_be_invented() {
    let linkers = linkers_by_span(0.0, 30.0).unwrap();
    assert!(linkers.iter().any(|linker| linker.name == "benzene"));
    assert!(linkers.windows(2).all(|pair| pair[0].span <= pair[1].span));
    let invented = invent_linker(7.0, true).unwrap();
    assert_eq!(invented.connection_count(), 2);
    assert!((invented.connection_points[1][0] - 7.0).abs() < 1e-12);
    assert!(invent_linker(0.5, false).is_err());
}

#[test]
fn rods_layers_and_merges_build() {
    let rod = build_rod_framework(2).unwrap();
    assert_eq!(rod.formula().matches("Zn").count(), 1);
    assert!(rod.sites().len() >= 10);
    let layered = build_layered_framework(2).unwrap();
    assert_eq!(layered.sites().len(), 6);
    let merged = merge_frameworks(&rod, &layered, [0.0, 0.0, 0.5]).unwrap();
    assert_eq!(
        merged.sites().len(),
        rod.sites().len() + layered.sites().len()
    );
    assert!(merge_frameworks(&rod, &layered, [f64::NAN, 0.0, 0.0]).is_err());
    assert!(build_rod_framework(0).is_err());
}

#[test]
fn decoration_and_pore_partition() {
    let rod = build_rod_framework(1).unwrap();
    let decorated = decorate_super_tetrahedron(&rod, 0).unwrap();
    assert_eq!(decorated.sites().len(), rod.sites().len() + 4);
    assert!(decorate_super_tetrahedron(&rod, 99).is_err());
    let volumes = partition_pore_space(std::slice::from_ref(&rod), 0.6, 1.2).unwrap();
    assert_eq!(volumes.len(), 1);
    assert!(volumes[0] >= 0.0);
    assert!(partition_pore_space(&[], 0.6, 1.2).is_err());
}
