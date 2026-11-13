//! Integration tests for Blatov interpenetration classes (Task 14.4.2).

#![allow(
    clippy::unwrap_used,
    clippy::expect_used,
    clippy::panic,
    clippy::indexing_slicing
)]

use xtalforge_analysis::interpenetration::{
    BlatovClass, classify, classify_with_symmetry, cycle_voltages, saturation_index,
};
use xtalforge_analysis::{NetEdge, QuotientGraph};
use xtalforge_core::SymOp;

fn graph(nodes: Vec<[f64; 3]>, edges: Vec<(usize, usize, [i32; 3])>) -> QuotientGraph {
    QuotientGraph {
        nodes,
        edges: edges
            .into_iter()
            .map(|(a, b, shift)| NetEdge { a, b, shift })
            .collect(),
        labels: Vec::new(),
    }
}

fn inversion() -> SymOp {
    SymOp::new([[-1, 0, 0], [0, -1, 0], [0, 0, -1]], [0.0, 0.0, 0.0]).unwrap()
}

#[test]
fn a_doubled_cell_is_class_ia_from_its_voltages() {
    let doubled = QuotientGraph::from_canonical_key("3:1:0,0,0,0,2").unwrap();
    assert_eq!(saturation_index(&cycle_voltages(&doubled)), 2);
    let report = classify(&doubled).unwrap();
    assert_eq!(report.class, BlatovClass::ClassIa);
    assert_eq!(report.class_number(), Some(1));
    assert_eq!(report.components, 1);
    assert_eq!(report.fold, 2);
    assert!(!report.needs_symmetry());
}

#[test]
fn translation_related_components_are_class_ia() {
    let translated = graph(
        vec![
            [0.0, 0.0, 0.0],
            [0.1, 0.0, 0.0],
            [0.3, 0.0, 0.0],
            [0.4, 0.0, 0.0],
        ],
        vec![
            (0, 1, [0, 0, 0]),
            (0, 0, [0, 0, 1]),
            (2, 3, [0, 0, 0]),
            (2, 2, [0, 0, 1]),
        ],
    );
    let report = classify(&translated).unwrap();
    assert_eq!(report.class, BlatovClass::ClassIa);
    assert_eq!(report.components, 2);
    assert_eq!(report.frameworks, 2);
    assert_eq!(report.fold, 2);
}

fn inversion_pair() -> QuotientGraph {
    graph(
        vec![
            [0.1, 0.0, 0.0],
            [0.2, 0.0, 0.0],
            [0.35, 0.0, 0.0],
            [0.9, 0.0, 0.0],
            [0.8, 0.0, 0.0],
            [0.65, 0.0, 0.0],
        ],
        vec![
            (0, 1, [0, 0, 0]),
            (1, 2, [0, 0, 0]),
            (0, 0, [0, 0, 1]),
            (3, 4, [0, 0, 0]),
            (4, 5, [0, 0, 0]),
            (3, 3, [0, 0, 1]),
        ],
    )
}

#[test]
fn inversion_related_components_need_symmetry_and_are_class_ii() {
    let pair = inversion_pair();
    let plain = classify(&pair).unwrap();
    assert_eq!(plain.class, BlatovClass::Undetermined);
    assert_eq!(plain.class_number(), None);
    assert!(plain.needs_symmetry());
    assert_eq!(plain.components, 2);

    let resolved = classify_with_symmetry(&pair, &[inversion()], 1e-4).unwrap();
    assert_eq!(resolved.class, BlatovClass::ClassII);
    assert_eq!(resolved.class_number(), Some(2));
}

#[test]
fn unrelated_components_are_class_iii() {
    let unrelated = graph(
        vec![
            [0.1, 0.0, 0.0],
            [0.2, 0.0, 0.0],
            [0.35, 0.0, 0.0],
            [0.9, 0.0, 0.0],
            [0.8, 0.0, 0.0],
            [0.7, 0.0, 0.0],
        ],
        vec![
            (0, 1, [0, 0, 0]),
            (1, 2, [0, 0, 0]),
            (0, 0, [0, 0, 1]),
            (3, 4, [0, 0, 0]),
            (4, 5, [0, 0, 0]),
            (3, 3, [0, 0, 1]),
        ],
    );
    let resolved = classify_with_symmetry(&unrelated, &[inversion()], 1e-4).unwrap();
    assert_eq!(resolved.class, BlatovClass::ClassIII);
    assert_eq!(resolved.class_number(), Some(3));
}

#[test]
fn a_single_net_is_not_interpenetrated() {
    let single =
        QuotientGraph::from_canonical_key("3:2:0,1,0,0,0;0,1,0,0,1;0,1,0,1,0;0,1,1,0,0").unwrap();
    let report = classify(&single).unwrap();
    assert_eq!(report.class, BlatovClass::None);
    assert_eq!(report.class_number(), Some(0));
    assert_eq!(report.components, 1);
    assert_eq!(report.fold, 1);
}

#[test]
fn interpenetration_errors_and_index_arithmetic() {
    let empty = QuotientGraph {
        nodes: Vec::new(),
        edges: Vec::new(),
        labels: Vec::new(),
    };
    assert_eq!(classify(&empty).unwrap_err().code(), "invalid_input");
    let single = QuotientGraph::from_canonical_key("3:1:0,0,0,0,1").unwrap();
    assert!(classify_with_symmetry(&single, &[], 0.0).is_err());

    assert_eq!(saturation_index(&[]), 1);
    assert_eq!(saturation_index(&[[0, 0, 2]]), 2);
    assert_eq!(saturation_index(&[[2, 0, 0], [0, 2, 0]]), 4);
    assert_eq!(saturation_index(&[[1, 0, 0], [0, 1, 0], [0, 0, 1]]), 1);
}
