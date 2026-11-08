//! Integration tests for the curated MOF builder (Phase 6 Tasks 6.11-6.16).

#![allow(
    clippy::unwrap_used,
    clippy::expect_used,
    clippy::panic,
    clippy::indexing_slicing
)]

use std::collections::BTreeSet;

use xtalforge_mof::{
    Block, BlockKind, BuiltMof, Error, MofRequest, Topology, blocks, build_mof, build_request_json,
    find_block, find_topology, node_oriented, topologies,
};

const EXPECTED_DEGREES: [(&str, usize); 6] = [
    ("pcu", 6),
    ("dia", 4),
    ("srs", 3),
    ("sql", 4),
    ("hcb", 3),
    ("bcu", 8),
];

#[test]
fn every_topology_loads_and_validates() {
    let all = topologies().unwrap();
    assert!(
        (6..=8).contains(&all.len()),
        "unexpected table size {}",
        all.len()
    );
    let names: BTreeSet<&str> = all.iter().map(|topology| topology.name.as_str()).collect();
    assert_eq!(names.len(), all.len(), "topology names must be unique");
    for (name, degree) in EXPECTED_DEGREES {
        let topology = find_topology(name)
            .unwrap()
            .unwrap_or_else(|| panic!("missing {name}"));
        topology.validate().unwrap();
        assert!(!topology.nodes.is_empty(), "{name} has no nodes");
        assert!(!topology.edges.is_empty(), "{name} has no edges");
        assert!(!topology.source.trim().is_empty(), "{name} has no source");
        assert!(topology.cell.is_some(), "{name} should carry a cell");
        for vertex in 0..topology.nodes.len() {
            assert_eq!(
                topology.degree(vertex),
                degree,
                "{name} vertex {vertex} degree"
            );
        }
    }
    assert!(
        find_topology("PCU").unwrap().is_some(),
        "lookup is case-insensitive"
    );
    assert!(find_topology("missing").unwrap().is_none());

    let default_cell = Topology {
        name: "bare".to_owned(),
        nodes: vec![[0.0, 0.0, 0.0]],
        edges: Vec::new(),
        cell: None,
        source: "test".to_owned(),
    };
    assert_eq!(
        default_cell.lattice().unwrap().lengths(),
        [20.0, 20.0, 20.0]
    );
}

#[test]
fn topology_validation_rejects_bad_edges() {
    let base = Topology {
        name: "bad".to_owned(),
        nodes: vec![[0.0, 0.0, 0.0]],
        edges: vec![[0, 2]],
        cell: Some([10.0, 10.0, 10.0, 90.0, 90.0, 90.0]),
        source: "test".to_owned(),
    };
    assert!(matches!(base.validate(), Err(Error::InvalidInput(_))));
    let duplicate = Topology {
        edges: vec![[0, 1], [1, 0]],
        nodes: vec![[0.0, 0.0, 0.0], [0.5, 0.0, 0.0]],
        ..base
    };
    assert!(matches!(duplicate.validate(), Err(Error::InvalidInput(_))));
}

#[test]
fn every_block_loads_and_validates() {
    let all = blocks().unwrap();
    assert_eq!(all.len(), 8, "unexpected block table size");
    let names: BTreeSet<&str> = all.iter().map(|block| block.name.as_str()).collect();
    assert_eq!(names.len(), all.len(), "block names must be unique");
    for block in &all {
        block.validate().unwrap();
        assert_eq!(block.points.len(), block.elements.len());
        assert!(block.atom_count() <= 60, "{} is too large", block.name);
        assert!(block.connection_count() >= 1);
    }
    let expected: [(&str, BlockKind, usize, usize); 8] = [
        ("Zn4O", BlockKind::Node, 6, 5),
        ("Cu2", BlockKind::Node, 4, 10),
        ("Cr3O", BlockKind::Node, 6, 10),
        ("benzene_tri", BlockKind::Node, 3, 18),
        ("benzene", BlockKind::Linker, 2, 16),
        ("biphenyl", BlockKind::Linker, 2, 26),
        ("stilbene", BlockKind::Linker, 2, 30),
        ("triphenylene", BlockKind::Linker, 3, 30),
    ];
    for (name, kind, connections, atoms) in expected {
        let block = find_block(name)
            .unwrap()
            .unwrap_or_else(|| panic!("missing {name}"));
        assert_eq!(block.kind, kind, "{name} kind");
        assert_eq!(block.connection_count(), connections, "{name} connections");
        assert_eq!(block.atom_count(), atoms, "{name} atoms");
    }
    assert!(find_block("BENZENE").unwrap().is_some());
    assert!(find_block("missing").unwrap().is_none());
}

#[test]
fn block_validation_rejects_bad_data() {
    let node = find_block("Zn4O").unwrap().unwrap();
    let short = Block {
        elements: node.elements[..2].to_vec(),
        ..node.clone()
    };
    assert!(matches!(short.validate(), Err(Error::InvalidInput(_))));
    let no_connections = Block {
        connection_points: Vec::new(),
        ..node
    };
    assert!(matches!(
        no_connections.validate(),
        Err(Error::InvalidInput(_))
    ));
}

#[test]
fn pcu_zn4o_benzene_builds_finite_atoms() {
    let built = build_mof(&MofRequest::new("pcu", "Zn4O", "benzene")).unwrap();
    assert!(
        built.structure.sites().len() > 20,
        "only {} atoms",
        built.structure.sites().len()
    );
    assert_eq!(built.topology, "pcu");
    assert_eq!(built.node, "Zn4O");
    assert_eq!(built.linker, "benzene");
    assert!(built.structure.is_p1());
    assert!((built.structure.lattice().lengths()[0] - 36.0).abs() < 1e-9);
    assert_eq!(built.structure.topology_bonds().len(), 81);
    for site in built.structure.sites() {
        assert!(site.frac().iter().all(|value| value.is_finite()));
        assert!(site.frac().iter().all(|value| (0.0..1.0).contains(value)));
    }
    assert_unique_positions(&built);
    assert!(
        built.warnings.is_empty(),
        "unexpected warnings: {:?}",
        built.warnings
    );
}

#[test]
fn sql_cu2_benzene_builds_without_warnings() {
    let built = build_mof(&MofRequest::new("sql", "Cu2", "benzene")).unwrap();
    assert!(built.structure.sites().len() > 100);
    assert_eq!(built.structure.topology_bonds().len(), 18);
    assert!(
        built.warnings.is_empty(),
        "unexpected warnings: {:?}",
        built.warnings
    );
    assert_unique_positions(&built);
}

#[test]
fn hcb_benzene_tri_triphenylene_builds() {
    let built = build_mof(&MofRequest::new("hcb", "benzene_tri", "triphenylene")).unwrap();
    assert!(built.structure.sites().len() > 20);
    assert_eq!(built.structure.topology_bonds().len(), 12);
    assert!(
        built.warnings.is_empty(),
        "unexpected warnings: {:?}",
        built.warnings
    );
}

#[test]
fn degree_mismatch_is_a_warning_not_an_error() {
    let built = build_mof(&MofRequest::new("pcu", "Cu2", "benzene")).unwrap();
    assert!(
        built
            .warnings
            .iter()
            .any(|warning| warning.contains("degree")),
        "expected a degree warning, got {:?}",
        built.warnings
    );
    assert!(!built.structure.sites().is_empty());
}

#[test]
fn unknown_names_are_typed_errors() {
    let topology_error = build_mof(&MofRequest::new("nope", "Zn4O", "benzene")).unwrap_err();
    assert!(matches!(topology_error, Error::UnknownTopology(_)));
    assert_eq!(topology_error.code(), "unknown_topology");

    let node_error = build_mof(&MofRequest::new("pcu", "nope", "benzene")).unwrap_err();
    assert!(matches!(node_error, Error::UnknownBlock(_)));
    assert_eq!(node_error.code(), "unknown_block");

    let kind_error = build_mof(&MofRequest::new("pcu", "benzene", "benzene")).unwrap_err();
    assert!(matches!(kind_error, Error::InvalidInput(_)));
    assert_eq!(kind_error.code(), "invalid_input");

    let json_error = build_request_json("{\"topology\": 7}").unwrap_err();
    assert!(matches!(json_error, Error::Json(_)));
    assert_eq!(json_error.code(), "json");
}

#[test]
fn request_json_round_trip_builds_the_same_framework() {
    let request = MofRequest::new("sql", "Cu2", "benzene");
    let text = serde_json::to_string(&request).unwrap();
    let from_json = build_request_json(&text).unwrap();
    let direct = build_mof(&request).unwrap();
    assert_eq!(from_json, direct);

    let serialized = from_json.to_json().unwrap();
    let decoded: BuiltMof = serde_json::from_str(&serialized).unwrap();
    assert_eq!(decoded, from_json);
}

#[test]
fn node_oriented_maps_the_reference_direction() {
    let block = find_block("benzene_tri").unwrap().unwrap();
    let reference = unit(block.connection_points[0]);
    let rotation = node_oriented(&block, [0.0, 0.0, 1.0]).unwrap();
    let mapped = apply(rotation, reference);
    assert!((mapped[0]).abs() < 1e-9, "{mapped:?}");
    assert!((mapped[1]).abs() < 1e-9, "{mapped:?}");
    assert!((mapped[2] - 1.0).abs() < 1e-9, "{mapped:?}");

    let opposite = node_oriented(&block, [-1.0, 0.0, 0.0]).unwrap();
    let mapped = apply(opposite, reference);
    assert!((mapped[0] + 1.0).abs() < 1e-9, "{mapped:?}");
    assert!((mapped[1]).abs() < 1e-9, "{mapped:?}");
    assert!((mapped[2]).abs() < 1e-9, "{mapped:?}");

    assert!(matches!(
        node_oriented(&block, [0.0, 0.0, 0.0]),
        Err(Error::InvalidInput(_))
    ));
    let empty = Block {
        connection_points: Vec::new(),
        ..block
    };
    assert!(matches!(
        node_oriented(&empty, [1.0, 0.0, 0.0]),
        Err(Error::InvalidInput(_))
    ));
}

fn assert_unique_positions(built: &BuiltMof) {
    let mut seen = BTreeSet::new();
    for site in built.structure.sites() {
        let [x, y, z] = site.frac();
        let key = (
            (snap(x) * 1e4).round() as i64,
            (snap(y) * 1e4).round() as i64,
            (snap(z) * 1e4).round() as i64,
        );
        assert!(seen.insert(key), "duplicate atom at {:?}", site.frac());
    }
}

fn snap(value: f64) -> f64 {
    if 1.0 - value < 1e-4 { 0.0 } else { value }
}

fn unit(vector: [f64; 3]) -> [f64; 3] {
    let [x, y, z] = vector;
    let length = (x * x + y * y + z * z).sqrt();
    [x / length, y / length, z / length]
}

fn apply(matrix: [[f64; 3]; 3], vector: [f64; 3]) -> [f64; 3] {
    let [x, y, z] = vector;
    matrix.map(|[a, b, c]| a * x + b * y + c * z)
}
