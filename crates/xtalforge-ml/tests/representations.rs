//! Graph, split and import tests (Tasks 21.1.1 and 21.2).

#![allow(
    clippy::unwrap_used,
    clippy::expect_used,
    clippy::panic,
    clippy::indexing_slicing
)]

use xtalforge_analysis::dataset::{Dataset, DatasetRecord};
use xtalforge_core::{Lattice, Site, StructureBuilder};
use xtalforge_ml::{
    GraphOptions, ImportFormat, ImportOptions, SplitOptions, SplitStrategy, crystal_graph,
    import_dataset, import_extxyz, split_dataset,
};

fn nacl() -> xtalforge_core::Structure {
    let lattice = Lattice::from_parameters(5.64, 5.64, 5.64, 90.0, 90.0, 90.0).expect("cell");
    StructureBuilder::new("NaCl", lattice)
        .site(Site::new("Na1", "Na", [0.0, 0.0, 0.0]).expect("Na"))
        .site(Site::new("Cl1", "Cl", [0.5, 0.5, 0.5]).expect("Cl"))
        .build()
        .expect("structure")
}

#[test]
fn graph_has_symmetric_edges_and_images() {
    let graph = crystal_graph(
        &nacl(),
        &GraphOptions {
            cutoff: 6.0,
            max_neighbors: None,
        },
    )
    .expect("graph");
    assert_eq!(graph.num_nodes(), 2);
    assert!(graph.num_edges() > 0);
    for [from, to] in graph.edge_index() {
        let forward = graph
            .edge_from
            .iter()
            .zip(&graph.edge_to)
            .position(|(a, b)| *a == from && *b == to)
            .expect("edge");
        assert!(graph.distances[forward] > 0.0);
        assert!(
            graph
                .edge_from
                .iter()
                .zip(&graph.edge_to)
                .any(|(a, b)| *a == to && *b == from),
            "reverse edge missing"
        );
    }
    assert!(!graph.truncated);
}

#[test]
fn graph_knn_truncates_per_node() {
    let graph = crystal_graph(
        &nacl(),
        &GraphOptions {
            cutoff: 6.0,
            max_neighbors: Some(2),
        },
    )
    .expect("graph");
    let mut counts = std::collections::BTreeMap::new();
    for from in &graph.edge_from {
        *counts.entry(*from).or_insert(0usize) += 1;
    }
    assert!(counts.values().all(|count| *count <= 2));
    assert!(graph.truncated);
}

#[test]
fn graph_rejects_bad_options() {
    assert!(
        crystal_graph(
            &nacl(),
            &GraphOptions {
                cutoff: 0.0,
                max_neighbors: None
            }
        )
        .is_err()
    );
    assert!(
        crystal_graph(
            &nacl(),
            &GraphOptions {
                cutoff: 4.0,
                max_neighbors: Some(0)
            }
        )
        .is_err()
    );
}

#[test]
fn splits_are_deterministic_and_cover_every_record() {
    let mut dataset = Dataset::new("split");
    for index in 0..10 {
        dataset
            .push(DatasetRecord::new(format!("r{index}"), nacl()).expect("record"))
            .expect("push");
    }
    let options = SplitOptions {
        train: 0.6,
        validation: 0.2,
        test: 0.2,
        seed: 7,
        strategy: SplitStrategy::Random,
    };
    let first = split_dataset(&dataset, &options).expect("split");
    let second = split_dataset(&dataset, &options).expect("split");
    assert_eq!(first, second);
    assert_eq!(first.labels.len(), 10);
    let summary = first.summary();
    assert_eq!(summary.train + summary.validation + summary.test, 10);
}

#[test]
fn composition_split_keeps_groups_whole() {
    let mut dataset = Dataset::new("composition");
    for index in 0..8 {
        let structure = nacl();
        dataset
            .push(DatasetRecord::new(format!("r{index}"), structure).expect("record"))
            .expect("push");
    }
    let assignment = split_dataset(
        &dataset,
        &SplitOptions {
            seed: 1,
            strategy: SplitStrategy::Composition,
            ..SplitOptions::default()
        },
    )
    .expect("split");
    // All records have the same composition, so they share one group.
    assert_eq!(assignment.n_groups, 1);
    assert!(assignment.labels.windows(2).all(|pair| pair[0] == pair[1]));
}

#[test]
fn extxyz_import_reads_frames_and_properties() {
    let text = "2\nLattice=\"5 0 0 0 5 0 0 0 5\" Properties=species:S:1:pos:R:3:energy:R:1 pbc=\"T T T\"\nNa 0 0 0 -1.5\nCl 2.5 2.5 2.5 -1.6\n2\nLattice=\"5 0 0 0 5 0 0 0 5\" Properties=species:S:1:pos:R:3:energy:R:1 pbc=\"T T T\"\nNa 0.1 0 0 -1.7\nCl 2.4 2.5 2.5 -1.8\n";
    let report = import_extxyz(text, "extxyz", None).expect("import");
    assert_eq!(report.n_records, 2);
    assert_eq!(report.skipped, 0);
    assert_eq!(report.dataset.len(), 2);
    let first = report.dataset.iter().next().expect("first");
    assert_eq!(first.id, "extxyz-0");
}

#[test]
fn json_and_optimade_imports_work() {
    let structure = nacl();
    let document = structure.to_document();
    let record = serde_json::json!({
        "id": "one",
        "structure": serde_json::from_str::<serde_json::Value>(&document.to_json().expect("json")).expect("value"),
        "properties": {"energy": -1.5, "label": "x"}
    });
    let text = serde_json::to_string(&vec![record]).expect("array");
    let report = import_dataset(&text, &ImportOptions::new(ImportFormat::Json, "json"))
        .expect("json import");
    assert_eq!(report.n_records, 1);

    let optimade = serde_json::json!({
        "data": [{
            "id": "mp-1",
            "attributes": {
                "lattice_vectors": [[5.0,0.0,0.0],[0.0,5.0,0.0],[0.0,0.0,5.0]],
                "cartesian_site_positions": [[0.0,0.0,0.0],[2.5,2.5,2.5]],
                "species_at_sites": ["Na","Cl"],
                "nelements": 2
            }
        }]
    });
    let report = import_dataset(
        &optimade.to_string(),
        &ImportOptions::new(ImportFormat::Optimade, "optimade"),
    )
    .expect("optimade");
    assert_eq!(report.n_records, 1);
    assert_eq!(report.dataset.len(), 1);
}

#[test]
fn invalid_imports_are_typed() {
    assert!(import_extxyz("not a frame", "bad", None).is_err());
    assert!(import_dataset("not json", &ImportOptions::new(ImportFormat::Json, "bad")).is_err());
}
