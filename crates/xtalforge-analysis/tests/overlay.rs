//! Integration tests for overlay data (Task 5.9.1).

#![allow(
    clippy::unwrap_used,
    clippy::expect_used,
    clippy::panic,
    clippy::indexing_slicing
)]

use xtalforge_analysis::overlay::{OverlaySet, OverlayStyle, overlay_net};
use xtalforge_analysis::{NetEdge, QuotientGraph};

fn pair() -> QuotientGraph {
    QuotientGraph {
        nodes: vec![[0.0, 0.0, 0.0], [0.5, 0.5, 0.0]],
        edges: vec![NetEdge {
            a: 0,
            b: 1,
            shift: [0; 3],
        }],
        labels: vec!["A".to_owned(), "B".to_owned()],
    }
}

#[test]
fn a_net_overlay_uses_cartesian_sites() {
    let cell = [[2.0, 0.0, 0.0], [0.0, 3.0, 0.0], [0.0, 0.0, 4.0]];
    let style = OverlayStyle {
        color: [10, 20, 30],
        radius: 0.4,
    };
    let overlay = overlay_net(&pair(), cell, style);
    assert_eq!(overlay.name, "net");
    assert_eq!(overlay.cell, cell);
    assert_eq!(overlay.sites.len(), 2);
    assert_eq!(overlay.sites[0].0, [0.0, 0.0, 0.0]);
    assert_eq!(overlay.sites[1].0, [1.0, 1.5, 0.0]);
    assert_eq!(overlay.sites[0].1, "A");
    assert_eq!(overlay.sites[1].1, "B");
    assert_eq!(overlay.style, style);
}

#[test]
fn overlay_sets_are_plain_json_data() {
    let mut set = OverlaySet::new();
    assert!(set.is_empty());
    set.push(overlay_net(
        &pair(),
        [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]],
        OverlayStyle {
            color: [255, 0, 0],
            radius: 0.25,
        },
    ));
    assert_eq!(set.len(), 1);
    let text = serde_json::to_string(&set).unwrap();
    let read_back: OverlaySet = serde_json::from_str(&text).unwrap();
    assert_eq!(read_back, set);
    assert!(text.contains("\"overlays\""));
    assert!(text.contains("\"color\":[255,0,0]"));
}

#[test]
fn missing_labels_get_a_stable_fallback() {
    let mut graph = pair();
    graph.labels.clear();
    let overlay = overlay_net(
        &graph,
        [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]],
        OverlayStyle {
            color: [0, 0, 0],
            radius: 0.1,
        },
    );
    assert_eq!(overlay.sites[0].1, "node 1");
    assert_eq!(overlay.sites[1].1, "node 2");
}
