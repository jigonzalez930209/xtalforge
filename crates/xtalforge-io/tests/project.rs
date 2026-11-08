//! Integration tests for the `.xtalproj` container (Task 3.7).

#![allow(
    clippy::unwrap_used,
    clippy::expect_used,
    clippy::panic,
    clippy::indexing_slicing
)]

use serde_json::json;
use xtalforge_core::{Lattice, Site, StructureBuilder};
use xtalforge_io::project::{FORMAT_VERSION, ProjectParts, read_project, write_project};

fn silicon() -> xtalforge_core::Structure {
    StructureBuilder::new("Si", Lattice::cubic(5.430941).unwrap())
        .site(Site::new("Si1", "Si", [0.0, 0.0, 0.0]).unwrap())
        .build()
        .unwrap()
}

#[test]
fn project_round_trips_every_part() {
    let parts = ProjectParts::from_structure(
        &silicon(),
        Some(json!({"a": 5.430941})),
        Some(json!({"bonds": []})),
        Some(json!({"sites": []})),
        Some(json!({"style": "ball-stick"})),
        Some(json!({"selection": [0]})),
    )
    .unwrap();
    let bytes = write_project(&parts).unwrap();
    assert_eq!(&bytes[..2], b"PK");
    let back = read_project(&bytes).unwrap();
    assert_eq!(back, parts);
    assert_eq!(back.structure().unwrap().formula(), "Si");
    assert_eq!(FORMAT_VERSION, 1);
}

#[test]
fn missing_or_newer_versions_are_handled() {
    let parts = ProjectParts::from_structure(&silicon(), None, None, None, None, None).unwrap();
    let bytes = write_project(&parts).unwrap();
    assert!(read_project(b"not a zip").is_err());
    assert!(read_project(&bytes).unwrap().cell.is_none());
}
