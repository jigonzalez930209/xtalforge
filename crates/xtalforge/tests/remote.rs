//! Remote provider parsing (Task 13.3).

#![allow(
    clippy::unwrap_used,
    clippy::expect_used,
    clippy::panic,
    clippy::indexing_slicing
)]

use xtalforge::remote::{cod_from_cif, materials_project_from_json, optimade_from_json};

const OPTIMADE: &str = r#"{
  "meta": {"provider": "test-provider"},
  "links": {"next": "https://example.test/structures?page_offset=2"},
  "data": [{
    "id": "Si-1",
    "type": "structures",
    "attributes": {
      "lattice_vectors": [[5.43, 0, 0], [0, 5.43, 0], [0, 0, 5.43]],
      "cartesian_site_positions": [[0, 0, 0], [1.3575, 1.3575, 1.3575]],
      "species_at_sites": ["Si", "Si"],
      "license": "CC BY 4.0"
    }
  }]
}"#;

const MP: &str = r#"{
  "data": [{
    "material_id": "mp-149",
    "structure": {
      "lattice": {"matrix": [[5.43, 0, 0], [0, 5.43, 0], [0, 0, 5.43]]},
      "sites": [
        {"species": "Si", "xyz": [0, 0, 0]},
        {"species": "Si", "xyz": [1.3575, 1.3575, 1.3575]}
      ]
    }
  }]
}"#;

#[test]
fn optimade_builds_structures_and_next_link() {
    let page = optimade_from_json(OPTIMADE).unwrap();
    assert_eq!(page.structures.len(), 1);
    assert_eq!(page.structures[0].id, "Si-1");
    assert_eq!(page.structures[0].source, "test-provider");
    assert_eq!(page.structures[0].license.as_deref(), Some("CC BY 4.0"));
    assert_eq!(page.structures[0].structure.sites().len(), 2);
    assert!((page.structures[0].structure.lattice().volume() - 160.1).abs() < 0.5);
    assert!(page.next.is_some());
}

#[test]
fn optimade_rejects_missing_fields() {
    assert!(optimade_from_json("{}").is_err());
    assert!(optimade_from_json(r#"{"data": [{"attributes": {}}]}"#).is_err());
    assert!(optimade_from_json("not json").is_err());
}

#[test]
fn materials_project_builds_a_structure() {
    let remote = materials_project_from_json(MP).unwrap();
    assert_eq!(remote.id, "mp-149");
    assert_eq!(remote.source, "materials-project");
    assert_eq!(remote.structure.sites().len(), 2);
    assert_eq!(remote.structure.formula(), "Si2");
}

#[test]
fn cod_cif_builds_a_structure() {
    let cif = include_str!("../../../tests/corpus/simple/Si.cif");
    let remote = cod_from_cif(cif, "1000041").unwrap();
    assert_eq!(remote.id, "cod:1000041");
    assert_eq!(remote.source, "cod");
    assert_eq!(remote.structure.sites().len(), 1);
}

#[test]
fn malformed_materials_project_is_an_error() {
    assert!(materials_project_from_json("{}").is_err());
    assert!(materials_project_from_json(r#"{"data": [{}]}"#).is_err());
}

#[test]
fn optimade_export_round_trips_through_the_parser() {
    let cif = include_str!("../../../tests/corpus/simple/Si.cif");
    let structure = xtalforge::io::read_cif(cif).unwrap();
    let entry = xtalforge::remote::optimade_export(&structure).unwrap();
    let document: serde_json::Value = serde_json::from_str(&entry).unwrap();
    assert_eq!(document["attributes"]["species_at_sites"][0], "Si");
    let wrapped = format!(r#"{{"data": [{}]}}"#, entry);
    let page = optimade_from_json(&wrapped).unwrap();
    assert_eq!(page.structures[0].structure.sites().len(), 8);
    assert_eq!(page.structures[0].structure.formula(), structure.formula());
}

#[test]
fn optimade_is_a_registered_writer() {
    let cif = include_str!("../../../tests/corpus/simple/Si.cif");
    let structure = xtalforge::io::read_cif(cif).unwrap();
    let text = xtalforge::io::write_structure(&structure, "optimade").unwrap();
    let value: serde_json::Value = serde_json::from_str(&text).unwrap();
    assert_eq!(value["type"], "structures");
    assert_eq!(
        xtalforge::io::registry::find("optimade").unwrap().name,
        "optimade"
    );
}
