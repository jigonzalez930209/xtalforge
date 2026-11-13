//! Integration tests for provenance embedding and extraction (Task 16.2.2).

#![allow(
    clippy::unwrap_used,
    clippy::expect_used,
    clippy::panic,
    clippy::indexing_slicing
)]

use serde_json::json;
use xtalforge_core::{Lattice, PROVENANCE_METADATA_KEY, Provenance, Site, StructureBuilder};
use xtalforge_io::project::{ProjectParts, read_project, write_project};
use xtalforge_io::{provenance, read_cif, read_xyz, write_cif, write_xyz};

const CIF: &str = "\
data_plain
_cell_length_a 4
_cell_length_b 4
_cell_length_c 4
_cell_angle_alpha 90
_cell_angle_beta 90
_cell_angle_gamma 90
loop_
_atom_site_label
_atom_site_type_symbol
_atom_site_fract_x
_atom_site_fract_y
_atom_site_fract_z
C1 C 0 0 0
";

const XYZ: &str = "\
1
Lattice=\"4.0 0.0 0.0 0.0 4.0 0.0 0.0 0.0 4.0\" Properties=species:S:1:pos:R:3
C 0.0 0.0 0.0
";

fn structure() -> xtalforge_core::Structure {
    StructureBuilder::new("Si", Lattice::cubic(5.43).unwrap())
        .site(Site::new("Si1", "Si", [0.0, 0.0, 0.0]).unwrap())
        .build()
        .unwrap()
}

fn record() -> Provenance {
    Provenance::new("relax")
        .with_option("optimizer", "lbfgs")
        .with_option("fmax", 0.01)
        .with_input_bytes(Some("input.cif".to_owned()), b"crystal")
        .with_bundle("uff", "data-2026.11.13")
        .with_seed(42)
        .with_timestamp("2026-11-13T12:00:00Z")
}

#[test]
fn cif_round_trips_provenance() {
    let record = record();
    let structure = record.attach(&structure()).unwrap();
    let text = write_cif(&structure).unwrap();
    assert!(text.contains("_xtalforge_provenance"));
    assert!(text.contains("\n;{\"schema_version\":1,"));

    let read = read_cif(&text).unwrap();
    assert_eq!(Provenance::from_structure(&read).unwrap().unwrap(), record);
    assert_eq!(provenance(&text).unwrap().unwrap(), record);
}

#[test]
fn cif_without_provenance_returns_none() {
    assert!(provenance(CIF).unwrap().is_none());
    let structure = read_cif(CIF).unwrap();
    assert!(Provenance::from_structure(&structure).unwrap().is_none());
    assert!(provenance(&structure).unwrap().is_none());
}

#[test]
fn cif_comment_is_read_even_when_it_is_the_only_carrier() {
    let record = record();
    let text = CIF.replacen(
        '\n',
        &format!("\n# xtalforge-provenance: {}\n", record.to_json().unwrap()),
        1,
    );
    assert_eq!(provenance(&text).unwrap().unwrap(), record);
}

#[test]
fn extxyz_round_trips_provenance() {
    let record = record();
    let structure = record.attach(&structure()).unwrap();
    let text = write_xyz(&structure).unwrap();
    assert!(text.contains("xtalforge_provenance="));
    assert!(!text.contains("\"schema_version\""));

    let read = read_xyz(&text).unwrap();
    assert_eq!(Provenance::from_structure(&read).unwrap().unwrap(), record);
    assert_eq!(provenance(&text).unwrap().unwrap(), record);
    assert!(provenance(XYZ).unwrap().is_none());
}

#[test]
fn extxyz_provenance_is_found_in_a_later_frame() {
    let plain = write_xyz(&structure()).unwrap();
    let with_record = write_xyz(&record().attach(&structure()).unwrap()).unwrap();
    let text = format!("{plain}{with_record}");
    assert_eq!(provenance(&text).unwrap().unwrap(), record());
}

#[test]
fn project_round_trips_provenance_and_keeps_the_session() {
    let record = record();
    let parts = ProjectParts::from_structure(
        &structure(),
        None,
        None,
        None,
        None,
        Some(json!({"selection": [0]})),
    )
    .unwrap()
    .with_provenance(&record)
    .unwrap();
    let bytes = write_project(&parts).unwrap();
    let read = read_project(&bytes).unwrap();
    assert_eq!(read.provenance().unwrap().unwrap(), record);
    assert_eq!(read.session.as_ref().unwrap()["selection"], json!([0]));
    assert_eq!(provenance(&bytes).unwrap().unwrap(), record);
    assert_eq!(provenance(bytes.as_slice()).unwrap().unwrap(), record);

    let plain = ProjectParts::from_structure(&structure(), None, None, None, None, None).unwrap();
    assert!(plain.provenance().unwrap().is_none());
    assert!(
        read_project(&write_project(&plain).unwrap())
            .unwrap()
            .provenance()
            .unwrap()
            .is_none()
    );
}

#[test]
fn structure_json_carries_provenance() {
    let record = record();
    let structure = record.attach(&structure()).unwrap();
    let document = structure.to_document().to_json().unwrap();
    assert_eq!(provenance(&document).unwrap().unwrap(), record);
    let bare = serde_json::to_string(&structure).unwrap();
    assert_eq!(provenance(&bare).unwrap().unwrap(), record);
}

#[test]
fn malformed_payloads_are_rejected() {
    let broken_xyz = XYZ.replace(
        "Properties=species:S:1:pos:R:3",
        "Properties=species:S:1:pos:R:3 xtalforge_provenance=***",
    );
    assert!(provenance(&broken_xyz).is_err());
    assert!(read_xyz(&broken_xyz).is_err());

    let broken_cif = CIF.replacen('\n', "\n_xtalforge_provenance\n;not json\n;\n", 1);
    assert!(provenance(&broken_cif).is_err());

    let invalid_session = ProjectParts::from_structure(
        &structure(),
        None,
        None,
        None,
        None,
        Some(json!({"provenance": {"schema_version": 99}})),
    )
    .unwrap();
    assert!(invalid_session.provenance().is_err());

    let non_object = ProjectParts::from_structure(
        &structure(),
        None,
        None,
        None,
        None,
        Some(json!("not an object")),
    )
    .unwrap();
    assert!(non_object.with_provenance(&record()).is_err());
    assert!(provenance("not a structure at all").is_err());
    assert!(provenance("").is_err());
}

#[test]
fn tampered_cif_input_is_detected_after_extraction() {
    let record =
        Provenance::new("read").with_input_bytes(Some("input.cif".to_owned()), b"original");
    let structure = record.attach(&structure()).unwrap();
    let text = write_cif(&structure).unwrap();
    let extracted = provenance(&text).unwrap().unwrap();
    extracted
        .verify_input(Some("input.cif"), b"original")
        .unwrap();
    assert!(
        extracted
            .verify_input(Some("input.cif"), b"original ")
            .is_err()
    );
    assert!(
        extracted
            .attach(&structure)
            .unwrap()
            .metadata()
            .get(PROVENANCE_METADATA_KEY)
            .is_some()
    );
}
