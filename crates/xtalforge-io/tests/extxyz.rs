//! Integration tests for extended XYZ properties and the plain XYZ behaviour
//! (Task 13.2.6).

#![allow(
    clippy::unwrap_used,
    clippy::expect_used,
    clippy::panic,
    clippy::indexing_slicing
)]

use serde_json::Value;
use xtalforge_io::registry::{find, read_with, sniff};
use xtalforge_io::{read_structure, read_xyz, write_xyz};

const EXTXYZ: &str = "\
2
Lattice=\"4.0 0.0 0.0 0.0 4.0 0.0 0.0 0.0 4.0\" Properties=species:S:1:pos:R:3:forces:R:3:charge:R:1:fixed:L:1 energy=-1.5 step=7 pbc=\"T T T\"
C 0.0 0.0 0.0 0.1 0.2 0.3 -0.4 T
O 1.0 0.0 0.0 0.0 0.0 0.0 0.4 F
";

const LABELS: &str = "\
2
Lattice=\"4.0 0.0 0.0 0.0 4.0 0.0 0.0 0.0 4.0\" Properties=species:S:1:pos:R:3:label:S:1
C 0.0 0.0 0.0 C_first
O 1.0 0.0 0.0 O_second
";

#[test]
fn extxyz_keeps_arbitrary_atom_properties() {
    let structure = read_xyz(EXTXYZ).unwrap();
    assert_eq!(structure.sites().len(), 2);
    assert_eq!(structure.sites()[0].element(), "C");
    let forces = structure.sites()[0].props().get("forces").unwrap();
    assert_eq!(forces.as_array().unwrap().len(), 3);
    assert!((forces[0].as_f64().unwrap() - 0.1).abs() < 1e-9);
    assert!((forces[2].as_f64().unwrap() - 0.3).abs() < 1e-9);
    assert!((structure.sites()[0].props()["charge"].as_f64().unwrap() + 0.4).abs() < 1e-9);
    assert_eq!(
        structure.sites()[0].props().get("fixed"),
        Some(&Value::Bool(true))
    );
    assert_eq!(
        structure.sites()[1].props().get("fixed"),
        Some(&Value::Bool(false))
    );
}

#[test]
fn extxyz_keeps_per_frame_metadata() {
    let structure = read_xyz(EXTXYZ).unwrap();
    assert_eq!(structure.metadata().get("energy").unwrap(), "-1.5");
    assert_eq!(structure.metadata().get("step").unwrap(), "7");
    assert_eq!(structure.metadata().get("pbc").unwrap(), "T T T");
    assert!((structure.lattice().lengths()[0] - 4.0).abs() < 1e-12);
}

#[test]
fn extxyz_round_trip_preserves_properties_and_metadata() {
    let structure = read_xyz(EXTXYZ).unwrap();
    let text = write_xyz(&structure).unwrap();
    assert!(text.contains("Properties=species:S:1:pos:R:3"));
    assert!(text.contains("forces:R:3"));
    assert!(text.contains("fixed:L:1"));
    assert!(text.contains("energy=-1.5"));
    let back = read_xyz(&text).unwrap();
    assert_eq!(back.sites().len(), 2);
    let forces = back.sites()[0].props().get("forces").unwrap();
    assert!((forces[1].as_f64().unwrap() - 0.2).abs() < 1e-9);
    assert!((back.sites()[1].props()["charge"].as_f64().unwrap() - 0.4).abs() < 1e-9);
    assert_eq!(
        back.sites()[1].props().get("fixed"),
        Some(&Value::Bool(false))
    );
    assert_eq!(back.metadata().get("energy").unwrap(), "-1.5");
    assert_eq!(back.metadata().get("pbc").unwrap(), "T T T");
}

#[test]
fn extxyz_label_column_sets_the_site_label() {
    let structure = read_xyz(LABELS).unwrap();
    assert_eq!(structure.sites()[0].label(), "C_first");
    assert_eq!(structure.sites()[1].label(), "O_second");
    let text = write_xyz(&structure).unwrap();
    let back = read_xyz(&text).unwrap();
    assert_eq!(back.sites()[0].label(), "C_first");
}

#[test]
fn plain_xyz_keeps_its_padded_box_and_comment() {
    let text = "2\nwater\nH 0.0 0.0 0.0\nH 1.0 0.0 0.0\n";
    let structure = read_xyz(text).unwrap();
    assert_eq!(structure.formula(), "H2");
    assert!((structure.lattice().lengths()[0] - 11.0).abs() < 1e-9);
    assert!(structure.metadata().is_empty());
    let written = write_xyz(&structure).unwrap();
    assert!(written.starts_with("2\nLattice="));
}

#[test]
fn extxyz_reports_bad_declarations_and_rows() {
    assert!(read_xyz("1\nProperties=species:S:1:pos:R:3 C 0 0 0\n").is_err());
    let bad_type = EXTXYZ.replace("forces:R:3", "forces:X:3");
    assert!(read_xyz(&bad_type).is_err());
    let zero = EXTXYZ.replace("charge:R:1", "charge:R:0");
    assert!(read_xyz(&zero).is_err());
    let short = EXTXYZ.replace("C 0.0 0.0 0.0 0.1 0.2 0.3 -0.4 T", "C 0.0 0.0 0.0 0.1");
    assert!(read_xyz(&short).is_err());
    let bad_logical = EXTXYZ.replace(" -0.4 T", " -0.4 maybe");
    assert!(read_xyz(&bad_logical).is_err());
}

#[test]
fn extxyz_is_registered_and_dispatch_sniffs_it() {
    let format = find(".extxyz").unwrap();
    assert_eq!(format.name, "extxyz");
    assert!(format.read.is_some());
    assert!(format.write.is_some());
    assert!(format.capabilities.cell);
    assert_eq!(sniff(EXTXYZ).unwrap().name, "extxyz");
    assert_eq!(sniff("2\nplain\nH 0 0 0\nH 1 0 0\n").unwrap().name, "xyz");
    assert_eq!(read_with(EXTXYZ, "extxyz").unwrap().sites().len(), 2);
    assert_eq!(read_structure(EXTXYZ, "").unwrap().sites().len(), 2);
}
