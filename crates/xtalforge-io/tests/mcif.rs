//! Integration tests for the magCIF reader and writer (Task 19.1.3).
//!
//! The fixtures are hand-written for these tests from the public magCIF tag
//! definitions; no third-party data is included. They are dedicated to the
//! public domain under CC0-1.0.

#![allow(
    clippy::unwrap_used,
    clippy::expect_used,
    clippy::panic,
    clippy::indexing_slicing
)]

use xtalforge_core::{Lattice, MagneticAxis, MagneticMoment, Site, StructureBuilder};
use xtalforge_io::registry::{find, sniff};
use xtalforge_io::{
    MAGNETIC_OPERATIONS_KEY, MagneticSymOp, read_mcif, read_mcif_block, read_structure, write_mcif,
    write_structure,
};

const MCIF: &str = r#"# Hand-written magCIF fixture, CC0-1.0.
data_mn_test
_chemical_formula_sum 'Mn2 O'
_cell_length_a 5.0
_cell_length_b 4.0
_cell_length_c 8.0
_cell_angle_alpha 90
_cell_angle_beta 90
_cell_angle_gamma 90
_space_group_name_H-M 'P 1'
_space_group_name_Hall 'P 1'
loop_
_space_group_symop_magn_operation_xyz
'x,y,z,+1'
'-x,-y,z,-1'
loop_
_atom_site_label
_atom_site_type_symbol
_atom_site_fract_x
_atom_site_fract_y
_atom_site_fract_z
_atom_site_occupancy
Mn1 Mn 0.0 0.0 0.0 1.0
Mn2 Mn 0.5 0.5 0.0 0.5
O1 O 0.25 0.25 0.25 1.0
loop_
_atom_site_moment_label
_atom_site_moment_crystalaxis_x
_atom_site_moment_crystalaxis_y
_atom_site_moment_crystalaxis_z
Mn1 0.0 0.0 4.0
Mn2 0.0 0.0 -4.0
O1 0.5 -0.25 1.0
"#;

const CENTERED: &str = r#"# Hand-written magCIF fixture, CC0-1.0.
data_centered
_cell_length_a 4.0
_cell_length_b 4.0
_cell_length_c 4.0
_cell_angle_alpha 90
_cell_angle_beta 90
_cell_angle_gamma 90
_space_group_name_Hall 'P 1'
loop_
_space_group_symop_magn_operation_xyz
'x,y,z,+1'
loop_
_space_group_symop_magn_centering_xyz
'0,0,0'
'1/2,1/2,1/2'
loop_
_atom_site_label
_atom_site_type_symbol
_atom_site_fract_x
_atom_site_fract_y
_atom_site_fract_z
Fe1 Fe 0.0 0.0 0.0
"#;

const CARTESIAN: &str = r#"# Hand-written magCIF fixture, CC0-1.0.
data_cartesian
_cell_length_a 4.0
_cell_length_b 2.0
_cell_length_c 8.0
_cell_angle_alpha 90
_cell_angle_beta 90
_cell_angle_gamma 90
_space_group_name_Hall 'P 1'
loop_
_atom_site_label
_atom_site_type_symbol
_atom_site_fract_x
_atom_site_fract_y
_atom_site_fract_z
C1 C 0.0 0.0 0.0
C2 C 0.5 0.5 0.5
loop_
_atom_site_moment_label
_atom_site_moment_cartn_x
_atom_site_moment_cartn_y
_atom_site_moment_cartn_z
C1 4.0 0.0 0.0
C2 0.0 2.0 0.0
"#;

fn magnetic_operations(text: &str) -> Vec<MagneticSymOp> {
    let structure = read_mcif(text).unwrap();
    serde_json::from_value(structure.metadata()[MAGNETIC_OPERATIONS_KEY].clone()).unwrap()
}

#[test]
fn reads_moments_operations_and_occupancies() {
    let structure = read_mcif(MCIF).unwrap();
    assert_eq!(structure.name(), "mn_test");
    assert_eq!(structure.sites().len(), 3);
    assert_eq!(
        structure.sites()[0].magnetic_moment(),
        Some(MagneticMoment::collinear(MagneticAxis::C, 4.0).unwrap())
    );
    assert_eq!(
        structure.sites()[1].magnetic_moment(),
        Some(MagneticMoment::collinear(MagneticAxis::C, -4.0).unwrap())
    );
    assert!((structure.sites()[1].occupancy() - 0.5).abs() < 1e-12);

    let oxygen = structure.sites()[2].magnetic_moment().unwrap();
    assert!(!oxygen.is_collinear());
    let [x, y, z] = oxygen.components();
    assert!((x - 0.5).abs() < 1e-12);
    assert!((y + 0.25).abs() < 1e-12);
    assert!((z - 1.0).abs() < 1e-12);

    let operations = magnetic_operations(MCIF);
    assert_eq!(operations.len(), 2);
    assert!(!operations[0].time_reversal);
    assert!(operations[1].time_reversal);
    assert_eq!(
        operations[1].operation.rotation(),
        [[-1, 0, 0], [0, -1, 0], [0, 0, 1]]
    );
}

#[test]
fn write_and_read_round_trip() {
    let structure = read_mcif(MCIF).unwrap();
    let text = write_mcif(&structure).unwrap();
    assert!(text.starts_with("data_mn_test"));
    assert!(text.contains("_space_group_symop_magn_operation_xyz"));
    assert!(text.contains("'-x,-y,z,-1'"));
    assert!(text.contains("_atom_site_moment_crystalaxis_z"));
    let back = read_mcif(&text).unwrap();
    assert_eq!(back.sites().len(), structure.sites().len());
    for (expected, actual) in structure.sites().iter().zip(back.sites()) {
        assert_eq!(expected.label(), actual.label());
        assert_eq!(expected.magnetic_moment(), actual.magnetic_moment());
        assert_eq!(expected.occupancy(), actual.occupancy());
    }
    let operations: Vec<MagneticSymOp> =
        serde_json::from_value(back.metadata()[MAGNETIC_OPERATIONS_KEY].clone()).unwrap();
    assert_eq!(
        operations,
        serde_json::from_value::<Vec<MagneticSymOp>>(
            structure.metadata()[MAGNETIC_OPERATIONS_KEY].clone()
        )
        .unwrap()
    );
}

#[test]
fn non_magnetic_structures_stay_plain_cif() {
    let structure = StructureBuilder::new("C", Lattice::cubic(4.0).unwrap())
        .site(Site::new("C1", "C", [0.0, 0.0, 0.0]).unwrap())
        .build()
        .unwrap();
    let text = write_mcif(&structure).unwrap();
    assert!(!text.contains("_space_group_symop_magn"));
    assert!(!text.contains("_atom_site_moment"));
    let back = read_mcif(&text).unwrap();
    assert_eq!(back.name(), structure.name());
    assert_eq!(back.lattice(), structure.lattice());
    assert_eq!(back.sites(), structure.sites());
    assert!(!back.metadata().contains_key(MAGNETIC_OPERATIONS_KEY));
}

#[test]
fn centering_translations_expand_the_operation_loop() {
    let operations = magnetic_operations(CENTERED);
    assert_eq!(operations.len(), 2);
    assert_eq!(operations[1].operation.translation(), [0.5, 0.5, 0.5]);
    assert!(!operations[1].time_reversal);
}

#[test]
fn ddlm_dot_spellings_are_accepted() {
    let dotted = MCIF
        .replace(
            "_space_group_symop_magn_operation_xyz",
            "_space_group_symop.magn_operation_xyz",
        )
        .replace(
            "_atom_site_moment_crystalaxis_",
            "_atom_site_moment.crystalaxis_",
        )
        .replace("_atom_site_moment_label", "_atom_site_moment.label");
    let structure = read_mcif(&dotted).unwrap();
    assert_eq!(structure.sites().len(), 3);
    assert_eq!(
        structure.sites()[0].magnetic_moment(),
        Some(MagneticMoment::collinear(MagneticAxis::C, 4.0).unwrap())
    );
}

#[test]
fn cartesian_moments_are_projected_onto_the_lattice_axes() {
    let structure = read_mcif(CARTESIAN).unwrap();
    assert_eq!(
        structure.sites()[0].magnetic_moment(),
        Some(MagneticMoment::collinear(MagneticAxis::A, 4.0).unwrap())
    );
    assert_eq!(
        structure.sites()[1].magnetic_moment(),
        Some(MagneticMoment::collinear(MagneticAxis::B, 2.0).unwrap())
    );
    let plain = CARTESIAN.replace("_atom_site_moment_cartn_", "_atom_site_moment_");
    let structure = read_mcif(&plain).unwrap();
    assert_eq!(
        structure.sites()[0].magnetic_moment(),
        Some(MagneticMoment::collinear(MagneticAxis::A, 4.0).unwrap())
    );
}

#[test]
fn malformed_magnetic_input_is_rejected() {
    assert!(read_mcif("").is_err());
    assert!(read_mcif("data_x\n_cell_length_a 4\n").is_err());
    let bad_operation = MCIF.replace("'-x,-y,z,-1'", "'-x,-y'");
    assert!(read_mcif(&bad_operation).is_err());
    let bad_moment = MCIF.replace("Mn1 0.0 0.0 4.0", "Mn1 0.0 0.0 nope");
    assert!(read_mcif(&bad_moment).is_err());
    let bad_centering = CENTERED.replace("'1/2,1/2,1/2'", "'1/2,1/2'");
    assert!(read_mcif(&bad_centering).is_err());
    assert_eq!(read_mcif_block(MCIF, "mn_test").unwrap().name(), "mn_test");
    assert!(
        read_mcif_block(MCIF, "missing")
            .unwrap_err()
            .message()
            .contains("not found")
    );
}

#[test]
fn mcif_is_registered_and_sniffed() {
    let format = find(".mcif").unwrap();
    assert_eq!(format.name, "mcif");
    assert!(format.read.is_some());
    assert!(format.write.is_some());
    assert_eq!(find("mcif").unwrap().name, "mcif");
    assert_eq!(find(".cif").unwrap().name, "cif");
    assert_eq!(sniff(MCIF).unwrap().name, "mcif");

    let structure = read_structure(MCIF, "mcif").unwrap();
    assert_eq!(
        structure.sites()[0].magnetic_moment(),
        Some(MagneticMoment::collinear(MagneticAxis::C, 4.0).unwrap())
    );
    let text = write_structure(&structure, "mcif").unwrap();
    assert!(text.contains("_atom_site_moment_crystalaxis_z"));
    assert!(
        read_structure(MCIF, "cif").unwrap().sites()[0]
            .magnetic_moment()
            .is_none()
    );
}
