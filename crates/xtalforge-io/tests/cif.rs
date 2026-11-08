//! Integration tests for the CIF reader and writer (Task 1.15).

#![allow(
    clippy::unwrap_used,
    clippy::expect_used,
    clippy::panic,
    clippy::indexing_slicing
)]

use xtalforge_core::{Lattice, Site, StructureBuilder};
use xtalforge_io::cif::{CifDocument, read_cif, read_cif_all, read_cif_block, write_cif};

const SI: &str = include_str!("../../../tests/corpus/simple/Si.cif");
const MOF5: &str = include_str!("../../../tests/corpus/cod/MOF-5.cif");
const MALFORMED: &str = include_str!("../../../tests/corpus/edge/malformed.cif");

#[test]
fn reads_silicon_and_expands_to_p1() {
    let structure = read_cif(SI).unwrap();
    assert_eq!(structure.name(), "9011998");
    assert!((structure.lattice().lengths()[0] - 5.430_941).abs() < 1e-6);
    assert_eq!(structure.formula(), "Si");
    let p1 = structure.p1().unwrap();
    assert_eq!(p1.n_atoms(), 8);
    assert!(structure.operations().len() > 1);
}

#[test]
fn reads_mof5_and_matches_the_golden_atom_count() {
    let structure = read_cif(MOF5).unwrap();
    let p1 = structure.p1().unwrap();
    assert!(p1.n_atoms() > 0);
    let golden: serde_json::Value =
        serde_json::from_str(include_str!("../../../tests/golden/phase-1/cod/MOF-5.json")).unwrap();
    let expected = golden["p1"]["n_atoms"].as_u64().unwrap() as usize;
    assert_eq!(p1.n_atoms(), expected);
}

#[test]
fn malformed_cif_returns_an_error() {
    let error = read_cif(MALFORMED).expect_err("malformed CIF must fail");
    assert_eq!(error.code(), "invalid_input");
}

#[test]
fn block_selection_and_multiple_blocks() {
    let text = "data_one\n_cell_length_a 4\n_cell_length_b 4\n_cell_length_c 4\n\
                _cell_angle_alpha 90\n_cell_angle_beta 90\n_cell_angle_gamma 90\n\
                loop_\n_atom_site_label\n_atom_site_type_symbol\n_atom_site_fract_x\n\
                _atom_site_fract_y\n_atom_site_fract_z\nC1 C 0 0 0\n\
                data_two\n_cell_length_a 5\n_cell_length_b 5\n_cell_length_c 5\n\
                _cell_angle_alpha 90\n_cell_angle_beta 90\n_cell_angle_gamma 90\n\
                loop_\n_atom_site_label\n_atom_site_type_symbol\n_atom_site_fract_x\n\
                _atom_site_fract_y\n_atom_site_fract_z\nC2 C 0 0 0\n";
    let all = read_cif_all(text).unwrap();
    assert_eq!(all.len(), 2);
    assert_eq!(read_cif_block(text, "two").unwrap().name(), "two");
    assert!(read_cif_block(text, "missing").is_err());
}

#[test]
fn parser_handles_comments_quotes_and_uncertainties() {
    let text = "data_test\n_cell_length_a 4.0(2) # comment\n\
                _cell_length_b 4.0\n_cell_length_c 4.0\n\
                _cell_angle_alpha 90\n_cell_angle_beta 90\n_cell_angle_gamma 90\n\
                loop_\n_atom_site_label\n_atom_site_type_symbol\n\
                _atom_site_fract_x\n_atom_site_fract_y\n_atom_site_fract_z\n\
                'C 1' C 0.0 0.0 0.0\n";
    let document = CifDocument::parse(text).unwrap();
    assert_eq!(document.blocks[0].item("_cell_length_a"), Some("4.0(2)"));
    let structure = read_cif(text).unwrap();
    assert_eq!(structure.formula(), "C");
}

#[test]
fn writer_round_trips_a_p1_structure() {
    let structure = StructureBuilder::new("C", Lattice::cubic(4.0).unwrap())
        .site(
            Site::new("C1", "C", [0.1, 0.2, 0.3])
                .unwrap()
                .with_occupancy(0.5)
                .unwrap()
                .with_u_iso(0.01)
                .unwrap()
                .with_charge(-1.0)
                .unwrap(),
        )
        .build()
        .unwrap();
    let text = write_cif(&structure).unwrap();
    let back = read_cif(&text).unwrap();
    assert_eq!(back.formula(), "C0.500");
    assert_eq!(back.sites().len(), 1);
    assert_eq!(back.sites()[0].occupancy(), 0.5);
    assert_eq!(back.sites()[0].charge(), Some(-1.0));
}
