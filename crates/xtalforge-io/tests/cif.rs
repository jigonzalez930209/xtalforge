//! Integration tests for the CIF reader and writer (Task 1.15).

#![allow(
    clippy::unwrap_used,
    clippy::expect_used,
    clippy::panic,
    clippy::indexing_slicing
)]

use xtalforge_core::{ExplicitBond, Lattice, Site, StructureBuilder};
use xtalforge_io::cif::{CifDocument, read_cif, read_cif_all, read_cif_block, write_cif};

const SI: &str = include_str!("../../../tests/corpus/simple/Si.cif");
const MOF5: &str = include_str!("../../../tests/corpus/cod/MOF-5.cif");
const MALFORMED: &str = include_str!("../../../tests/corpus/edge/malformed.cif");

#[test]
fn reads_silicon_and_expands_to_p1() {
    let structure = read_cif(SI).unwrap();
    assert_eq!(structure.name(), "9011998");
    assert!((structure.lattice().lengths()[0] - 5.430_941).abs() < 1e-6);
    assert_eq!(structure.formula(), "Si8");
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
fn symmetry_priority_uses_the_hall_symbol_when_operations_are_missing() {
    let text = "data_warn\n_cell_length_a 4\n_cell_length_b 4\n_cell_length_c 4\n\
                _cell_angle_alpha 90\n_cell_angle_beta 90\n_cell_angle_gamma 90\n\
                _symmetry_space_group_name_H-M 'P 21/c'\n\
                loop_\n_atom_site_label\n_atom_site_type_symbol\n\
                _atom_site_fract_x\n_atom_site_fract_y\n_atom_site_fract_z\n\
                C1 C 0 0 0\n";
    let structure = read_cif(text).unwrap();
    assert_eq!(structure.hall(), "P 21/c");
    assert_eq!(structure.operations().len(), 4);
    assert_eq!(
        structure
            .metadata()
            .get("symmetry_source")
            .and_then(|value| value.as_str()),
        Some("Hall symbol")
    );
    assert!(!structure.metadata().contains_key("symmetry_warning"));
}

#[test]
fn an_unknown_hall_symbol_still_warns() {
    let text = "data_warn\n_cell_length_a 4\n_cell_length_b 4\n_cell_length_c 4\n\
                _cell_angle_alpha 90\n_cell_angle_beta 90\n_cell_angle_gamma 90\n\
                _symmetry_space_group_name_Hall 'not-a-space-group'\n\
                loop_\n_atom_site_label\n_atom_site_type_symbol\n\
                _atom_site_fract_x\n_atom_site_fract_y\n_atom_site_fract_z\n\
                C1 C 0 0 0\n";
    let structure = read_cif(text).unwrap();
    assert!(structure.metadata().contains_key("symmetry_warning"));
}

#[test]
fn malformed_inputs_never_panic() {
    for text in [
        "",
        "data_",
        "loop_",
        "loop_\n_atom_site_label\nC1",
        "data_x\n_cell_length_a 4\n",
        "data_x\nloop_\n_a\n_b\n1\n",
        "data_x\n_cell_length_a 4\n_cell_length_b 4\n_cell_length_c 4\n\
         _cell_angle_alpha 90\n_cell_angle_beta 90\n_cell_angle_gamma 90\n\
         loop_\n_atom_site_label\n_atom_site_type_symbol\n_atom_site_fract_x\n\
         _atom_site_fract_y\n_atom_site_fract_z\nC1 C 0 0\n",
    ] {
        let _ = read_cif(text);
    }
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

#[test]
fn out_of_range_occupancy_is_clamped_and_reported() {
    let text = "data_clamp\n_cell_length_a 4\n_cell_length_b 4\n_cell_length_c 4\n\
                _cell_angle_alpha 90\n_cell_angle_beta 90\n_cell_angle_gamma 90\n\
                loop_\n_atom_site_label\n_atom_site_type_symbol\n\
                _atom_site_fract_x\n_atom_site_fract_y\n_atom_site_fract_z\n\
                _atom_site_occupancy\nC1 C 0 0 0 1.8\n";
    let structure = read_cif(text).unwrap();
    assert_eq!(structure.sites()[0].occupancy(), 1.0);
    assert_eq!(structure.metadata()["occupancy_clamped"].as_u64(), Some(1));
}

#[test]
fn round_trips_aniso_disorder_and_assembly() {
    let text = "data_aniso\n_cell_length_a 10\n_cell_length_b 10\n_cell_length_c 10\n\
                _cell_angle_alpha 90\n_cell_angle_beta 90\n_cell_angle_gamma 90\n\
                loop_\n_atom_site_label\n_atom_site_type_symbol\n\
                _atom_site_fract_x\n_atom_site_fract_y\n_atom_site_fract_z\n\
                _atom_site_occupancy\n_atom_site_disorder_group\n\
                _atom_site_disorder_assembly\n\
                C1 C 0 0 0 0.5 1 A\nO1 O 0.2 0 0 1 . .\n\
                loop_\n_atom_site_aniso_label\n_atom_site_aniso_U_11\n\
                _atom_site_aniso_U_22\n_atom_site_aniso_U_33\n\
                _atom_site_aniso_U_12\n_atom_site_aniso_U_13\n\
                _atom_site_aniso_U_23\nC1 0.01 0.02 0.03 0 0 0\n";
    let structure = read_cif(text).unwrap();
    assert_eq!(structure.sites()[0].disorder_group(), Some(1));
    assert_eq!(
        structure.sites()[0].props().get("disorder_assembly"),
        Some(&serde_json::Value::String("A".to_owned()))
    );
    assert!(structure.sites()[0].u_aniso().is_some());
    assert!(structure.sites()[1].u_aniso().is_none());
    let written = write_cif(&structure).unwrap();
    assert!(written.contains("_atom_site_aniso_label"));
    assert!(written.contains("_atom_site_disorder_group"));
    assert!(written.contains("_atom_site_disorder_assembly"));
    let back = read_cif(&written).unwrap();
    assert_eq!(back.sites()[0].disorder_group(), Some(1));
    assert_eq!(
        back.sites()[0].u_aniso(),
        Some([0.01, 0.02, 0.03, 0.0, 0.0, 0.0])
    );
    assert!(back.sites()[1].u_aniso().is_none());
}

#[test]
fn oxidation_number_fills_a_missing_site_charge() {
    let text = "data_oxidation\n_cell_length_a 10\n_cell_length_b 10\n_cell_length_c 10\n\
                _cell_angle_alpha 90\n_cell_angle_beta 90\n_cell_angle_gamma 90\n\
                loop_\n_atom_type_symbol\n_atom_type_oxidation_number\nFe3+ 3\n\
                loop_\n_atom_site_label\n_atom_site_type_symbol\n\
                _atom_site_fract_x\n_atom_site_fract_y\n_atom_site_fract_z\n\
                Fe1 Fe3+ 0 0 0\n";
    let structure = read_cif(text).unwrap();
    assert_eq!(structure.sites()[0].charge(), Some(3.0));
}

#[test]
fn round_trips_the_explicit_bond_loop() {
    let structure = StructureBuilder::new("explicit", Lattice::cubic(10.0).unwrap())
        .site(Site::new("C1", "C", [0.0, 0.0, 0.0]).unwrap())
        .site(Site::new("C2", "C", [0.1, 0.0, 0.0]).unwrap())
        .explicit_bond(ExplicitBond {
            i: 0,
            j: 1,
            image: [0, 0, 0],
            order: 2.0,
            remove: false,
        })
        .build()
        .unwrap();
    let written = write_cif(&structure).unwrap();
    assert!(written.contains("_geom_bond_atom_site_label_1"));
    assert!(written.contains("_xtal_bond_kind"));
    let back = read_cif(&written).unwrap();
    let bonds = back.explicit_bonds();
    assert_eq!(bonds.len(), 1);
    assert_eq!((bonds[0].i, bonds[0].j), (0, 1));
    assert_eq!(bonds[0].image, [0, 0, 0]);
    assert_eq!(bonds[0].order, 2.0);
    assert!(!bonds[0].remove);
}

#[test]
fn non_p1_explicit_bonds_map_to_the_golden_count() {
    let structure = read_cif(include_str!("../../../tests/corpus/MIL53.cif")).unwrap();
    let golden: serde_json::Value =
        serde_json::from_str(include_str!("../../../tests/golden/phase-1/MIL53.json")).unwrap();
    let stated = golden["bonds"]
        .as_array()
        .unwrap()
        .iter()
        .filter(|bond| bond["stated"].as_bool().unwrap_or(false))
        .count();
    assert_eq!(structure.explicit_bonds().len(), stated);
    assert_eq!(structure.explicit_bonds().len(), 36);
    assert!(
        structure
            .explicit_bonds()
            .iter()
            .all(|bond| (bond.order - 1.5).abs() < 1e-12)
    );
    let graph = structure.bond_graph().unwrap();
    assert_eq!(graph.bonds.len(), golden["bonds"].as_array().unwrap().len());
}

#[test]
fn an_operation_loop_without_a_space_group_symbol_is_not_p1() {
    let text = "data_test\n_cell_length_a 5\n_cell_length_b 5\n_cell_length_c 5\n\
                _cell_angle_alpha 90\n_cell_angle_beta 90\n_cell_angle_gamma 90\n\
                loop_\n_symmetry_equiv_pos_as_xyz\nx,y,z\n-x,-y,z\n\
                loop_\n_atom_site_label\n_atom_site_type_symbol\n\
                _atom_site_fract_x\n_atom_site_fract_y\n_atom_site_fract_z\n\
                Si1 Si 0.1 0.2 0.3\n";
    let structure = read_cif(text).unwrap();
    assert!(!structure.is_p1());
    assert_eq!(structure.p1().unwrap().n_atoms(), 2);
    assert_eq!(structure.formula(), "Si2");
}

#[test]
fn a_hall_symbol_without_an_operation_loop_expands_silicon() {
    let text = "data_Si\n_cell_length_a 5.430941\n_cell_length_b 5.430941\n\
                _cell_length_c 5.430941\n_cell_angle_alpha 90\n\
                _cell_angle_beta 90\n_cell_angle_gamma 90\n\
                _space_group_name_Hall 'F 4d 2 3 -1d'\n\
                loop_\n_atom_site_label\n_atom_site_type_symbol\n\
                _atom_site_fract_x\n_atom_site_fract_y\n_atom_site_fract_z\n\
                Si1 Si 0 0 0\n";
    let structure = read_cif(text).unwrap();
    assert!(!structure.is_p1());
    assert_eq!(structure.p1().unwrap().n_atoms(), 8);
    assert_eq!(structure.formula(), "Si8");
    assert!((structure.density().unwrap() - 2.329_142).abs() < 1e-5);
}

#[test]
fn an_operation_loop_count_disagreement_warns() {
    let text = "data_warn\n_cell_length_a 4\n_cell_length_b 4\n_cell_length_c 4\n\
                _cell_angle_alpha 90\n_cell_angle_beta 90\n_cell_angle_gamma 90\n\
                _symmetry_space_group_name_H-M 'P 21/c'\n\
                loop_\n_symmetry_equiv_pos_as_xyz\nx,y,z\n\
                loop_\n_atom_site_label\n_atom_site_type_symbol\n\
                _atom_site_fract_x\n_atom_site_fract_y\n_atom_site_fract_z\n\
                C1 C 0 0 0\n";
    let structure = read_cif(text).unwrap();
    let warning = structure
        .metadata()
        .get("symmetry_warning")
        .and_then(|value| value.as_str())
        .expect("a mismatch warning");
    assert!(warning.contains("operation loop"));
    assert!(warning.contains("P 21/c"));
}

#[test]
fn non_p1_explicit_bonds_round_trip_as_p1_indices() {
    let structure = read_cif(include_str!("../../../tests/corpus/MIL53.cif")).unwrap();
    assert_eq!(structure.explicit_bonds().len(), 36);
    let written = write_cif(&structure).unwrap();
    assert!(written.contains("_xtal_bond_p1_count 104"));
    assert!(written.contains("_xtal_bond_atom_1"));
    assert!(!written.contains("_geom_bond_atom_site_label_1"));
    let back = read_cif(&written).unwrap();
    assert_eq!(back.p1().unwrap().n_atoms(), 104);
    assert_eq!(back.explicit_bonds().len(), 36);
    assert_eq!(back.bond_graph().unwrap().bonds.len(), 126);
}

#[test]
fn a_mismatched_operation_set_warns() {
    let text = "data_warn\n_cell_length_a 4\n_cell_length_b 4\n_cell_length_c 4\n\
                _cell_angle_alpha 90\n_cell_angle_beta 90\n_cell_angle_gamma 90\n\
                _symmetry_space_group_name_H-M 'P 21/c'\n\
                loop_\n_symmetry_equiv_pos_as_xyz\nx,y,z\nx,y,z\nx,y,z\nx,y,z\n\
                loop_\n_atom_site_label\n_atom_site_type_symbol\n\
                _atom_site_fract_x\n_atom_site_fract_y\n_atom_site_fract_z\n\
                C1 C 0 0 0\n";
    let structure = read_cif(text).unwrap();
    let warning = structure
        .metadata()
        .get("symmetry_warning")
        .and_then(|value| value.as_str())
        .expect("a mismatch warning");
    assert!(warning.contains("does not match"));
}

#[test]
fn chemical_formula_is_kept_for_the_disorder_planner() {
    let text = "data_formula\n_chemical_formula_sum 'C1 N1'\n\
                _cell_length_a 10\n_cell_length_b 10\n_cell_length_c 10\n\
                _cell_angle_alpha 90\n_cell_angle_beta 90\n_cell_angle_gamma 90\n\
                loop_\n_atom_site_label\n_atom_site_type_symbol\n\
                _atom_site_fract_x\n_atom_site_fract_y\n_atom_site_fract_z\n\
                C1 C 0 0 0\n";
    let structure = read_cif(text).unwrap();
    assert_eq!(
        structure
            .metadata()
            .get("chemical_formula")
            .and_then(|value| value.as_str()),
        Some("C1 N1")
    );
}
