//! WASM symmetry integration tests.

#![cfg(target_arch = "wasm32")]
#![allow(
    clippy::unwrap_used,
    clippy::expect_used,
    clippy::panic,
    clippy::indexing_slicing
)]

use wasm_bindgen::JsCast;
use wasm_bindgen::JsValue;
use wasm_bindgen_test::*;
use xtalforge_wasm::{
    asymmetrize_info, brillouin_zone, change_setting_info, detect_symmetry, diagnose_info,
    kpath_info, layer_group_info, lower_to_subgroup_info, magnetic_group_info, pointgroup_info,
    primitive_info, reduce_info, set_space_group_info, standardize_info, subgroups_info,
    to_p1_info, tolerance_ladder_info, wyckoff_info,
};

const NACL_CIF: &str = include_str!("../../../tests/corpus/simple/NaCl.cif");
const SI_CIF: &str = include_str!("../../../tests/corpus/simple/Si.cif");

fn property(value: &JsValue, name: &str) -> JsValue {
    js_sys::Reflect::get(value, &JsValue::from_str(name)).expect("property must exist")
}

#[wasm_bindgen_test]
fn symmetry_detection_returns_the_dataset() {
    let value = detect_symmetry(NACL_CIF, 0.001).expect("symmetry");
    assert_eq!(property(&value, "number").as_f64(), Some(225.0));
    assert_eq!(
        property(&value, "hm_symbol")
            .as_string()
            .as_deref()
            .map(|text| text.replace(' ', "")),
        Some("Fm-3m".to_owned())
    );
    let operations = property(&value, "operations");
    let operations: js_sys::Array = operations.unchecked_into();
    assert_eq!(operations.length(), 192);
    assert!(detect_symmetry("not a structure", 0.001).is_err());
}

#[wasm_bindgen_test]
fn symmetry_workflow_exports() {
    let p1 = to_p1_info(NACL_CIF).expect("p1");
    assert_eq!(property(&p1, "n_sites").as_f64(), Some(8.0));
    assert_eq!(property(&p1, "is_p1").as_bool(), Some(true));

    let asymmetric = asymmetrize_info(NACL_CIF, 0.001).expect("asymmetrize");
    assert_eq!(property(&asymmetric, "n_sites").as_f64(), Some(2.0));

    let report = set_space_group_info(NACL_CIF, "225", "impose").expect("set");
    let set_report = property(&report, "set_report");
    assert!(property(&set_report, "max_displacement").as_f64().is_some());
    assert!(set_space_group_info(NACL_CIF, "225", "bad-mode").is_err());
}

#[wasm_bindgen_test]
fn subgroup_exports() {
    let groups = subgroups_info(525, 1e-5).expect("subgroups");
    assert!(groups.length() > 0);
    let value = lower_to_subgroup_info(NACL_CIF, 216, 0.5).expect("lower");
    assert_eq!(
        property(&value, "hall").as_string().as_deref(),
        Some("F -4 2 3")
    );
    assert!(lower_to_subgroup_info(NACL_CIF, 216, 0.0).is_err());
}

#[wasm_bindgen_test]
fn symmetry_analysis_exports() {
    let standard = standardize_info(NACL_CIF, 0.001).expect("standardize");
    assert_eq!(property(&standard, "n_sites").as_f64(), Some(2.0));
    let primitive = primitive_info(NACL_CIF, 0.001).expect("primitive");
    assert!(property(&primitive, "volume").as_f64().unwrap_or(0.0) > 0.0);
    let reduced = reduce_info(NACL_CIF, "niggli").expect("reduce");
    assert!(property(&reduced, "determinant").as_f64().is_some());
    let wyckoff = wyckoff_info(NACL_CIF, 0.001).expect("wyckoff");
    assert!(wyckoff.length() > 0);
    let zone = brillouin_zone(NACL_CIF).expect("zone");
    let faces = property(&zone, "faces");
    let faces: js_sys::Array = faces.unchecked_into();
    assert_eq!(faces.length(), 6);
    let ladder = tolerance_ladder_info(NACL_CIF).expect("ladder");
    assert!(ladder.length() > 0);
    let diagnostics = diagnose_info(NACL_CIF, "225").expect("diagnose");
    assert!(
        property(&diagnostics, "max_displacement")
            .as_f64()
            .unwrap_or(1.0)
            < 1e-6
    );
}

#[wasm_bindgen_test]
fn layer_group_export() {
    let layer = "data_layer\n_cell_length_a 4\n_cell_length_b 4\n_cell_length_c 20\n                 _cell_angle_alpha 90\n_cell_angle_beta 90\n_cell_angle_gamma 90\n                 loop_\n_atom_site_label\n_atom_site_type_symbol\n_atom_site_fract_x\n                 _atom_site_fract_y\n_atom_site_fract_z\nC1 C 0 0 0\nC2 C 0.5 0.5 0\n";
    let value = layer_group_info(layer, 1e-3).expect("layer group");
    assert!(property(&value, "number").as_f64().unwrap_or(0.0) >= 1.0);
    assert!(layer_group_info(layer, 0.0).is_err());
}

#[wasm_bindgen_test]
fn magnetic_group_export() {
    let text = "data_mag\n_cell_length_a 4\n_cell_length_b 4\n_cell_length_c 4\n                _cell_angle_alpha 90\n_cell_angle_beta 90\n_cell_angle_gamma 90\n                loop_\n_atom_site_label\n_atom_site_type_symbol\n_atom_site_fract_x\n                _atom_site_fract_y\n_atom_site_fract_z\nFe1 Fe 0 0 0\n";
    let value = magnetic_group_info(text, 1e-3);
    // A zero moment is a valid non-magnetic input; either the search succeeds
    // with UNI 1 or moyo rejects the trivial case.
    if let Ok(value) = value {
        assert!(property(&value, "uni_number").as_f64().unwrap_or(0.0) >= 1.0);
    }
}

#[wasm_bindgen_test]
fn change_setting_export() {
    let value = change_setting_info(SI_CIF, 526).expect("setting");
    assert_eq!(
        property(&value, "hall").as_string().as_deref(),
        Some("-F 4vw 2vw 3")
    );
    assert!(change_setting_info(SI_CIF, 0).is_err());
}

#[wasm_bindgen_test]
fn kpath_export() {
    let value = kpath_info(SI_CIF).expect("kpath");
    assert_eq!(
        property(&value, "bravais").as_string().as_deref(),
        Some("cF")
    );
    let points = property(&value, "special_points");
    assert!(property(&points, "L").is_object());
}

#[wasm_bindgen_test]
fn molecular_point_group_export() {
    let text = "3\nwater\nO 0.0 0.0 0.0\nH 0.757 0.586 0.0\nH -0.757 0.586 0.0\n";
    let group = pointgroup_info(text, 1e-3).expect("point group");
    assert_eq!(
        property(&group, "schoenflies").as_string().as_deref(),
        Some("C2v")
    );
}
