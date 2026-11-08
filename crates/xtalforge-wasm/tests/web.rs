//! WASM integration tests, run with `wasm-pack test --node`.

#![cfg(target_arch = "wasm32")]
#![allow(
    clippy::unwrap_used,
    clippy::expect_used,
    clippy::panic,
    clippy::indexing_slicing
)]

use serde::Deserialize;
use wasm_bindgen::JsCast;
use wasm_bindgen::JsValue;
use wasm_bindgen_test::*;
use xtalforge_wasm::{
    angle, bonds, cell_volume, convert_structure, detect_symmetry, distance, element_info, groups,
    lattice_info, lattice_minimum_image, limits_bounds, limits_check, neighbors, read_cif,
    read_structure, rings, slab_info, slab_text, structure_info, supercell_info, supercell_text,
    torsion, version,
};

const NACL_CIF: &str = include_str!("../../../tests/corpus/simple/NaCl.cif");

/// A planar regular hexagon, one C-C bond of 1.4 A per side.
const HEXAGON_XYZ: &str = "6
hexagon
C 3.4 2.0 0.0
C 2.7 3.212435565 0.0
C 1.3 3.212435565 0.0
C 0.6 2.0 0.0
C 1.3 0.787564435 0.0
C 2.7 0.787564435 0.0
";

/// Methanol, with one hydroxyl group.
const METHANOL_XYZ: &str = "3
methanol
C 0.0 0.0 0.0
O 1.43 0.0 0.0
H 1.90 0.0 0.0
";

#[derive(Debug, Deserialize)]
struct Vectors {
    cell_volume: Vec<Case>,
}

#[derive(Debug, Deserialize)]
struct Case {
    a: f64,
    b: f64,
    c: f64,
    alpha: f64,
    beta: f64,
    gamma: f64,
    expected: f64,
    tolerance: f64,
}

#[derive(Debug, Deserialize)]
struct ElementVectors {
    elements: Vec<ElementCase>,
}

#[derive(Debug, Deserialize)]
struct ElementCase {
    input: String,
    symbol: String,
    z: u8,
}

#[wasm_bindgen_test]
fn version_is_not_empty() {
    let value = version();
    assert!(!value.is_empty());
    assert_eq!(value, xtalforge::core::version());
}

#[wasm_bindgen_test]
fn cell_volume_shared_vectors() {
    let vectors: Vectors =
        serde_json::from_str(include_str!("../../../tests/vectors/hello_world.json"))
            .expect("shared hello-world vectors must be valid JSON");
    for case in vectors.cell_volume {
        let got = cell_volume(case.a, case.b, case.c, case.alpha, case.beta, case.gamma)
            .expect("valid vector case");
        assert!((got - case.expected).abs() <= case.tolerance);
    }
}

#[wasm_bindgen_test]
fn invalid_cell_returns_an_error() {
    assert!(cell_volume(0.0, 1.0, 1.0, 90.0, 90.0, 90.0).is_err());
}

fn property(value: &JsValue, name: &str) -> JsValue {
    js_sys::Reflect::get(value, &JsValue::from_str(name)).expect("property must exist")
}

#[wasm_bindgen_test]
fn element_info_shared_vectors() {
    let vectors: ElementVectors =
        serde_json::from_str(include_str!("../../../tests/vectors/elements.json"))
            .expect("element vectors must be valid JSON");
    for case in vectors.elements {
        let value = element_info(&case.input).expect("valid element input");
        assert_eq!(
            property(&value, "symbol").as_string().as_deref(),
            Some(case.symbol.as_str())
        );
        assert_eq!(property(&value, "z").as_f64(), Some(f64::from(case.z)));
    }
}

#[wasm_bindgen_test]
fn invalid_element_returns_an_error() {
    assert!(element_info("Zz").is_err());
}

#[wasm_bindgen_test]
fn lattice_shared_vectors() {
    let vectors: serde_json::Value =
        serde_json::from_str(include_str!("../../../tests/vectors/lattice.json"))
            .expect("lattice vectors must be valid JSON");
    for case in vectors["cases"].as_array().expect("cases array") {
        let number = |name: &str| case[name].as_f64().expect("numeric field");
        let info = lattice_info(
            number("a"),
            number("b"),
            number("c"),
            number("alpha"),
            number("beta"),
            number("gamma"),
        )
        .expect("valid cell");
        let volume = property(&info, "volume").as_f64().expect("volume");
        assert!((volume - number("volume")).abs() <= 1e-9);

        let delta: Vec<f64> = case["minimum_image"]["delta"]
            .as_array()
            .expect("delta")
            .iter()
            .map(|value| value.as_f64().expect("number"))
            .collect();
        let result = lattice_minimum_image(
            number("a"),
            number("b"),
            number("c"),
            number("alpha"),
            number("beta"),
            number("gamma"),
            delta,
        )
        .expect("valid cell");
        let distance = property(&result, "distance").as_f64().expect("distance");
        assert!(
            (distance
                - case["minimum_image"]["distance"]
                    .as_f64()
                    .expect("expected"))
            .abs()
                <= 1e-9
        );
    }
}

#[wasm_bindgen_test]
fn structure_info_vector() {
    let text = include_str!("../../../tests/vectors/structure.json");
    let value = structure_info(text).expect("valid structure JSON");
    assert_eq!(
        property(&value, "name").as_string().as_deref(),
        Some("NaCl")
    );
    assert_eq!(
        property(&value, "formula").as_string().as_deref(),
        Some("Na Cl")
    );
    assert_eq!(property(&value, "n_sites").as_f64(), Some(2.0));
    assert!(property(&value, "is_p1").as_bool().unwrap_or(false));
}

#[wasm_bindgen_test]
fn invalid_lattice_returns_an_error() {
    assert!(lattice_info(0.0, 1.0, 1.0, 90.0, 90.0, 90.0).is_err());
}

#[wasm_bindgen_test]
fn limits_surface() {
    let bounds = limits_bounds("standard", "supercell").expect("bounds");
    assert_eq!(bounds.get(0).as_f64(), Some(50_000.0));
    assert_eq!(bounds.get(1).as_f64(), Some(200_000.0));
    let warn = limits_check("standard", "supercell", 80_000.0).expect("warn");
    assert_eq!(
        property(&warn, "level").as_string().as_deref(),
        Some("warn")
    );
    let refuse = limits_check("standard", "supercell", 300_000.0).expect("refuse");
    assert_eq!(
        property(&refuse, "level").as_string().as_deref(),
        Some("refuse")
    );
    assert!(property(&refuse, "refused").as_bool().unwrap_or(false));
}

fn as_array(value: &JsValue) -> js_sys::Array {
    value
        .clone()
        .dyn_into::<js_sys::Array>()
        .expect("the export must return a JavaScript array")
}

#[wasm_bindgen_test]
fn read_structure_summary_from_cif_corpus() {
    let reference = xtalforge::io::read_cif(NACL_CIF).expect("NaCl corpus must parse");
    let expected = reference.summary().expect("summary");
    let value = read_structure(NACL_CIF, "cif").expect("CIF must parse");
    assert_eq!(
        property(&value, "name").as_string().as_deref(),
        Some(expected.name.as_str())
    );
    assert_eq!(
        property(&value, "formula").as_string().as_deref(),
        Some(expected.formula.as_str())
    );
    assert_eq!(
        property(&value, "hall").as_string().as_deref(),
        Some(expected.hall.as_str())
    );
    assert_eq!(
        property(&value, "n_sites").as_f64(),
        Some(expected.n_sites as f64)
    );
    let volume = property(&value, "volume").as_f64().expect("volume");
    assert!((volume - expected.volume).abs() <= 1e-9);
    assert!((volume - 177.5).abs() < 1.0);
    assert!(!property(&value, "is_p1").as_bool().unwrap_or(true));
}

#[wasm_bindgen_test]
fn read_structure_empty_hint_detects_the_format() {
    let value = read_structure(HEXAGON_XYZ, "").expect("XYZ must be detected");
    assert_eq!(property(&value, "n_sites").as_f64(), Some(6.0));
    assert!(property(&value, "is_p1").as_bool().unwrap_or(false));
}

#[wasm_bindgen_test]
fn supercell_info_matches_the_core_api() {
    let reference = xtalforge::io::read_cif(NACL_CIF).expect("NaCl corpus must parse");
    let expected = reference
        .supercell(2, 1, 1)
        .expect("supercell")
        .summary()
        .expect("summary");
    let value = supercell_info(NACL_CIF, 2, 1, 1).expect("supercell summary");
    assert_eq!(
        property(&value, "n_sites").as_f64(),
        Some(expected.n_sites as f64)
    );
    assert_eq!(property(&value, "n_sites").as_f64(), Some(4.0));
    let volume = property(&value, "volume").as_f64().expect("volume");
    assert!((volume - expected.volume).abs() <= 1e-9);
    assert!(property(&value, "is_p1").as_bool().unwrap_or(false));
    let repeats = as_array(&property(&value, "repeats"));
    assert_eq!(repeats.get(0).as_f64(), Some(2.0));
    assert_eq!(repeats.get(2).as_f64(), Some(1.0));
}

#[wasm_bindgen_test]
fn slab_info_matches_the_core_api() {
    let reference = xtalforge::io::read_cif(NACL_CIF).expect("NaCl corpus must parse");
    let expected = xtalforge::core::make_slab(&reference, [1, 0, 0], 2, 10.0, 0.0).expect("slab");
    let value = slab_info(NACL_CIF, 1, 0, 0, 2, 10.0).expect("slab summary");
    assert_eq!(
        property(&value, "n_sites").as_f64(),
        Some(expected.structure.sites().len() as f64)
    );
    assert!(property(&value, "is_p1").as_bool().unwrap_or(false));
    let hkl = as_array(&property(&value, "hkl"));
    assert_eq!(hkl.get(0).as_f64(), Some(1.0));
    assert_eq!(hkl.get(1).as_f64(), Some(0.0));
    assert_eq!(property(&value, "layers").as_f64(), Some(2.0));
    assert_eq!(property(&value, "vacuum").as_f64(), Some(10.0));
    assert!(property(&value, "thickness").as_f64().unwrap_or(0.0) > 0.0);
}

#[wasm_bindgen_test]
fn bonds_from_cif_corpus_match_the_core_api() {
    let reference = xtalforge::io::read_cif(NACL_CIF).expect("NaCl corpus must parse");
    let expected = reference.bond_graph().expect("bond graph");
    let value = bonds(NACL_CIF).expect("bond list");
    let array = as_array(&value);
    assert_eq!(array.length() as usize, expected.bonds.len());
    assert!(array.length() > 0);
    let first = array.get(0);
    assert_eq!(
        property(&first, "i").as_f64(),
        expected.bonds.first().map(|bond| bond.i as f64)
    );
    let image = as_array(&property(&first, "image"));
    assert_eq!(image.length(), 3);
    assert!(
        property(&first, "distance")
            .as_f64()
            .is_some_and(|distance| distance > 0.0)
    );
    assert!(!property(&first, "explicit").as_bool().unwrap_or(true));
}

#[wasm_bindgen_test]
fn rings_find_the_hexagon() {
    let reference = xtalforge::io::read_xyz(HEXAGON_XYZ).expect("hexagon must parse");
    let expected = xtalforge::core::rings(&reference, 8).expect("rings");
    let value = rings(HEXAGON_XYZ, 8).expect("ring list");
    let array = as_array(&value);
    assert_eq!(array.length() as usize, expected.len());
    assert_eq!(array.length(), 1);
    let ring = array.get(0);
    assert_eq!(property(&ring, "size").as_f64(), Some(6.0));
    assert_eq!(as_array(&property(&ring, "atoms")).length(), 6);
}

#[wasm_bindgen_test]
fn groups_find_the_hydroxyl() {
    let reference = xtalforge::io::read_xyz(METHANOL_XYZ).expect("methanol must parse");
    let expected = xtalforge::core::find_groups(&reference).expect("groups");
    let value = groups(METHANOL_XYZ).expect("group list");
    let array = as_array(&value);
    assert_eq!(array.length() as usize, expected.len());
    assert!(array.length() >= 1);
    let group = array.get(0);
    assert_eq!(
        property(&group, "kind").as_string().as_deref(),
        Some("hydroxyl")
    );
    assert_eq!(as_array(&property(&group, "atoms")).length(), 2);
}

#[wasm_bindgen_test]
fn neighbors_from_cif_corpus_match_the_core_api() {
    let reference = xtalforge::io::read_cif(NACL_CIF).expect("NaCl corpus must parse");
    let expected = xtalforge::core::neighbors(&reference, 3.0).expect("neighbors");
    let value = neighbors(NACL_CIF, 3.0).expect("neighbor list");
    let array = as_array(&value);
    assert_eq!(array.length() as usize, expected.len());
    assert!(array.length() > 0);
    for (index, expected) in expected.iter().enumerate().take(16) {
        let item = array.get(index as u32);
        assert_eq!(property(&item, "i").as_f64(), Some(expected.i as f64));
        assert_eq!(property(&item, "j").as_f64(), Some(expected.j as f64));
        assert!(
            property(&item, "distance")
                .as_f64()
                .is_some_and(|distance| (distance - expected.distance).abs() <= 1e-12)
        );
    }
}

#[wasm_bindgen_test]
fn measures_match_the_core_api() {
    let reference = xtalforge::io::read_xyz(HEXAGON_XYZ).expect("hexagon must parse");
    let d = distance(HEXAGON_XYZ, 0, 1).expect("distance");
    let d_expected = xtalforge::core::distance(&reference, 0, 1).expect("core distance");
    assert!((d - d_expected).abs() <= 1e-12);
    assert!((d - 1.4).abs() <= 1e-9);

    let a = angle(HEXAGON_XYZ, 0, 1, 2).expect("angle");
    let a_expected = xtalforge::core::angle(&reference, 0, 1, 2).expect("core angle");
    assert!((a - a_expected).abs() <= 1e-12);
    assert!((a - 120.0).abs() <= 1e-6);

    let t = torsion(HEXAGON_XYZ, 0, 1, 2, 3).expect("torsion");
    let t_expected = xtalforge::core::torsion(&reference, 0, 1, 2, 3).expect("core torsion");
    assert!((t - t_expected).abs() <= 1e-12);
}

#[wasm_bindgen_test]
fn phase_one_exports_reject_invalid_input() {
    assert!(read_structure(NACL_CIF, "unknown-format").is_err());
    assert!(bonds("not a structure").is_err());
    assert!(rings(HEXAGON_XYZ, 2).is_err());
    assert!(neighbors(HEXAGON_XYZ, 0.0).is_err());
    assert!(distance(HEXAGON_XYZ, 99, 0).is_err());
    assert!(supercell_info(NACL_CIF, 0, 1, 1).is_err());
    assert!(slab_info(NACL_CIF, 0, 0, 0, 1, 10.0).is_err());
}

#[wasm_bindgen_test]
fn phase_one_writers_return_text() {
    let summary = read_cif(NACL_CIF).expect("read_cif");
    assert_eq!(
        property(&summary, "formula").as_string().as_deref(),
        Some("Na Cl")
    );
    let xyz = convert_structure(NACL_CIF, "cif", "xyz").expect("convert");
    assert!(xyz.contains("Lattice="));
    let cell = supercell_text(NACL_CIF, 2, 1, 1).expect("supercell text");
    assert!(cell.starts_with("data_"));
    let slab = slab_text(NACL_CIF, 1, 0, 0, 1, 10.0).expect("slab text");
    assert!(slab.starts_with("data_"));
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
