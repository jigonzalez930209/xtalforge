//! WASM force-field tests (Task 4.12.2).

#![cfg(target_arch = "wasm32")]
#![allow(
    clippy::unwrap_used,
    clippy::expect_used,
    clippy::panic,
    clippy::indexing_slicing
)]

use wasm_bindgen::JsCast;
use wasm_bindgen_test::*;
use xtalforge_wasm::{
    ff_charges_info, ff_energy, ff_ir_info, ff_optimize_info, ff_thermo_info, ff_types_info,
    ff_vibrations_info,
};

const SI_CIF: &str = include_str!("../../../tests/corpus/simple/Si.cif");
const NACL_CIF: &str = include_str!("../../../tests/corpus/simple/NaCl.cif");

fn property(value: &wasm_bindgen::JsValue, name: &str) -> wasm_bindgen::JsValue {
    js_sys::Reflect::get(value, &wasm_bindgen::JsValue::from_str(name)).expect("property")
}

#[wasm_bindgen_test]
fn force_field_exports() {
    let energy = ff_energy(SI_CIF).expect("energy");
    assert!(energy.is_finite());
    let types = ff_types_info(SI_CIF).expect("types");
    assert!(types.length() > 0);
    let vibrations = ff_vibrations_info(SI_CIF, 0.005).expect("vibrations");
    let frequencies = property(&vibrations, "frequencies");
    let frequencies: js_sys::Array = frequencies.unchecked_into();
    assert_eq!(frequencies.length(), 24);
    let ir = ff_ir_info(SI_CIF, 10.0).expect("ir");
    assert!(property(&ir, "peaks").is_object());
    let thermo = ff_thermo_info(SI_CIF, 300.0, 1.0).expect("thermo");
    assert!(
        property(&thermo, "gibbs")
            .as_f64()
            .unwrap_or(f64::NAN)
            .is_finite()
    );
    let optimized = ff_optimize_info(SI_CIF, "fire", 2, 0.05, false).expect("optimize");
    assert!(
        property(&optimized, "energy")
            .as_f64()
            .unwrap_or(f64::NAN)
            .is_finite()
    );
    let charges = ff_charges_info(NACL_CIF, "qeq", 0.0).expect("charges");
    assert!(property(&charges, "total").as_f64().unwrap_or(1.0).abs() < 1e-6);
    assert!(ff_charges_info(NACL_CIF, "bad", 0.0).is_err());
    assert!(ff_energy("not a structure").is_err());
}
