//! Phase 19 exports: magnetism, disorder, diffraction and spectroscopy.

use wasm_bindgen::prelude::*;

use crate::common::{parse_any, to_js_error};

fn json(value: &impl serde::Serialize) -> Result<JsValue, JsValue> {
    let text = serde_json::to_string(value)
        .map_err(|error| JsValue::from_str(&format!("phase 19 JSON: {error}")))?;
    js_sys::JSON::parse(&text).map_err(|_| JsValue::from_str("phase 19 result is not JSON"))
}

/// Enumerates symmetry-ranked collinear magnetic orderings.
#[wasm_bindgen]
pub fn magnetic_orderings_info(
    text: &str,
    supercell: u32,
    max_orderings: usize,
) -> Result<JsValue, JsValue> {
    let structure = parse_any(text).map_err(to_js_error)?;
    let options = xtalforge::analysis::MagneticOrderingOptions {
        supercell: [supercell, supercell, supercell],
        max_orderings,
        ..xtalforge::analysis::MagneticOrderingOptions::default()
    };
    let result =
        xtalforge::analysis::magnetic_orderings(&structure, &options).map_err(to_js_error)?;
    json(&result)
}

/// Reports disorder groups, composition and density.
#[wasm_bindgen]
pub fn disorder_groups_info(text: &str) -> Result<JsValue, JsValue> {
    let structure = parse_any(text).map_err(to_js_error)?;
    let groups = xtalforge::analysis::disorder_groups(&structure).map_err(to_js_error)?;
    let composition = xtalforge::analysis::disorder_composition(&structure).map_err(to_js_error)?;
    let value = serde_json::json!({
        "formula": composition.formula(),
        "density": composition.density(structure.lattice().volume()).map_err(to_js_error)?,
        "groups": groups,
    });
    json(&value)
}

/// Enumerates symmetry-distinct derivative structures.
#[wasm_bindgen]
pub fn derivative_structures_info(
    text: &str,
    max_index: u32,
    max_structures: usize,
) -> Result<JsValue, JsValue> {
    let structure = parse_any(text).map_err(to_js_error)?;
    let options = xtalforge::analysis::DerivativeOptions {
        max_index,
        max_structures,
        ..xtalforge::analysis::DerivativeOptions::default()
    };
    let result =
        xtalforge::analysis::derivative_structures(&structure, &options).map_err(to_js_error)?;
    json(&result)
}

/// Calculates a nuclear neutron powder pattern.
#[wasm_bindgen]
pub fn neutron_powder_info(text: &str, wavelength: f64, step: f64) -> Result<JsValue, JsValue> {
    let structure = parse_any(text).map_err(to_js_error)?;
    let options = xtalforge::analysis::NeutronOptions {
        geometry: xtalforge::analysis::NeutronGeometry::ConstantWavelength { wavelength },
        step,
        ..xtalforge::analysis::NeutronOptions::default()
    };
    let pattern = xtalforge::analysis::neutron_powder(&structure, &options).map_err(to_js_error)?;
    json(&pattern)
}

/// Calculates the magnetic neutron powder pattern.
#[wasm_bindgen]
pub fn magnetic_powder_info(text: &str, wavelength: f64, step: f64) -> Result<JsValue, JsValue> {
    let structure = parse_any(text).map_err(to_js_error)?;
    let options = xtalforge::analysis::MagneticOptions {
        profile: xtalforge::analysis::InstrumentProfile {
            radiation: xtalforge::analysis::RadiationSource {
                kind: xtalforge::analysis::RadiationKind::Neutron,
                wavelength,
            },
            ..xtalforge::analysis::InstrumentProfile::default()
        },
        step,
        ..xtalforge::analysis::MagneticOptions::default()
    };
    let pattern = xtalforge::analysis::magnetic_diffraction::magnetic_powder(&structure, &options)
        .map_err(to_js_error)?;
    json(&pattern)
}

/// Calculates a kinematic SAED pattern for a zone axis.
#[wasm_bindgen]
pub fn electron_diffraction_info(
    text: &str,
    zone: &str,
    voltage: f64,
    size: usize,
) -> Result<JsValue, JsValue> {
    let structure = parse_any(text).map_err(to_js_error)?;
    let zone = parse_zone(zone)?;
    let options = xtalforge::analysis::SaedOptions {
        acceleration_kv: voltage,
        zone_axis: zone,
        image_size: size,
        ..xtalforge::analysis::SaedOptions::default()
    };
    let pattern =
        xtalforge::analysis::zone_axis_pattern(&structure, &options).map_err(to_js_error)?;
    let value = serde_json::json!({
        "spots": pattern.spots,
        "max_intensity": pattern.max_intensity(),
    });
    json(&value)
}

/// Returns the powder average of every 3x3 tensor in a JSON list.
#[wasm_bindgen]
pub fn raman_powder_average_info(tensors_json: &str) -> Result<JsValue, JsValue> {
    let tensors: Vec<[[f64; 3]; 3]> = serde_json::from_str(tensors_json)
        .map_err(|error| JsValue::from_str(&format!("tensors: {error}")))?;
    let values = tensors
        .iter()
        .map(xtalforge::ff::powder_average)
        .collect::<Vec<_>>();
    json(&values)
}

fn parse_zone(text: &str) -> Result<[i32; 3], JsValue> {
    let values = text
        .split(',')
        .map(|value| {
            value
                .trim()
                .parse::<i32>()
                .map_err(|_| JsValue::from_str(&format!("bad zone axis {text:?}")))
        })
        .collect::<Result<Vec<_>, _>>()?;
    match values.as_slice() {
        [u, v, w] => Ok([*u, *v, *w]),
        _ => Err(JsValue::from_str("zone axis needs three integers")),
    }
}
