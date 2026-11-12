//! WebAssembly exports for the MOF platform (Phase 10).

use std::collections::BTreeMap;

use serde::Deserialize;
use wasm_bindgen::prelude::*;

use crate::common::to_js_error;
use crate::phase9_export::{document, json, load, options};

#[derive(Debug, Default, Deserialize)]
#[serde(default)]
struct SolventOptions {
    remove_ions: bool,
}

/// MOF quality report as JSON.
#[wasm_bindgen]
pub fn mof_check(text: &str, format_hint: &str) -> Result<String, JsValue> {
    json(&xtalforge::mof::check_mof(&load(text, format_hint)?).map_err(to_js_error)?)
}

/// Solvent removal report as JSON.
#[wasm_bindgen]
pub fn mof_solvent_remove(
    text: &str,
    format_hint: &str,
    options_json: &str,
) -> Result<String, JsValue> {
    let options: SolventOptions = options(options_json)?;
    let report = xtalforge::mof::remove_solvent(
        &load(text, format_hint)?,
        &xtalforge::mof::SolventOptions {
            remove_ions: options.remove_ions,
            ..xtalforge::mof::SolventOptions::default()
        },
    )
    .map_err(to_js_error)?;
    json(&report)
}

/// Disorder resolution; returns structure JSON.
#[wasm_bindgen]
pub fn mof_disorder(text: &str, format_hint: &str) -> Result<String, JsValue> {
    let structure = xtalforge::mof::resolve_disorder(
        &load(text, format_hint)?,
        &xtalforge::mof::DisorderOptions::default(),
    )
    .map_err(to_js_error)?;
    document(&structure)
}

/// Charge-balance report as JSON.
#[wasm_bindgen]
pub fn mof_charge_balance(text: &str, format_hint: &str) -> Result<String, JsValue> {
    json(&xtalforge::mof::charge_balance(&load(text, format_hint)?).map_err(to_js_error)?)
}

/// Deconstruction as JSON.
#[wasm_bindgen]
pub fn mof_deconstruct(text: &str, format_hint: &str) -> Result<String, JsValue> {
    json(&xtalforge::mof::deconstruct(&load(text, format_hint)?).map_err(to_js_error)?)
}

/// MOFid and MOFkey as JSON.
#[wasm_bindgen]
pub fn mof_mofid(text: &str, format_hint: &str) -> Result<String, JsValue> {
    json(&xtalforge::mof::mofid(&load(text, format_hint)?).map_err(to_js_error)?)
}

/// Open metal sites as JSON.
#[wasm_bindgen]
pub fn mof_open_metal_sites(text: &str, format_hint: &str) -> Result<String, JsValue> {
    json(
        &xtalforge::mof::open_metal_sites(
            &load(text, format_hint)?,
            &xtalforge::mof::OpenMetalOptions::default(),
        )
        .map_err(to_js_error)?,
    )
}

/// Framework descriptors as JSON.
#[wasm_bindgen]
pub fn mof_descriptors(text: &str, format_hint: &str) -> Result<String, JsValue> {
    json(&xtalforge::mof::framework_descriptors(&load(text, format_hint)?).map_err(to_js_error)?)
}

#[derive(Debug, Default, Deserialize)]
#[serde(default)]
struct PoreOptions {
    spacing: f64,
    probe: f64,
}

/// Pore descriptors as JSON.
#[wasm_bindgen]
pub fn mof_pores(text: &str, format_hint: &str, options_json: &str) -> Result<String, JsValue> {
    let options: PoreOptions = options(options_json)?;
    let report = xtalforge::mof::pore_descriptors(
        &load(text, format_hint)?,
        &xtalforge::mof::PoreOptions {
            spacing: if options.spacing > 0.0 {
                options.spacing
            } else {
                0.8
            },
            probe_radius: if options.probe > 0.0 {
                options.probe
            } else {
                1.2
            },
            channel_probe: if options.probe > 0.0 {
                options.probe
            } else {
                1.2
            },
        },
    )
    .map_err(to_js_error)?;
    json(&report)
}

/// Embedded guests as JSON.
#[wasm_bindgen]
pub fn mof_guests() -> Result<String, JsValue> {
    json(&xtalforge::mof::guests().map_err(to_js_error)?)
}

/// Host-guest force field as JSON.
#[wasm_bindgen]
pub fn mof_forcefield(text: &str, format_hint: &str) -> Result<String, JsValue> {
    json(&xtalforge::mof::force_field(&load(text, format_hint)?).map_err(to_js_error)?)
}

#[derive(Debug, Default, Deserialize)]
#[serde(default)]
struct WidomOptions {
    guest: String,
    temperature: f64,
    samples: usize,
    seed: u64,
}

/// Widom insertion as JSON.
#[wasm_bindgen]
pub fn mof_widom(text: &str, format_hint: &str, options_json: &str) -> Result<String, JsValue> {
    let options: WidomOptions = options(options_json)?;
    json(
        &xtalforge::mof::widom_insertion(
            &load(text, format_hint)?,
            if options.guest.is_empty() {
                "CO2"
            } else {
                &options.guest
            },
            &xtalforge::mof::WidomOptions {
                temperature: if options.temperature > 0.0 {
                    options.temperature
                } else {
                    298.0
                },
                samples: if options.samples > 0 {
                    options.samples
                } else {
                    1000
                },
                seed: options.seed,
            },
        )
        .map_err(to_js_error)?,
    )
}

#[derive(Debug, Default, Deserialize)]
#[serde(default)]
struct GcmcOptions {
    guest: String,
    temperature: f64,
    pressures: Vec<f64>,
    steps: usize,
    seed: u64,
}

/// GCMC isotherm as JSON.
#[wasm_bindgen]
pub fn mof_gcmc(text: &str, format_hint: &str, options_json: &str) -> Result<String, JsValue> {
    let options: GcmcOptions = options(options_json)?;
    let pressures = if options.pressures.is_empty() {
        vec![1.0, 5.0]
    } else {
        options.pressures
    };
    json(
        &xtalforge::mof::gcmc_isotherm(
            &load(text, format_hint)?,
            if options.guest.is_empty() {
                "CO2"
            } else {
                &options.guest
            },
            &xtalforge::mof::GcmcOptions {
                temperature: if options.temperature > 0.0 {
                    options.temperature
                } else {
                    298.0
                },
                pressures,
                steps: if options.steps > 0 {
                    options.steps
                } else {
                    1000
                },
                seed: options.seed,
            },
        )
        .map_err(to_js_error)?,
    )
}

#[derive(Debug, Default, Deserialize)]
#[serde(default)]
struct GridOptions {
    guest: String,
    spacing: f64,
}

/// Energy grid summary as JSON.
#[wasm_bindgen]
pub fn mof_energy_grid(
    text: &str,
    format_hint: &str,
    options_json: &str,
) -> Result<String, JsValue> {
    let options: GridOptions = options(options_json)?;
    let grid = xtalforge::mof::energy_grid(
        &load(text, format_hint)?,
        if options.guest.is_empty() {
            "CO2"
        } else {
            &options.guest
        },
        if options.spacing > 0.0 {
            options.spacing
        } else {
            1.0
        },
    )
    .map_err(to_js_error)?;
    json(&serde_json::json!({
        "guest": grid.guest,
        "shape": grid.grid.shape,
        "minimum": grid.minimum,
        "maximum": grid.maximum,
    }))
}

#[derive(Debug, Default, Deserialize)]
#[serde(default)]
struct RaspaOptions {
    guest: String,
    temperature: f64,
    pressure: f64,
}

/// RASPA input set as JSON.
#[wasm_bindgen]
pub fn mof_raspa(text: &str, format_hint: &str, options_json: &str) -> Result<String, JsValue> {
    let options: RaspaOptions = options(options_json)?;
    json(
        &xtalforge::mof::raspa_input(
            &load(text, format_hint)?,
            if options.guest.is_empty() {
                "CO2"
            } else {
                &options.guest
            },
            if options.temperature > 0.0 {
                options.temperature
            } else {
                298.0
            },
            if options.pressure > 0.0 {
                options.pressure
            } else {
                1.0
            },
        )
        .map_err(to_js_error)?,
    )
}

#[derive(Debug, Default, Deserialize)]
#[serde(default)]
struct FunctionalizeOptions {
    group: String,
    positions: Vec<usize>,
    fraction: f64,
    seed: u64,
}

/// Functionalisation; returns structure JSON.
#[wasm_bindgen]
pub fn mof_functionalize(
    text: &str,
    format_hint: &str,
    options_json: &str,
) -> Result<String, JsValue> {
    let options: FunctionalizeOptions = options(options_json)?;
    let group = match options.group.as_str() {
        "hydroxyl" | "oh" => xtalforge::mof::FunctionalGroup::Hydroxyl,
        "methyl" | "ch3" => xtalforge::mof::FunctionalGroup::Methyl,
        "fluoro" | "f" => xtalforge::mof::FunctionalGroup::Fluorine,
        _ => xtalforge::mof::FunctionalGroup::Amino,
    };
    let structure = xtalforge::mof::functionalize(
        &load(text, format_hint)?,
        &xtalforge::mof::FunctionalizeOptions {
            group,
            positions: options.positions,
            seed: options.seed,
            fraction: if options.fraction > 0.0 {
                options.fraction
            } else {
                1.0
            },
        },
    )
    .map_err(to_js_error)?;
    document(&structure)
}

#[derive(Debug, Default, Deserialize)]
#[serde(default)]
struct DefectOptions {
    linker_fraction: f64,
    cluster_fraction: f64,
    cap: String,
    seed: u64,
}

/// Missing-linker/cluster defects as JSON.
#[wasm_bindgen]
pub fn mof_defects(text: &str, format_hint: &str, options_json: &str) -> Result<String, JsValue> {
    let options: DefectOptions = options(options_json)?;
    let cap = match options.cap.as_str() {
        "water" => xtalforge::mof::CapSpecies::Water,
        "hydroxide" | "oh" => xtalforge::mof::CapSpecies::Hydroxide,
        _ => xtalforge::mof::CapSpecies::Formate,
    };
    let report = xtalforge::mof::missing_defects(
        &load(text, format_hint)?,
        &xtalforge::mof::DefectOptions {
            linker_fraction: options.linker_fraction,
            cluster_fraction: options.cluster_fraction,
            cap,
            seed: options.seed,
        },
    )
    .map_err(to_js_error)?;
    json(&serde_json::json!({
        "removed_linkers": report.removed_linkers,
        "removed_clusters": report.removed_clusters,
        "caps": report.caps,
        "structure": serde_json::from_str::<serde_json::Value>(&document(&report.structure)?)
            .map_err(|error| JsValue::from_str(&error.to_string()))?,
    }))
}

#[derive(Debug, Default, Deserialize)]
#[serde(default)]
struct VariantOptions {
    folds: u32,
    layers: u32,
    offset: [f64; 2],
    map: BTreeMap<String, String>,
    ordered: bool,
    fraction: f64,
    seed: u64,
}

/// n-fold interpenetration; returns structure JSON.
#[wasm_bindgen]
pub fn mof_interpenetrate(
    text: &str,
    format_hint: &str,
    options_json: &str,
) -> Result<String, JsValue> {
    let options: VariantOptions = options(options_json)?;
    let structure = xtalforge::mof::interpenetrate(
        &load(text, format_hint)?,
        if options.folds >= 2 { options.folds } else { 2 },
    )
    .map_err(to_js_error)?;
    document(&structure)
}

/// Layer stacking; returns structure JSON.
#[wasm_bindgen]
pub fn mof_stack_layers(
    text: &str,
    format_hint: &str,
    options_json: &str,
) -> Result<String, JsValue> {
    let options: VariantOptions = options(options_json)?;
    let structure = xtalforge::mof::stack_layers(
        &load(text, format_hint)?,
        if options.layers > 0 {
            options.layers
        } else {
            2
        },
        options.offset,
    )
    .map_err(to_js_error)?;
    document(&structure)
}

/// Mixed-metal/linker framework; returns structure JSON.
#[wasm_bindgen]
pub fn mof_mixed(text: &str, format_hint: &str, options_json: &str) -> Result<String, JsValue> {
    let options: VariantOptions = options(options_json)?;
    let structure = xtalforge::mof::mixed_framework(
        &load(text, format_hint)?,
        &xtalforge::mof::MixtureOptions {
            map: options.map,
            ordered: options.ordered,
            fraction: if options.fraction > 0.0 {
                options.fraction
            } else {
                1.0
            },
            seed: options.seed,
        },
    )
    .map_err(to_js_error)?;
    document(&structure)
}
