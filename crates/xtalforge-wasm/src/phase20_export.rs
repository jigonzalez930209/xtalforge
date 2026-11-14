//! Phase 20 exports: generation, prototypes and global search.

use std::collections::BTreeMap;

use wasm_bindgen::prelude::*;

use crate::common::{parse_any, to_js_error};

fn json(value: &impl serde::Serialize) -> Result<JsValue, JsValue> {
    let text = serde_json::to_string(value)
        .map_err(|error| JsValue::from_str(&format!("phase 20 JSON: {error}")))?;
    js_sys::JSON::parse(&text).map_err(|_| JsValue::from_str("phase 20 result is not JSON"))
}

fn composition(text: &str) -> Result<Vec<(String, u32)>, JsValue> {
    let map: BTreeMap<String, u32> = serde_json::from_str(text)
        .map_err(|error| JsValue::from_str(&format!("composition JSON: {error}")))?;
    if map.is_empty() {
        return Err(JsValue::from_str("the composition is empty"));
    }
    Ok(map.into_iter().collect())
}

fn groups(text: &str) -> Result<Vec<u32>, JsValue> {
    let values: Vec<u32> = serde_json::from_str(text)
        .map_err(|error| JsValue::from_str(&format!("space-group JSON: {error}")))?;
    if values.is_empty() {
        return Err(JsValue::from_str("no space groups were given"));
    }
    Ok(values)
}

fn functional_group(name: &str) -> Result<xtalforge::mof::FunctionalGroup, JsValue> {
    match name.to_ascii_lowercase().as_str() {
        "amino" | "nh2" => Ok(xtalforge::mof::FunctionalGroup::Amino),
        "hydroxyl" | "oh" => Ok(xtalforge::mof::FunctionalGroup::Hydroxyl),
        "methyl" | "ch3" => Ok(xtalforge::mof::FunctionalGroup::Methyl),
        "fluorine" | "f" => Ok(xtalforge::mof::FunctionalGroup::Fluorine),
        other => Err(JsValue::from_str(&format!(
            "unknown functional group `{other}`"
        ))),
    }
}

/// Generates a Wyckoff-constrained random crystal.
#[wasm_bindgen]
pub fn generate_crystal_info(
    composition_json: &str,
    space_group: u32,
    seed: f64,
    volume_per_atom: f64,
    min_distance: f64,
    attempts: usize,
) -> Result<JsValue, JsValue> {
    let options = xtalforge::generation::CrystalOptions {
        space_group,
        composition: composition(composition_json)?,
        volume_per_atom,
        min_distance,
        max_attempts: attempts,
        seed: seed as u64,
    };
    let crystal = xtalforge::generation::generate_crystal(&options).map_err(to_js_error)?;
    json(&crystal)
}

/// Generates a layer crystal.
#[wasm_bindgen]
pub fn generate_layer_crystal_info(
    composition_json: &str,
    layer_group: u32,
    seed: f64,
    area_per_atom: f64,
    min_distance: f64,
    attempts: usize,
) -> Result<JsValue, JsValue> {
    let options = xtalforge::generation::LayerCrystalOptions {
        layer_group,
        composition: composition(composition_json)?,
        area_per_atom,
        min_distance,
        max_attempts: attempts,
        seed: seed as u64,
    };
    let crystal = xtalforge::generation::generate_layer_crystal(&options).map_err(to_js_error)?;
    json(&crystal)
}

/// Generates a rod crystal.
#[wasm_bindgen]
pub fn generate_rod_crystal_info(
    composition_json: &str,
    rod_group: u32,
    seed: f64,
    repeat: f64,
    min_distance: f64,
    attempts: usize,
) -> Result<JsValue, JsValue> {
    let options = xtalforge::generation::RodCrystalOptions {
        rod_group,
        composition: composition(composition_json)?,
        repeat,
        min_distance,
        max_attempts: attempts,
        seed: seed as u64,
    };
    let crystal = xtalforge::generation::generate_rod_crystal(&options).map_err(to_js_error)?;
    json(&crystal)
}

/// Places a rigid molecule on Wyckoff positions.
#[wasm_bindgen]
pub fn generate_molecule_crystal_info(
    text: &str,
    space_group: u32,
    copies: u32,
    seed: f64,
    volume_per_atom: f64,
    min_distance: f64,
    attempts: usize,
) -> Result<JsValue, JsValue> {
    let molecule = parse_any(text).map_err(to_js_error)?;
    let options = xtalforge::generation::MolecularCrystalOptions {
        space_group,
        molecule,
        copies,
        volume_per_atom,
        min_distance,
        tolerance: 0.25,
        max_attempts: attempts,
        seed: seed as u64,
    };
    let crystal =
        xtalforge::generation::generate_molecular_crystal(&options).map_err(to_js_error)?;
    json(&crystal)
}

/// Lists the 80 layer groups.
#[wasm_bindgen]
pub fn layer_groups_info() -> Result<JsValue, JsValue> {
    let groups = xtalforge::symmetry::layer_groups().map_err(to_js_error)?;
    json(&groups)
}

/// Lists the prototype library.
#[wasm_bindgen]
pub fn prototype_library_info() -> Result<JsValue, JsValue> {
    let library = xtalforge::generation::prototype_library().map_err(to_js_error)?;
    json(&library)
}

/// Builds a prototype by label.
#[wasm_bindgen]
pub fn build_prototype_info(label: &str, scale: f64) -> Result<JsValue, JsValue> {
    let structure = xtalforge::generation::build_prototype(
        label,
        &xtalforge::generation::PrototypeOptions {
            scale,
            ..xtalforge::generation::PrototypeOptions::default()
        },
    )
    .map_err(to_js_error)?;
    json(&structure)
}

/// Identifies the prototype label of a structure.
#[wasm_bindgen]
pub fn identify_prototype_info(text: &str, stol: f64) -> Result<JsValue, JsValue> {
    let structure = parse_any(text).map_err(to_js_error)?;
    let identification = xtalforge::generation::identify_prototype(
        &structure,
        &xtalforge::generation::PrototypeOptions {
            matcher: xtalforge::symmetry::MatcherOptions {
                stol,
                ..Default::default()
            },
            ..xtalforge::generation::PrototypeOptions::default()
        },
    )
    .map_err(to_js_error)?;
    json(&identification)
}

/// Proposes chemically similar substitutions.
#[wasm_bindgen]
pub fn substitution_candidates_info(
    text: &str,
    radius_tolerance: f64,
    limit: usize,
    same_charge: bool,
) -> Result<JsValue, JsValue> {
    let structure = parse_any(text).map_err(to_js_error)?;
    let report = xtalforge::generation::substitution_report(
        &structure,
        &xtalforge::generation::SubstitutionOptions {
            radius_tolerance,
            same_charge,
            limit,
        },
    )
    .map_err(to_js_error)?;
    json(&report)
}

/// Runs an AIRSS-style random search.
#[wasm_bindgen]
pub fn random_search_info(
    composition_json: &str,
    groups_json: &str,
    structures: usize,
    seed: f64,
    relax: bool,
    cell: bool,
    steps: usize,
) -> Result<JsValue, JsValue> {
    let options = xtalforge::generation::RandomSearchOptions {
        composition: composition(composition_json)?,
        space_groups: groups(groups_json)?,
        n_structures: structures,
        seed: seed as u64,
        volume_per_atom: 20.0,
        min_distance: 1.0,
        max_attempts: 32,
        relax,
        relax_cell: cell,
        max_steps: steps,
        deduplicate: true,
        matcher: Default::default(),
        minimizer: Default::default(),
    };
    let report = xtalforge::generation::random_search(&options).map_err(to_js_error)?;
    json(&report)
}

/// Runs an evolutionary search.
#[wasm_bindgen]
pub fn evolutionary_search_info(
    composition_json: &str,
    groups_json: &str,
    population: usize,
    generations: usize,
    seed: f64,
    cell: bool,
    steps: usize,
) -> Result<JsValue, JsValue> {
    let options = xtalforge::generation::EvolutionaryOptions {
        population,
        generations,
        relax_cell: cell,
        max_steps: steps,
        ..xtalforge::generation::EvolutionaryOptions::new(
            composition(composition_json)?,
            groups(groups_json)?,
            seed as u64,
        )
    };
    let report = xtalforge::generation::evolutionary_search(&options).map_err(to_js_error)?;
    json(&report)
}

/// Emits a checkpointed search workflow document.
#[wasm_bindgen]
pub fn search_workflow_info(
    kind: &str,
    seed_path: &str,
    options_json: &str,
    checkpoint: Option<String>,
    output: &str,
    name: &str,
) -> Result<JsValue, JsValue> {
    let kind = match kind {
        "random" => xtalforge::generation::SearchKind::Random,
        "evolutionary" => xtalforge::generation::SearchKind::Evolutionary,
        other => {
            return Err(JsValue::from_str(&format!(
                "search kind `{other}` is not random or evolutionary"
            )));
        }
    };
    let options: serde_json::Value = serde_json::from_str(options_json)
        .map_err(|error| JsValue::from_str(&format!("search options: {error}")))?;
    let document = xtalforge::generation::search_workflow_document(
        name,
        kind,
        seed_path,
        options,
        checkpoint.as_deref(),
        output,
    );
    json(&document)
}

/// Enumerates hypothetical MOF combinations.
#[wasm_bindgen]
pub fn mof_dataset_info(
    topologies_json: &str,
    nodes_json: &str,
    linkers_json: &str,
    limit: usize,
    relax: bool,
    descriptors: bool,
    seed: f64,
) -> Result<JsValue, JsValue> {
    let parse_list = |text: &str| -> Result<Vec<String>, JsValue> {
        serde_json::from_str(text)
            .map_err(|error| JsValue::from_str(&format!("name list JSON: {error}")))
    };
    let report =
        xtalforge::generation::generate_mof_dataset(&xtalforge::generation::MofDatasetOptions {
            topologies: parse_list(topologies_json)?,
            nodes: parse_list(nodes_json)?,
            linkers: parse_list(linkers_json)?,
            limit,
            relax,
            max_steps: 100,
            descriptors,
            pore_spacing: 1.0,
            seed: seed as u64,
        })
        .map_err(to_js_error)?;
    json(&report)
}

/// Functionalises a framework over groups and fractions.
#[wasm_bindgen]
pub fn functionalize_sweep_info(
    text: &str,
    groups_json: &str,
    fractions_json: &str,
    seed: f64,
    descriptors: bool,
) -> Result<JsValue, JsValue> {
    let structure = parse_any(text).map_err(to_js_error)?;
    let names: Vec<String> = serde_json::from_str(groups_json)
        .map_err(|error| JsValue::from_str(&format!("group JSON: {error}")))?;
    let fractions: Vec<f64> = serde_json::from_str(fractions_json)
        .map_err(|error| JsValue::from_str(&format!("fraction JSON: {error}")))?;
    let groups: Vec<xtalforge::mof::FunctionalGroup> = names
        .iter()
        .map(|name| functional_group(name))
        .collect::<Result<_, _>>()?;
    let sweep = xtalforge::generation::functionalization_sweep(
        &structure,
        &groups,
        &fractions,
        seed as u64,
        descriptors,
    )
    .map_err(to_js_error)?;
    json(&sweep)
}
