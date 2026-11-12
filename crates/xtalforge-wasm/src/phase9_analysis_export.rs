//! WebAssembly exports for surfaces, defects, interfaces and mechanics
//! (Tasks 9.4-9.6).

use serde::Deserialize;
use wasm_bindgen::prelude::*;

use crate::phase9_export::{document, json, load, options};

fn document_value(structure: &xtalforge::core::Structure) -> Result<serde_json::Value, JsValue> {
    let text = document(structure)?;
    serde_json::from_str(&text)
        .map_err(|error| JsValue::from_str(&format!("could not serialize structure: {error}")))
}

#[derive(Debug, Default, Deserialize)]
#[serde(default)]
pub(crate) struct MillerOptions {
    pub max_index: i32,
    pub symprec: f64,
}

/// Symmetry-unique Miller indices as JSON.
///
/// # Errors
///
/// Returns a JavaScript `Error` for invalid text or a bad index bound.
#[wasm_bindgen]
pub fn unique_miller_indices(
    text: &str,
    format_hint: &str,
    options_json: &str,
) -> Result<String, JsValue> {
    let options: MillerOptions = options(options_json)?;
    let indices = xtalforge::analysis::unique_miller_indices(
        &load(text, format_hint)?,
        if options.max_index > 0 {
            options.max_index
        } else {
            2
        },
        if options.symprec > 0.0 {
            options.symprec
        } else {
            1e-5
        },
    )
    .map_err(crate::common::to_js_error)?;
    json(&indices)
}

#[derive(Debug, Default, Deserialize)]
#[serde(default)]
pub(crate) struct TerminationOptions {
    pub h: i32,
    pub k: i32,
    pub l: i32,
    pub layers: u32,
    pub vacuum: f64,
    pub symmetric: bool,
}

/// Every termination of a surface as JSON with structure documents.
///
/// # Errors
///
/// Returns a JavaScript `Error` for invalid text or a degenerate plane.
#[wasm_bindgen]
pub fn slab_terminations(
    text: &str,
    format_hint: &str,
    options_json: &str,
) -> Result<String, JsValue> {
    let options: TerminationOptions = options(options_json)?;
    let report = xtalforge::analysis::slab_terminations(
        &load(text, format_hint)?,
        [options.h, options.k, options.l],
        if options.layers > 0 {
            options.layers
        } else {
            2
        },
        if options.vacuum > 0.0 {
            options.vacuum
        } else {
            10.0
        },
        options.symmetric,
    )
    .map_err(crate::common::to_js_error)?;
    let mut entries = Vec::with_capacity(report.terminations.len());
    for termination in &report.terminations {
        entries.push(serde_json::json!({
            "index": termination.index,
            "shift": termination.shift,
            "thickness": termination.thickness,
            "dipole": termination.dipole,
            "polar": termination.polar,
            "structure": document_value(&termination.structure)?,
        }));
    }
    json(&serde_json::json!({
        "hkl": report.hkl,
        "n_terminations": report.n_terminations,
        "terminations": entries,
    }))
}

#[derive(Debug, Default, Deserialize)]
#[serde(default)]
pub(crate) struct WulffOptions {
    pub surfaces: Vec<(i32, i32, i32, f64)>,
}

/// Wulff construction as JSON.
///
/// # Errors
///
/// Returns a JavaScript `Error` for invalid text or surface energies.
#[wasm_bindgen]
pub fn wulff(text: &str, format_hint: &str, options_json: &str) -> Result<String, JsValue> {
    let options: WulffOptions = options(options_json)?;
    let structure = load(text, format_hint)?;
    let energies: Vec<xtalforge::core::SurfaceEnergy> = options
        .surfaces
        .into_iter()
        .map(|(h, k, l, energy)| xtalforge::core::SurfaceEnergy {
            hkl: [h, k, l],
            energy,
        })
        .collect();
    let shape = xtalforge::core::wulff(structure.lattice(), &energies)
        .map_err(crate::common::to_js_error)?;
    json(&shape)
}

/// Adsorption sites as JSON.
///
/// # Errors
///
/// Returns a JavaScript `Error` for invalid text or tolerance.
#[wasm_bindgen]
pub fn adsorption_sites(text: &str, format_hint: &str, tolerance: f64) -> Result<String, JsValue> {
    let sites = xtalforge::core::adsorption_sites(
        &load(text, format_hint)?,
        if tolerance > 0.0 { tolerance } else { 0.05 },
    )
    .map_err(crate::common::to_js_error)?;
    json(&sites)
}

#[derive(Debug, Default, Deserialize)]
#[serde(default)]
pub(crate) struct PlaceOptions {
    pub site_index: usize,
    pub offset: f64,
}

/// Places a molecule on an adsorption site and returns structure JSON.
///
/// # Errors
///
/// Returns a JavaScript `Error` for invalid text or an out-of-range site.
#[wasm_bindgen]
pub fn place_molecule(
    slab: &str,
    molecule: &str,
    format_hint: &str,
    options_json: &str,
) -> Result<String, JsValue> {
    let options: PlaceOptions = options(options_json)?;
    let slab = load(slab, format_hint)?;
    let molecule = load(molecule, format_hint)?;
    let sites =
        xtalforge::core::adsorption_sites(&slab, 0.05).map_err(crate::common::to_js_error)?;
    let site = sites
        .get(options.site_index)
        .ok_or_else(|| JsValue::from_str("adsorption site index is out of range"))?;
    let placed = xtalforge::core::place_molecule(
        &slab,
        &molecule,
        site,
        if options.offset != 0.0 {
            options.offset
        } else {
            2.0
        },
    )
    .map_err(crate::common::to_js_error)?;
    document(&placed)
}

#[derive(Debug, Default, Deserialize)]
#[serde(default)]
pub(crate) struct DefectOptions {
    pub species: String,
    pub min_image: f64,
    pub cutoff: f64,
    pub clearance: f64,
}

fn defects_to_json(defects: &[xtalforge::analysis::PointDefect]) -> Result<String, JsValue> {
    let mut entries = Vec::with_capacity(defects.len());
    for defect in defects {
        entries.push(serde_json::json!({
            "kind": format!("{:?}", defect.kind),
            "site_index": defect.site_index,
            "original": defect.original,
            "species": defect.species,
            "frac": defect.frac,
            "supercell": defect.supercell,
            "image_distance": defect.image_distance,
            "structure": document_value(&defect.structure)?,
        }));
    }
    json(&serde_json::json!({ "n_defects": defects.len(), "defects": entries }))
}

/// Symmetry-unique vacancies as JSON.
///
/// # Errors
///
/// Returns a JavaScript `Error` for invalid text or distance.
#[wasm_bindgen]
pub fn vacancies(text: &str, format_hint: &str, min_image: f64) -> Result<String, JsValue> {
    let defects = xtalforge::analysis::vacancies(
        &load(text, format_hint)?,
        if min_image > 0.0 { min_image } else { 10.0 },
    )
    .map_err(crate::common::to_js_error)?;
    defects_to_json(&defects)
}

/// Symmetry-unique substitutions as JSON.
///
/// # Errors
///
/// Returns a JavaScript `Error` for invalid text or species.
#[wasm_bindgen]
pub fn substitutions(
    text: &str,
    format_hint: &str,
    species: &str,
    min_image: f64,
) -> Result<String, JsValue> {
    let defects = xtalforge::analysis::substitutions(
        &load(text, format_hint)?,
        species,
        if min_image > 0.0 { min_image } else { 10.0 },
    )
    .map_err(crate::common::to_js_error)?;
    defects_to_json(&defects)
}

/// Interstitial candidates as JSON.
///
/// # Errors
///
/// Returns a JavaScript `Error` for invalid text or options.
#[wasm_bindgen]
pub fn interstitials(
    text: &str,
    format_hint: &str,
    species: &str,
    options_json: &str,
) -> Result<String, JsValue> {
    let options: DefectOptions = options(options_json)?;
    let defects = xtalforge::analysis::interstitials(
        &load(text, format_hint)?,
        species,
        if options.min_image > 0.0 {
            options.min_image
        } else {
            10.0
        },
        if options.cutoff > 0.0 {
            options.cutoff
        } else {
            5.0
        },
        if options.clearance > 0.0 {
            options.clearance
        } else {
            1.0
        },
    )
    .map_err(crate::common::to_js_error)?;
    defects_to_json(&defects)
}

#[derive(Debug, Default, Deserialize)]
#[serde(default)]
pub(crate) struct OrderingOptions {
    pub supercell: [u32; 3],
    pub max_orderings: usize,
    pub seed: u64,
    pub cutoff: f64,
    pub iterations: usize,
}

/// Orderings of partial occupancies as JSON.
///
/// # Errors
///
/// Returns a JavaScript `Error` for a structure without partial occupancy.
#[wasm_bindgen]
pub fn orderings(text: &str, format_hint: &str, options_json: &str) -> Result<String, JsValue> {
    let options: OrderingOptions = options(options_json)?;
    let results = xtalforge::analysis::orderings(
        &load(text, format_hint)?,
        &xtalforge::analysis::OrderingOptions {
            supercell: if options.supercell.contains(&0) {
                [1, 1, 1]
            } else {
                options.supercell
            },
            max_orderings: if options.max_orderings > 0 {
                options.max_orderings
            } else {
                64
            },
            ..xtalforge::analysis::OrderingOptions::default()
        },
    )
    .map_err(crate::common::to_js_error)?;
    let mut entries = Vec::with_capacity(results.len());
    for ordering in &results {
        entries.push(serde_json::json!({
            "energy": ordering.energy,
            "assignment": ordering.assignment,
            "structure": document_value(&ordering.structure)?,
        }));
    }
    json(&serde_json::json!({
        "n_orderings": results.len(),
        "orderings": entries,
    }))
}

/// SQS generation as JSON.
///
/// # Errors
///
/// Returns a JavaScript `Error` for a structure without binary disorder.
#[wasm_bindgen]
pub fn special_quasirandom(
    text: &str,
    format_hint: &str,
    options_json: &str,
) -> Result<String, JsValue> {
    let options: OrderingOptions = options(options_json)?;
    let result = xtalforge::analysis::special_quasirandom(
        &load(text, format_hint)?,
        &xtalforge::analysis::SqsOptions {
            supercell: if options.supercell.contains(&0) {
                [2, 2, 2]
            } else {
                options.supercell
            },
            seed: options.seed,
            cutoff: if options.cutoff > 0.0 {
                options.cutoff
            } else {
                6.0
            },
            iterations: if options.iterations > 0 {
                options.iterations
            } else {
                5000
            },
            ..xtalforge::analysis::SqsOptions::default()
        },
    )
    .map_err(crate::common::to_js_error)?;
    json(&serde_json::json!({
        "assignment": result.assignment,
        "shells": result.shells,
        "triplet_correlation": result.triplet_correlation,
        "objective": result.objective,
        "accepted": result.accepted,
        "attempted": result.attempted,
        "structure": document_value(&result.structure)?,
    }))
}

#[derive(Debug, Default, Deserialize)]
#[serde(default)]
pub(crate) struct InterfaceOptions {
    pub max_index: i32,
    pub max_area: f64,
    pub tolerance: f64,
    pub gap: f64,
    pub vacuum: f64,
    pub shift: [f64; 3],
}

/// Commensurate 2D lattice matches as JSON.
///
/// # Errors
///
/// Returns a JavaScript `Error` for invalid text or options.
#[wasm_bindgen]
pub fn lattice_matches(
    film: &str,
    substrate: &str,
    format_hint: &str,
    options_json: &str,
) -> Result<String, JsValue> {
    let options: InterfaceOptions = options(options_json)?;
    let matches = xtalforge::analysis::lattice_matches(
        &load(film, format_hint)?,
        &load(substrate, format_hint)?,
        &xtalforge::analysis::InterfaceOptions {
            max_index: if options.max_index > 0 {
                options.max_index
            } else {
                3
            },
            max_area: if options.max_area > 0.0 {
                options.max_area
            } else {
                200.0
            },
            length_tolerance: if options.tolerance > 0.0 {
                options.tolerance
            } else {
                0.05
            },
            ..xtalforge::analysis::InterfaceOptions::default()
        },
    )
    .map_err(crate::common::to_js_error)?;
    json(&matches)
}

/// Stacks a film on a substrate; returns structure JSON.
///
/// # Errors
///
/// Returns a JavaScript `Error` for empty slabs or invalid options.
#[wasm_bindgen]
pub fn stack_slabs(
    film: &str,
    substrate: &str,
    format_hint: &str,
    options_json: &str,
) -> Result<String, JsValue> {
    let options: InterfaceOptions = options(options_json)?;
    let result = xtalforge::analysis::stack_slabs(
        &load(film, format_hint)?,
        &load(substrate, format_hint)?,
        if options.gap >= 0.0 { options.gap } else { 2.0 },
        options.shift,
        if options.vacuum >= 0.0 {
            options.vacuum
        } else {
            15.0
        },
    )
    .map_err(crate::common::to_js_error)?;
    document(&result.structure)
}
