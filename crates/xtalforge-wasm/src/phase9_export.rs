#![allow(clippy::indexing_slicing)] // fixed three-coordinate chunks

//! WebAssembly exports for comparison, edits and local environments
//! (Tasks 9.1-9.3).

use serde::de::DeserializeOwned;
use serde::{Deserialize, Serialize};
use wasm_bindgen::prelude::*;

use xtalforge::core::Structure;

use crate::common::{parse_any, to_js_error};

/// Parses a structure from text with an optional format hint.
pub(crate) fn load(text: &str, hint: &str) -> Result<Structure, JsValue> {
    if hint.trim().is_empty() {
        parse_any(text).map_err(to_js_error)
    } else {
        xtalforge::io::read_structure(text, hint).map_err(to_js_error)
    }
}

/// Serializes any value to a JSON string.
pub(crate) fn json<T: Serialize>(value: &T) -> Result<String, JsValue> {
    serde_json::to_string(value)
        .map_err(|error| JsValue::from_str(&format!("could not serialize result: {error}")))
}

/// Parses an options JSON string, treating an empty string as defaults.
pub(crate) fn options<T: DeserializeOwned + Default>(text: &str) -> Result<T, JsValue> {
    if text.trim().is_empty() {
        Ok(T::default())
    } else {
        serde_json::from_str(text)
            .map_err(|error| JsValue::from_str(&format!("invalid options JSON: {error}")))
    }
}

/// A structure as its versioned JSON document.
pub(crate) fn document(structure: &Structure) -> Result<String, JsValue> {
    structure.to_document().to_json().map_err(to_js_error)
}

#[derive(Debug, Default, Deserialize)]
#[serde(default)]
pub(crate) struct MatcherOptions {
    pub ltol: f64,
    pub stol: f64,
    pub angle_tol: f64,
    pub allow_supercell: bool,
}

impl MatcherOptions {
    fn build(&self) -> xtalforge::symmetry::MatcherOptions {
        let defaults = xtalforge::symmetry::MatcherOptions::default();
        let mut options = defaults;
        if self.ltol > 0.0 {
            options.ltol = self.ltol;
        }
        if self.stol > 0.0 {
            options.stol = self.stol;
        }
        if self.angle_tol > 0.0 {
            options.angle_tol = self.angle_tol;
        }
        options.allow_supercell = self.allow_supercell;
        options
    }
}

/// Matches two structures and returns the match record as JSON.
///
/// # Errors
///
/// Returns a JavaScript `Error` for invalid text or tolerances.
#[wasm_bindgen]
pub fn match_structures(
    first: &str,
    second: &str,
    format_hint: &str,
    options_json: &str,
) -> Result<String, JsValue> {
    let options: MatcherOptions = options(options_json)?;
    let outcome = xtalforge::symmetry::match_structures(
        &load(first, format_hint)?,
        &load(second, format_hint)?,
        &options.build(),
    )
    .map_err(to_js_error)?;
    json(&outcome)
}

/// Periodic RMSD for an explicit atom mapping.
///
/// # Errors
///
/// Returns a JavaScript `Error` for invalid text or an incomplete mapping.
#[wasm_bindgen]
pub fn periodic_rmsd(
    first: &str,
    second: &str,
    mapping: Vec<usize>,
    format_hint: &str,
) -> Result<f64, JsValue> {
    xtalforge::core::periodic_rmsd(
        &load(first, format_hint)?,
        &load(second, format_hint)?,
        &mapping,
    )
    .map_err(to_js_error)
}

/// Kabsch alignment; points are a flat `xyzxyz...` array.
///
/// # Errors
///
/// Returns a JavaScript `Error` for mismatched or non-finite point sets.
#[wasm_bindgen]
pub fn kabsch_alignment(points: Vec<f64>, reference: Vec<f64>) -> Result<String, JsValue> {
    let points = flat_points(&points)?;
    let reference = flat_points(&reference)?;
    let alignment = xtalforge::core::kabsch_alignment(&points, &reference).map_err(to_js_error)?;
    json(&alignment)
}

fn flat_points(values: &[f64]) -> Result<Vec<[f64; 3]>, JsValue> {
    if !values.len().is_multiple_of(3) {
        return Err(JsValue::from_str("points need three coordinates per atom"));
    }
    Ok(values
        .chunks(3)
        .map(|chunk| [chunk[0], chunk[1], chunk[2]])
        .collect())
}

/// Equality classes of a list of structure texts.
///
/// # Errors
///
/// Returns a JavaScript `Error` for invalid text.
#[wasm_bindgen]
pub fn group_structures(
    structures: Vec<String>,
    format_hint: &str,
    options_json: &str,
) -> Result<String, JsValue> {
    let options: MatcherOptions = options(options_json)?;
    let mut parsed = Vec::with_capacity(structures.len());
    for text in &structures {
        parsed.push(load(text, format_hint)?);
    }
    let groups =
        xtalforge::symmetry::group_structures(&parsed, &options.build()).map_err(to_js_error)?;
    json(&serde_json::json!({
        "n_structures": parsed.len(),
        "n_groups": groups.len(),
        "groups": groups,
    }))
}

#[derive(Debug, Default, Deserialize)]
#[serde(default)]
pub(crate) struct DeformOptions {
    pub voigt: Option<[f64; 6]>,
    pub gradient: Option<[[f64; 3]; 3]>,
    pub volume: Option<f64>,
    pub density: Option<f64>,
}

/// Applies a strain, gradient, volume or density and returns structure JSON.
///
/// # Errors
///
/// Returns a JavaScript `Error` for invalid text or deformation options.
#[wasm_bindgen]
pub fn deform(text: &str, format_hint: &str, options_json: &str) -> Result<String, JsValue> {
    let options: DeformOptions = options(options_json)?;
    let mut structure = load(text, format_hint)?;
    if let Some(voigt) = options.voigt {
        let tensor = xtalforge::core::voigt_to_tensor(voigt).map_err(to_js_error)?;
        structure = xtalforge::core::strain(&structure, tensor).map_err(to_js_error)?;
    }
    if let Some(gradient) = options.gradient {
        structure = xtalforge::core::deform(&structure, gradient).map_err(to_js_error)?;
    }
    if let Some(volume) = options.volume {
        structure = xtalforge::core::scale_to_volume(&structure, volume).map_err(to_js_error)?;
    }
    if let Some(density) = options.density {
        structure = xtalforge::core::scale_to_density(&structure, density).map_err(to_js_error)?;
    }
    document(&structure)
}

#[derive(Debug, Default, Deserialize)]
#[serde(default)]
pub(crate) struct RattleOptions {
    pub amplitude: f64,
    pub seed: u64,
    pub bounded: bool,
}

/// Random displacements with a fixed seed; returns structure JSON.
///
/// # Errors
///
/// Returns a JavaScript `Error` for invalid text or options.
#[wasm_bindgen]
pub fn rattle(text: &str, format_hint: &str, options_json: &str) -> Result<String, JsValue> {
    let options: RattleOptions = options(options_json)?;
    let report = xtalforge::core::rattle(
        &load(text, format_hint)?,
        &xtalforge::core::RattleOptions {
            amplitude: if options.amplitude > 0.0 {
                options.amplitude
            } else {
                0.05
            },
            seed: options.seed,
            bounded: options.bounded,
            ..xtalforge::core::RattleOptions::default()
        },
    )
    .map_err(to_js_error)?;
    document(&report.structure)
}

#[derive(Debug, Default, Deserialize)]
#[serde(default)]
pub(crate) struct EditOptions {
    pub add: Vec<(String, [f64; 3])>,
    pub remove: Vec<usize>,
    pub replace: Vec<(usize, String)>,
    pub move_sites: Vec<(usize, [f64; 3])>,
    pub species: std::collections::BTreeMap<String, String>,
    pub fraction: f64,
    pub seed: u64,
    pub sort: Option<String>,
    pub wrap: bool,
}

/// Site edits; returns structure JSON.
///
/// # Errors
///
/// Returns a JavaScript `Error` for invalid text or edit options.
#[wasm_bindgen]
pub fn edit_structure(
    text: &str,
    format_hint: &str,
    options_json: &str,
) -> Result<String, JsValue> {
    let options: EditOptions = options(options_json)?;
    let mut structure = load(text, format_hint)?;
    if !options.species.is_empty() {
        structure = xtalforge::core::replace_species(
            &structure,
            &xtalforge::core::SpeciesReplacement {
                map: options.species,
                fraction: if options.fraction > 0.0 {
                    options.fraction
                } else {
                    1.0
                },
                seed: options.seed,
            },
        )
        .map_err(to_js_error)?;
    }
    if !options.add.is_empty() {
        structure = xtalforge::core::add_sites(&structure, &options.add).map_err(to_js_error)?;
    }
    if !options.remove.is_empty() {
        structure =
            xtalforge::core::remove_sites(&structure, &options.remove).map_err(to_js_error)?;
    }
    for (index, element) in options.replace {
        structure =
            xtalforge::core::replace_site(&structure, index, &element).map_err(to_js_error)?;
    }
    for (index, frac) in options.move_sites {
        structure = xtalforge::core::move_site(&structure, index, frac).map_err(to_js_error)?;
    }
    if let Some(sort) = options.sort {
        let key = match sort.as_str() {
            "element" => xtalforge::core::SiteSort::Element,
            "electronegativity" | "en" => xtalforge::core::SiteSort::Electronegativity,
            "mass" => xtalforge::core::SiteSort::Mass,
            "label" => xtalforge::core::SiteSort::Label,
            other => {
                return Err(JsValue::from_str(&format!("unknown sort key '{other}'")));
            }
        };
        structure = xtalforge::core::sort_sites(&structure, key).map_err(to_js_error)?;
    }
    if options.wrap {
        structure = xtalforge::core::wrap_sites(&structure).map_err(to_js_error)?;
    }
    document(&structure)
}

/// Linear interpolation between matched structures; returns a JSON array of
/// structure documents.
///
/// # Errors
///
/// Returns a JavaScript `Error` when the structures do not match.
#[wasm_bindgen]
pub fn interpolate(
    first: &str,
    second: &str,
    images: u32,
    format_hint: &str,
    options_json: &str,
) -> Result<String, JsValue> {
    let options: MatcherOptions = options(options_json)?;
    let frames = xtalforge::symmetry::interpolate_matched(
        &load(first, format_hint)?,
        &load(second, format_hint)?,
        images,
        &options.build(),
    )
    .map_err(to_js_error)?;
    let mut documents = Vec::with_capacity(frames.len());
    for frame in &frames {
        documents.push(frame.to_document().to_json().map_err(to_js_error)?);
    }
    json(&documents)
}

/// Supercell closest to a cube; returns structure JSON.
///
/// # Errors
///
/// Returns a JavaScript `Error` for an invalid cell count.
#[wasm_bindgen]
pub fn optimal_supercell(text: &str, cells: u32, format_hint: &str) -> Result<String, JsValue> {
    let choice = xtalforge::core::optimal_supercell(&load(text, format_hint)?, cells)
        .map_err(to_js_error)?;
    document(&choice.structure)
}
