//! Simulation-engine exports (Phase 18).

use wasm_bindgen::prelude::*;

use crate::common::{parse_any, set, to_js_error};

fn object_from_json(text: &str) -> Result<JsValue, JsValue> {
    js_sys::JSON::parse(text).map_err(|_| JsValue::from_str("engine result is not JSON"))
}

/// Lists the input-set presets and descriptions.
///
/// # Errors
///
/// Never fails in the current build.
#[wasm_bindgen]
pub fn engine_presets_info() -> Result<JsValue, JsValue> {
    let array = js_sys::Array::new();
    for name in xtalforge::engines::list_presets() {
        let item = js_sys::Object::new();
        set(&item, "name", JsValue::from_str(name));
        set(
            &item,
            "description",
            JsValue::from_str(xtalforge::engines::preset_description(name).unwrap_or("")),
        );
        array.push(&item.into());
    }
    Ok(array.into())
}

/// Writes one input set and returns `{preset, files:[{name,text}]}`.
///
/// # Errors
///
/// Returns a JavaScript `Error` for an unknown preset, bad overrides or an
/// invalid structure.
#[wasm_bindgen]
pub fn engine_input_set_info(
    text: &str,
    format_hint: &str,
    preset: &str,
    overrides_json: &str,
) -> Result<JsValue, JsValue> {
    let structure = if format_hint.trim().is_empty() {
        parse_any(text).map_err(to_js_error)?
    } else {
        xtalforge::io::read_structure(text, format_hint).map_err(to_js_error)?
    };
    let overrides = serde_json::from_str(overrides_json)
        .map_err(|error| JsValue::from_str(&format!("overrides: {error}")))?;
    let summary = xtalforge::engines::write_input_set(
        &structure,
        &xtalforge::engines::InputSetRequest {
            preset: preset.to_owned(),
            overrides,
        },
    )
    .map_err(|error| JsValue::from_str(&error.to_string()))?;
    let text = serde_json::to_string(&summary)
        .map_err(|error| JsValue::from_str(&format!("input set: {error}")))?;
    object_from_json(&text)
}

/// Parses one engine output into the typed result model.
///
/// # Errors
///
/// Returns a JavaScript `Error` for an unknown code or malformed output.
#[wasm_bindgen]
pub fn engine_parse_info(text: &str, code: &str) -> Result<JsValue, JsValue> {
    let result = parse(code, text)?;
    let text = serde_json::to_string(&result)
        .map_err(|error| JsValue::from_str(&format!("engine result: {error}")))?;
    object_from_json(&text)
}

/// Parses a band file and returns the band structure and gap.
///
/// # Errors
///
/// Returns a JavaScript `Error` for an unknown code, missing bands or an
/// invalid electronic structure.
#[wasm_bindgen]
pub fn engine_bands_info(text: &str, code: &str) -> Result<JsValue, JsValue> {
    let result = parse(code, text)?;
    let bands = result
        .bands
        .ok_or_else(|| JsValue::from_str(&format!("{code} has no band data")))?;
    let electronic = electronic_from_bands(&bands);
    let structure = xtalforge::engines::band_structure(&electronic).map_err(to_js_error)?;
    let gap = xtalforge::engines::band_gap(&electronic).map_err(to_js_error)?;
    let text = serde_json::to_string(&serde_json::json!({"bands": structure, "gap": gap}))
        .map_err(|error| JsValue::from_str(&format!("bands: {error}")))?;
    object_from_json(&text)
}

/// Parses a DOS file, applies Gaussian smearing and returns the d-band centre.
///
/// # Errors
///
/// Returns a JavaScript `Error` for an unknown code, missing DOS data or an
/// invalid smearing width.
#[wasm_bindgen]
pub fn engine_dos_info(text: &str, code: &str, sigma: f64) -> Result<JsValue, JsValue> {
    let result = parse(code, text)?;
    let dos = result
        .dos
        .ok_or_else(|| JsValue::from_str(&format!("{code} has no DOS data")))?;
    let smeared =
        xtalforge::engines::smear_dos(&dos.energies, &dos.total, sigma).map_err(to_js_error)?;
    let centre = xtalforge::engines::d_band_centre(
        &dos.energies,
        dos.channels
            .iter()
            .find(|channel| channel.label == "d")
            .map_or(&dos.total, |channel| &channel.values),
    )
    .map_err(to_js_error)?;
    let text = serde_json::to_string(&serde_json::json!({
        "energies": dos.energies,
        "total": dos.total,
        "smeared": smeared,
        "channels": dos.channels,
        "fermi": dos.fermi,
        "d_band_centre": centre,
    }))
    .map_err(|error| JsValue::from_str(&format!("DOS: {error}")))?;
    object_from_json(&text)
}

/// Writes a SLURM/PBS submission script.
///
/// # Errors
///
/// Returns a JavaScript `Error` for an unknown scheduler or invalid options.
#[wasm_bindgen]
pub fn engine_script_info(
    scheduler: &str,
    command: &str,
    queue: &str,
    nodes: usize,
    walltime: &str,
) -> Result<String, JsValue> {
    let scheduler = match scheduler.to_ascii_lowercase().as_str() {
        "slurm" => xtalforge::engines::Scheduler::Slurm,
        "pbs" => xtalforge::engines::Scheduler::Pbs,
        other => {
            return Err(JsValue::from_str(&format!(
                "unknown scheduler {other:?}; use slurm or pbs"
            )));
        }
    };
    let options = xtalforge::engines::JobOptions {
        queue: queue.to_owned(),
        nodes,
        walltime: walltime.to_owned(),
        ..xtalforge::engines::JobOptions::new(scheduler, command)
    };
    xtalforge::engines::scheduler_script(&options)
        .map_err(|error| JsValue::from_str(&error.to_string()))
}

/// Lists the known failure patterns and corrective actions.
///
/// # Errors
///
/// Never fails in the current build.
#[wasm_bindgen]
pub fn engine_failures_info() -> Result<JsValue, JsValue> {
    let text = serde_json::to_string(&xtalforge::engines::failure_catalogue())
        .map_err(|error| JsValue::from_str(&format!("failures: {error}")))?;
    object_from_json(&text)
}

fn parse(code: &str, text: &str) -> Result<xtalforge::engines::CalculationResult, JsValue> {
    let result = match code {
        "vasp-vasprun" | "vasprun" => xtalforge::engines::parse_vasp_vasprun(text),
        "vasp-outcar" | "outcar" => xtalforge::engines::parse_vasp_outcar(text),
        "vasp-oszicar" | "oszicar" => xtalforge::engines::parse_vasp_oszicar(text),
        "vasp-eigenval" | "eigenval" => xtalforge::engines::parse_vasp_eigenval(text),
        "vasp-doscar" | "doscar" => xtalforge::engines::parse_vasp_doscar(text),
        "vasp-procar" | "procar" => xtalforge::engines::parse_vasp_procar(text),
        "qe" | "pw" => xtalforge::engines::parse_qe(text),
        "qe-xml" => xtalforge::engines::parse_qe_xml(text),
        "qe-bands" => xtalforge::engines::parse_qe_bands(text),
        "qe-projwfc" => xtalforge::engines::parse_qe_projwfc(text),
        "cp2k" => xtalforge::engines::parse_cp2k(text),
        "lammps" => xtalforge::engines::parse_lammps(text),
        "lammps-dump" => xtalforge::engines::parse_lammps_dump(text),
        "gulp" => xtalforge::engines::parse_gulp(text),
        other => {
            return Err(JsValue::from_str(&format!(
                "unknown engine output code {other:?}"
            )));
        }
    };
    result.map_err(|error| JsValue::from_str(&error.to_string()))
}

fn electronic_from_bands(
    bands: &xtalforge::engines::BandData,
) -> xtalforge::engines::ElectronicStructure {
    let n_bands = bands.energies.first().map_or(0, Vec::len);
    let mut energies = Vec::with_capacity(n_bands);
    for band in 0..n_bands {
        energies.push(
            bands
                .energies
                .iter()
                .map(|row| row.get(band).copied().unwrap_or(f64::NAN))
                .collect(),
        );
    }
    let mut structure =
        xtalforge::engines::ElectronicStructure::new(bands.kpoints.clone(), energies);
    if !bands.occupations.is_empty() {
        let mut occupations = Vec::with_capacity(n_bands);
        for band in 0..n_bands {
            occupations.push(
                bands
                    .occupations
                    .iter()
                    .map(|row| row.get(band).copied().unwrap_or(0.0))
                    .collect(),
            );
        }
        structure = structure.with_occupations(occupations);
    }
    if let Some(fermi) = bands.fermi {
        structure = structure.with_fermi(fermi);
    }
    structure
}
