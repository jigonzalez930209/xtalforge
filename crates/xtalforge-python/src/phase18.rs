//! Simulation-engine bindings (Phase 18).

use pyo3::prelude::*;
use pyo3::types::{PyDict, PyList};

use xtalforge::engines::{
    BandData, CalculationResult, ElectronicStructure, InputSetRequest, JobOptions, Scheduler,
    band_gap, band_structure, d_band_centre, failure_catalogue, list_presets, parse_cp2k,
    parse_gulp, parse_lammps, parse_lammps_dump, parse_qe, parse_qe_bands, parse_qe_projwfc,
    parse_qe_xml, parse_vasp_doscar, parse_vasp_eigenval, parse_vasp_oszicar, parse_vasp_outcar,
    parse_vasp_procar, parse_vasp_vasprun, preset_description, scheduler_script, smear_dos,
    write_input_set,
};

fn json_value(py: Python<'_>, value: &impl serde::Serialize) -> PyResult<Py<PyAny>> {
    let text = serde_json::to_string(value)
        .map_err(|error| crate::invalid_input(format!("engine JSON: {error}")))?;
    let module = py.import("json")?;
    module
        .call_method1("loads", (text,))
        .map(|value| value.unbind())
}

/// Lists every input-set preset with its description.
#[pyfunction]
pub(crate) fn engine_presets() -> PyResult<Vec<(String, String)>> {
    Ok(list_presets()
        .into_iter()
        .map(|name| {
            (
                name.to_owned(),
                preset_description(name).unwrap_or("").to_owned(),
            )
        })
        .collect())
}

/// Writes one input set and returns its named files.
#[pyfunction]
#[pyo3(signature = (text, preset, overrides = None, format_hint = None))]
pub(crate) fn engine_input_set<'py>(
    py: Python<'py>,
    text: &str,
    preset: &str,
    overrides: Option<Bound<'py, PyAny>>,
    format_hint: Option<&str>,
) -> PyResult<Bound<'py, PyDict>> {
    let structure = crate::structures::load_structure(text, format_hint)?;
    let overrides = match overrides {
        Some(value) => {
            let json = py.import("json")?;
            let text = json.call_method1("dumps", (value,))?.extract::<String>()?;
            serde_json::from_str(&text)
                .map_err(|error| crate::invalid_input(format!("overrides: {error}")))?
        }
        None => serde_json::Map::new(),
    };
    let summary = write_input_set(
        &structure,
        &InputSetRequest {
            preset: preset.to_owned(),
            overrides,
        },
    )
    .map_err(|error| crate::invalid_input(error.to_string()))?;
    let dict = PyDict::new(py);
    dict.set_item("preset", summary.preset)?;
    let files = PyList::empty(py);
    for file in &summary.files {
        let item = PyDict::new(py);
        item.set_item("name", &file.name)?;
        item.set_item("text", &file.text)?;
        files.append(item)?;
    }
    dict.set_item("files", files)?;
    Ok(dict)
}

/// Parses one engine output into the typed result model.
#[pyfunction]
pub(crate) fn engine_parse<'py>(
    py: Python<'py>,
    text: &str,
    code: &str,
) -> PyResult<Bound<'py, PyAny>> {
    let result = parse(code, text)?;
    Ok(json_value(py, &result)?.into_bound(py))
}

/// Parses a band file and returns the gap and effective masses.
#[pyfunction]
pub(crate) fn engine_bands<'py>(
    py: Python<'py>,
    text: &str,
    code: &str,
) -> PyResult<Bound<'py, PyAny>> {
    let result = parse(code, text)?;
    let bands = result
        .bands
        .ok_or_else(|| crate::invalid_input(format!("{code} has no band data")))?;
    let electronic = electronic_from_bands(&bands);
    let structure = band_structure(&electronic).map_err(crate::map_error)?;
    let gap = band_gap(&electronic).map_err(crate::map_error)?;
    let value = serde_json::json!({"bands": structure, "gap": gap});
    Ok(json_value(py, &value)?.into_bound(py))
}

/// Parses a DOS file, applies Gaussian smearing and returns the d-band centre.
#[pyfunction]
#[pyo3(signature = (text, code, sigma = 0.1))]
pub(crate) fn engine_dos<'py>(
    py: Python<'py>,
    text: &str,
    code: &str,
    sigma: f64,
) -> PyResult<Bound<'py, PyAny>> {
    let result = parse(code, text)?;
    let dos = result
        .dos
        .ok_or_else(|| crate::invalid_input(format!("{code} has no DOS data")))?;
    let smeared = smear_dos(&dos.energies, &dos.total, sigma).map_err(crate::map_error)?;
    let centre = d_band_centre(
        &dos.energies,
        dos.channels
            .iter()
            .find(|channel| channel.label == "d")
            .map_or(&dos.total, |channel| &channel.values),
    )
    .map_err(crate::map_error)?;
    let value = serde_json::json!({
        "energies": dos.energies,
        "total": dos.total,
        "smeared": smeared,
        "channels": dos.channels,
        "fermi": dos.fermi,
        "d_band_centre": centre,
    });
    Ok(json_value(py, &value)?.into_bound(py))
}

/// Writes a SLURM/PBS submission script.
#[pyfunction]
#[pyo3(signature = (scheduler, command, queue = "main", nodes = 1, walltime = "01:00:00"))]
pub(crate) fn engine_script(
    scheduler: &str,
    command: &str,
    queue: &str,
    nodes: usize,
    walltime: &str,
) -> PyResult<String> {
    let scheduler = parse_scheduler(scheduler)?;
    let options = JobOptions {
        queue: queue.to_owned(),
        nodes,
        walltime: walltime.to_owned(),
        ..JobOptions::new(scheduler, command)
    };
    scheduler_script(&options).map_err(|error| crate::invalid_input(error.to_string()))
}

/// Lists the known engine failure patterns.
#[pyfunction]
pub(crate) fn engine_failures() -> Vec<(String, String, String)> {
    failure_catalogue()
        .into_iter()
        .map(|pattern| (pattern.code, pattern.pattern, pattern.action))
        .collect()
}

fn parse(code: &str, text: &str) -> PyResult<CalculationResult> {
    let result = match code {
        "vasp-vasprun" | "vasprun" => parse_vasp_vasprun(text),
        "vasp-outcar" | "outcar" => parse_vasp_outcar(text),
        "vasp-oszicar" | "oszicar" => parse_vasp_oszicar(text),
        "vasp-eigenval" | "eigenval" => parse_vasp_eigenval(text),
        "vasp-doscar" | "doscar" => parse_vasp_doscar(text),
        "vasp-procar" | "procar" => parse_vasp_procar(text),
        "qe" | "pw" => parse_qe(text),
        "qe-xml" => parse_qe_xml(text),
        "qe-bands" => parse_qe_bands(text),
        "qe-projwfc" => parse_qe_projwfc(text),
        "cp2k" => parse_cp2k(text),
        "lammps" => parse_lammps(text),
        "lammps-dump" => parse_lammps_dump(text),
        "gulp" => parse_gulp(text),
        other => {
            return Err(crate::invalid_input(format!(
                "unknown engine output code {other:?}"
            )));
        }
    };
    result.map_err(|error| crate::invalid_input(error.to_string()))
}

fn parse_scheduler(name: &str) -> PyResult<Scheduler> {
    match name.to_ascii_lowercase().as_str() {
        "slurm" => Ok(Scheduler::Slurm),
        "pbs" => Ok(Scheduler::Pbs),
        other => Err(crate::invalid_input(format!(
            "unknown scheduler {other:?}; use slurm or pbs"
        ))),
    }
}

fn electronic_from_bands(bands: &BandData) -> ElectronicStructure {
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
    let mut structure = ElectronicStructure::new(bands.kpoints.clone(), energies);
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

/// Registers the Phase 18 functions on the extension module.
pub(crate) fn register(module: &Bound<'_, PyModule>) -> PyResult<()> {
    module.add_function(wrap_pyfunction!(engine_presets, module)?)?;
    module.add_function(wrap_pyfunction!(engine_input_set, module)?)?;
    module.add_function(wrap_pyfunction!(engine_parse, module)?)?;
    module.add_function(wrap_pyfunction!(engine_bands, module)?)?;
    module.add_function(wrap_pyfunction!(engine_dos, module)?)?;
    module.add_function(wrap_pyfunction!(engine_script, module)?)?;
    module.add_function(wrap_pyfunction!(engine_failures, module)?)?;
    Ok(())
}
