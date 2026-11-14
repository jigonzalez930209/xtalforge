//! Phase 20 bindings: generation, prototypes and global search.

use std::collections::BTreeMap;

use pyo3::prelude::*;

use xtalforge::generation::{
    CrystalOptions, EvolutionaryOptions, LayerCrystalOptions, MofDatasetOptions,
    MolecularCrystalOptions, PrototypeOptions, RandomSearchOptions, RodCrystalOptions, SearchKind,
    SubstitutionOptions, build_prototype, evolutionary_search, functionalization_sweep,
    generate_crystal, generate_crystal_range, generate_layer_crystal, generate_mof_dataset,
    generate_molecular_crystal, generate_rod_crystal, identify_prototype, prototype_library,
    random_search, search_workflow_document, substitution_report,
};

fn json_value(py: Python<'_>, value: &impl serde::Serialize) -> PyResult<Py<PyAny>> {
    let text = serde_json::to_string(value)
        .map_err(|error| crate::invalid_input(format!("phase 20 JSON: {error}")))?;
    let module = py.import("json")?;
    module
        .call_method1("loads", (text,))
        .map(|value| value.unbind())
}

fn composition_map(values: &BTreeMap<String, u32>) -> Vec<(String, u32)> {
    values
        .iter()
        .map(|(element, count)| (element.clone(), *count))
        .collect()
}

fn space_group_list(values: Vec<u32>) -> Vec<u32> {
    if values.is_empty() { vec![225] } else { values }
}

fn group(name: &str) -> PyResult<xtalforge::mof::FunctionalGroup> {
    match name.to_ascii_lowercase().as_str() {
        "amino" | "nh2" => Ok(xtalforge::mof::FunctionalGroup::Amino),
        "hydroxyl" | "oh" => Ok(xtalforge::mof::FunctionalGroup::Hydroxyl),
        "methyl" | "ch3" => Ok(xtalforge::mof::FunctionalGroup::Methyl),
        "fluorine" | "f" => Ok(xtalforge::mof::FunctionalGroup::Fluorine),
        other => Err(crate::invalid_input(format!(
            "unknown functional group `{other}`"
        ))),
    }
}

/// Generates a Wyckoff-constrained random crystal.
#[pyfunction(name = "generate_crystal")]
#[pyo3(signature = (composition, space_group = 225, seed = 0, volume_per_atom = 20.0, min_distance = 1.0, attempts = 64))]
pub(crate) fn generation_crystal(
    py: Python<'_>,
    composition: BTreeMap<String, u32>,
    space_group: u32,
    seed: u64,
    volume_per_atom: f64,
    min_distance: f64,
    attempts: usize,
) -> PyResult<Py<PyAny>> {
    let options = CrystalOptions {
        space_group,
        composition: composition_map(&composition),
        volume_per_atom,
        min_distance,
        max_attempts: attempts,
        seed,
    };
    let crystal = generate_crystal(&options).map_err(crate::map_error)?;
    json_value(py, &crystal)
}

/// Generates one crystal per space group.
#[pyfunction(name = "generate_crystal_range")]
#[pyo3(signature = (composition, groups, seed = 0, volume_per_atom = 20.0, min_distance = 1.0, attempts = 64))]
pub(crate) fn generation_crystal_range(
    py: Python<'_>,
    composition: BTreeMap<String, u32>,
    groups: Vec<u32>,
    seed: u64,
    volume_per_atom: f64,
    min_distance: f64,
    attempts: usize,
) -> PyResult<Py<PyAny>> {
    let options = CrystalOptions {
        space_group: 1,
        composition: composition_map(&composition),
        volume_per_atom,
        min_distance,
        max_attempts: attempts,
        seed,
    };
    let results = generate_crystal_range(&groups, &options).map_err(crate::map_error)?;
    let reports: Vec<serde_json::Value> = results
        .into_iter()
        .zip(groups)
        .map(|(result, group)| match result {
            Ok(crystal) => serde_json::json!({
                "space_group": group,
                "status": "ok",
                "crystal": crystal,
            }),
            Err(error) => serde_json::json!({
                "space_group": group,
                "status": "failed",
                "error": error.message(),
            }),
        })
        .collect();
    json_value(py, &reports)
}

#[allow(clippy::too_many_arguments)] // one Python keyword per option
/// Places a rigid molecule on Wyckoff positions.
#[pyfunction(name = "generate_molecule_crystal")]
#[pyo3(signature = (text, space_group = 1, copies = 1, seed = 0, volume_per_atom = 20.0, min_distance = 1.0, attempts = 64, format_hint = None))]
pub(crate) fn generation_molecule(
    py: Python<'_>,
    text: &str,
    space_group: u32,
    copies: u32,
    seed: u64,
    volume_per_atom: f64,
    min_distance: f64,
    attempts: usize,
    format_hint: Option<&str>,
) -> PyResult<Py<PyAny>> {
    let molecule = crate::structures::load_structure(text, format_hint)?;
    let options = MolecularCrystalOptions {
        space_group,
        molecule,
        copies,
        volume_per_atom,
        min_distance,
        tolerance: 0.25,
        max_attempts: attempts,
        seed,
    };
    let crystal = generate_molecular_crystal(&options).map_err(crate::map_error)?;
    json_value(py, &crystal)
}

/// Generates a layer crystal.
#[pyfunction(name = "generate_layer_crystal")]
#[pyo3(signature = (composition, layer_group = 1, seed = 0, area_per_atom = 10.0, min_distance = 1.0, attempts = 64))]
pub(crate) fn generation_layer(
    py: Python<'_>,
    composition: BTreeMap<String, u32>,
    layer_group: u32,
    seed: u64,
    area_per_atom: f64,
    min_distance: f64,
    attempts: usize,
) -> PyResult<Py<PyAny>> {
    let options = LayerCrystalOptions {
        layer_group,
        composition: composition_map(&composition),
        area_per_atom,
        min_distance,
        max_attempts: attempts,
        seed,
    };
    let crystal = generate_layer_crystal(&options).map_err(crate::map_error)?;
    json_value(py, &crystal)
}

/// Generates a rod crystal.
#[pyfunction(name = "generate_rod_crystal")]
#[pyo3(signature = (composition, rod_group = 1, seed = 0, repeat = 5.0, min_distance = 1.0, attempts = 64))]
pub(crate) fn generation_rod(
    py: Python<'_>,
    composition: BTreeMap<String, u32>,
    rod_group: u32,
    seed: u64,
    repeat: f64,
    min_distance: f64,
    attempts: usize,
) -> PyResult<Py<PyAny>> {
    let options = RodCrystalOptions {
        rod_group,
        composition: composition_map(&composition),
        repeat,
        min_distance,
        max_attempts: attempts,
        seed,
    };
    let crystal = generate_rod_crystal(&options).map_err(crate::map_error)?;
    json_value(py, &crystal)
}

/// Lists the 80 layer groups.
#[pyfunction(name = "layer_groups")]
pub(crate) fn generation_layer_groups(py: Python<'_>) -> PyResult<Py<PyAny>> {
    let groups = xtalforge::symmetry::layer_groups().map_err(crate::map_error)?;
    json_value(py, &groups)
}

/// Lists the prototype library.
#[pyfunction(name = "prototype_library")]
pub(crate) fn generation_prototypes(py: Python<'_>) -> PyResult<Py<PyAny>> {
    let library = prototype_library().map_err(crate::map_error)?;
    json_value(py, &library)
}

/// Builds a prototype by label.
#[pyfunction(name = "build_prototype")]
#[pyo3(signature = (label, scale = 1.0))]
pub(crate) fn generation_prototype(py: Python<'_>, label: &str, scale: f64) -> PyResult<Py<PyAny>> {
    let structure = build_prototype(
        label,
        &PrototypeOptions {
            scale,
            ..PrototypeOptions::default()
        },
    )
    .map_err(crate::map_error)?;
    json_value(py, &structure)
}

/// Identifies the prototype label of a structure.
#[pyfunction(name = "identify_prototype")]
#[pyo3(signature = (text, stol = 0.3, format_hint = None))]
pub(crate) fn generation_identify(
    py: Python<'_>,
    text: &str,
    stol: f64,
    format_hint: Option<&str>,
) -> PyResult<Py<PyAny>> {
    let structure = crate::structures::load_structure(text, format_hint)?;
    let identification = identify_prototype(
        &structure,
        &PrototypeOptions {
            matcher: xtalforge::symmetry::MatcherOptions {
                stol,
                ..Default::default()
            },
            ..PrototypeOptions::default()
        },
    )
    .map_err(crate::map_error)?;
    json_value(py, &identification)
}

/// Proposes chemically similar substitutions.
#[pyfunction(name = "substitution_candidates")]
#[pyo3(signature = (text, radius_tolerance = 0.15, limit = 32, same_charge = true, format_hint = None))]
pub(crate) fn generation_substitute(
    py: Python<'_>,
    text: &str,
    radius_tolerance: f64,
    limit: usize,
    same_charge: bool,
    format_hint: Option<&str>,
) -> PyResult<Py<PyAny>> {
    let structure = crate::structures::load_structure(text, format_hint)?;
    let report = substitution_report(
        &structure,
        &SubstitutionOptions {
            radius_tolerance,
            same_charge,
            limit,
        },
    )
    .map_err(crate::map_error)?;
    json_value(py, &report)
}

#[allow(clippy::too_many_arguments)] // one Python keyword per option
/// Runs an AIRSS-style random search.
#[pyfunction(name = "random_search")]
#[pyo3(signature = (composition, space_groups, structures = 8, seed = 0, relax = true, cell = false, steps = 200))]
pub(crate) fn generation_search(
    py: Python<'_>,
    composition: BTreeMap<String, u32>,
    space_groups: Vec<u32>,
    structures: usize,
    seed: u64,
    relax: bool,
    cell: bool,
    steps: usize,
) -> PyResult<Py<PyAny>> {
    let options = RandomSearchOptions {
        composition: composition_map(&composition),
        space_groups: space_group_list(space_groups),
        n_structures: structures,
        seed,
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
    let report = random_search(&options).map_err(crate::map_error)?;
    json_value(py, &report)
}

#[allow(clippy::too_many_arguments)] // one Python keyword per option
/// Runs an evolutionary search.
#[pyfunction(name = "evolutionary_search")]
#[pyo3(signature = (composition, space_groups, population = 8, generations = 5, seed = 0, cell = false, steps = 150))]
pub(crate) fn generation_evolve(
    py: Python<'_>,
    composition: BTreeMap<String, u32>,
    space_groups: Vec<u32>,
    population: usize,
    generations: usize,
    seed: u64,
    cell: bool,
    steps: usize,
) -> PyResult<Py<PyAny>> {
    let options = EvolutionaryOptions {
        population,
        generations,
        relax_cell: cell,
        max_steps: steps,
        ..EvolutionaryOptions::new(
            composition_map(&composition),
            space_group_list(space_groups),
            seed,
        )
    };
    let report = evolutionary_search(&options).map_err(crate::map_error)?;
    json_value(py, &report)
}

/// Emits a checkpointed search workflow document.
#[pyfunction(name = "search_workflow")]
#[pyo3(signature = (kind, seed_path, options, checkpoint = None, output = "search-results.json", name = "search"))]
pub(crate) fn generation_workflow(
    py: Python<'_>,
    kind: &str,
    seed_path: &str,
    options: &Bound<'_, PyAny>,
    checkpoint: Option<&str>,
    output: &str,
    name: &str,
) -> PyResult<Py<PyAny>> {
    let kind = match kind {
        "random" => SearchKind::Random,
        "evolutionary" => SearchKind::Evolutionary,
        other => {
            return Err(crate::invalid_input(format!(
                "search kind `{other}` is not random or evolutionary"
            )));
        }
    };
    let module = py.import("json")?;
    let text: String = module.call_method1("dumps", (options,))?.extract()?;
    let options: serde_json::Value = serde_json::from_str(&text)
        .map_err(|error| crate::invalid_input(format!("search options: {error}")))?;
    let document = search_workflow_document(name, kind, seed_path, options, checkpoint, output);
    json_value(py, &document)
}

#[allow(clippy::too_many_arguments)] // one Python keyword per option
/// Enumerates hypothetical MOF combinations.
#[pyfunction(name = "mof_dataset")]
#[pyo3(signature = (topologies = Vec::new(), nodes = Vec::new(), linkers = Vec::new(), limit = 32, relax = false, descriptors = true, seed = 0))]
pub(crate) fn generation_mof_dataset(
    py: Python<'_>,
    topologies: Vec<String>,
    nodes: Vec<String>,
    linkers: Vec<String>,
    limit: usize,
    relax: bool,
    descriptors: bool,
    seed: u64,
) -> PyResult<Py<PyAny>> {
    let report = generate_mof_dataset(&MofDatasetOptions {
        topologies,
        nodes,
        linkers,
        limit,
        relax,
        max_steps: 100,
        descriptors,
        pore_spacing: 1.0,
        seed,
    })
    .map_err(crate::map_error)?;
    json_value(py, &report)
}

/// Functionalises a framework over groups and fractions.
#[pyfunction(name = "functionalize_sweep")]
#[pyo3(signature = (text, groups = vec!["amino".to_owned()], fractions = vec![1.0], seed = 0, descriptors = true, format_hint = None))]
pub(crate) fn generation_functionalize(
    py: Python<'_>,
    text: &str,
    groups: Vec<String>,
    fractions: Vec<f64>,
    seed: u64,
    descriptors: bool,
    format_hint: Option<&str>,
) -> PyResult<Py<PyAny>> {
    let structure = crate::structures::load_structure(text, format_hint)?;
    let groups: Vec<xtalforge::mof::FunctionalGroup> = groups
        .iter()
        .map(|name| group(name))
        .collect::<PyResult<_>>()?;
    let sweep = functionalization_sweep(&structure, &groups, &fractions, seed, descriptors)
        .map_err(crate::map_error)?;
    json_value(py, &sweep)
}

/// Registers the Phase 20 functions.
pub(crate) fn register(module: &Bound<'_, PyModule>) -> PyResult<()> {
    module.add_function(wrap_pyfunction!(generation_crystal, module)?)?;
    module.add_function(wrap_pyfunction!(generation_crystal_range, module)?)?;
    module.add_function(wrap_pyfunction!(generation_molecule, module)?)?;
    module.add_function(wrap_pyfunction!(generation_layer, module)?)?;
    module.add_function(wrap_pyfunction!(generation_rod, module)?)?;
    module.add_function(wrap_pyfunction!(generation_layer_groups, module)?)?;
    module.add_function(wrap_pyfunction!(generation_prototypes, module)?)?;
    module.add_function(wrap_pyfunction!(generation_prototype, module)?)?;
    module.add_function(wrap_pyfunction!(generation_identify, module)?)?;
    module.add_function(wrap_pyfunction!(generation_substitute, module)?)?;
    module.add_function(wrap_pyfunction!(generation_search, module)?)?;
    module.add_function(wrap_pyfunction!(generation_evolve, module)?)?;
    module.add_function(wrap_pyfunction!(generation_workflow, module)?)?;
    module.add_function(wrap_pyfunction!(generation_mof_dataset, module)?)?;
    module.add_function(wrap_pyfunction!(generation_functionalize, module)?)?;
    Ok(())
}
