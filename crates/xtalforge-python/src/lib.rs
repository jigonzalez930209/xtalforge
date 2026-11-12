//! Python bindings for `xtalforge`.
//!
//! The module layout is `xtalforge._xtalforge`, re-exported by the Python
//! package in `python/xtalforge`. The public surface is kept at parity with
//! the Rust, CLI and WASM APIs.

#![allow(unsafe_code)] // pyo3 macros expand to FFI code; this crate writes none.

use pyo3::create_exception;
use pyo3::exceptions::PyException;
use pyo3::prelude::*;

mod agent;
mod analysis;
mod classes;
mod lebail;
mod ml;
mod phase10;
mod phase3;
mod phase4;
mod phase5;
mod phase6;
mod phase6_extra;
mod phase7;
mod phase9;
mod phase9_analysis;
mod phase9_local;
mod plugins;
mod structures;
mod symmetry;
mod symmetry_analysis;
mod view;

use analysis::{angle, bonds, distance, groups, neighbors, pxrd, ring_paths, rings, torsion};
use classes::{Lattice, P1Cell, Site, Structure};
use lebail::powder_lebail;
use phase3::{
    Document, param_schema, project_pack, project_read, workspace_archive, workspace_open,
    workspace_read_archive, workspace_read_run, workspace_runs, workspace_write_run,
};
use phase4::{
    ff_charges, ff_energy, ff_ir, ff_optimize, ff_prepare, ff_thermo, ff_types, ff_vibrations,
};
use phase5::{
    add_hydrogens, net_invariants, net_overlay, porosity, rcsr_names, zeo_parse, zeo_run,
};
use phase6::{build_mof, mof_blocks, mof_topologies};
use phase6_extra::{
    build_c20, build_carbon, build_polymer, build_schwarzite, build_ztc, ribbon, twist,
};
use phase7::{
    powder_automatic, powder_cells, powder_index, powder_pawley, powder_refine,
    powder_supported_groups, rietx,
};
use phase9::{
    deform, edit_structure, group_structures, interpolate, interpolate_structures,
    kabsch_alignment, match_structures, optimal_supercell, periodic_rmsd, rattle,
};
use phase9_analysis::{
    adsorption_sites, convex_hull, elastic_tensor, eos_fit, interstitials, lattice_matches,
    orderings, place_molecule, slab_terminations, special_quasirandom, stack_slabs, substitutions,
    unique_miller_indices, vacancies, wulff,
};
use phase9_local::{
    bond_valence_sums, coordination, fingerprint_similarity, oxidation_states, radial_distribution,
    structure_fingerprint, voronoi,
};
use phase10::{
    mof_charge_balance, mof_check, mof_deconstruct, mof_deduplicate, mof_defects, mof_descriptors,
    mof_disorder, mof_energy_grid, mof_forcefield, mof_functionalize, mof_gcmc, mof_guests,
    mof_interpenetrate, mof_mixed, mof_mofid, mof_open_metal_sites, mof_pores, mof_raspa,
    mof_read_database, mof_screen, mof_solvent_remove, mof_stack_layers, mof_widom,
};
use structures::{
    build_c60, build_graphene, build_nanotube, read_cif, read_structure, slab, supercell,
};
use symmetry::{
    asymmetrize, detect_symmetry, lower_to_subgroup, set_space_group, subgroups, to_p1,
};
use symmetry_analysis::{
    band_path, brillouin_zone, change_setting, detect_layer_symmetry, detect_magnetic_symmetry,
    diagnose, molecular_point_group, primitive_cell, reduce_structure, standardize,
    tolerance_ladder, wyckoff_assignments,
};
use view::{mesh_stl, scene_data, scene_summary};

create_exception!(
    _xtalforge,
    XtalforgeError,
    PyException,
    "Base class for every xtalforge error."
);
create_exception!(
    _xtalforge,
    InputError,
    XtalforgeError,
    "A geometry or value does not satisfy the API preconditions."
);

/// Version of the workspace.
#[pyfunction]
fn version() -> &'static str {
    xtalforge::core::version()
}

fn map_error(error: xtalforge::core::Error) -> PyErr {
    let args = (error.code(), error.message());
    match error {
        xtalforge::core::Error::InvalidInput(_) => InputError::new_err(args),
        _ => XtalforgeError::new_err(args),
    }
}

fn map_mof_error(error: xtalforge::mof::Error) -> PyErr {
    let args = (error.code(), error.message());
    match error {
        xtalforge::mof::Error::InvalidInput(_) => InputError::new_err(args),
        _ => XtalforgeError::new_err(args),
    }
}

fn map_agent_error(error: xtalforge::agent::Error) -> PyErr {
    let args = (error.code(), error.to_string());
    match error {
        xtalforge::agent::Error::Argument { .. } => InputError::new_err(args),
        _ => XtalforgeError::new_err(args),
    }
}

fn map_ff_error(error: xtalforge::ff::Error) -> PyErr {
    let args = (error.code(), error.message());
    match error {
        xtalforge::ff::Error::InvalidInput(_) => InputError::new_err(args),
        _ => XtalforgeError::new_err(args),
    }
}

pub(crate) fn invalid_input(message: impl Into<String>) -> PyErr {
    InputError::new_err(("invalid_input", message.into()))
}

/// Volume of a cell in cubic Angstrom.
#[pyfunction]
#[pyo3(signature = (a, b, c, alpha, beta, gamma))]
fn cell_volume(a: f64, b: f64, c: f64, alpha: f64, beta: f64, gamma: f64) -> PyResult<f64> {
    xtalforge::core::cell_volume(a, b, c, alpha, beta, gamma).map_err(map_error)
}

/// Unit-cell information from the six parameters.
#[pyfunction]
#[pyo3(signature = (a, b, c, alpha, beta, gamma))]
fn lattice_info<'py>(
    py: Python<'py>,
    a: f64,
    b: f64,
    c: f64,
    alpha: f64,
    beta: f64,
    gamma: f64,
) -> PyResult<Bound<'py, pyo3::types::PyDict>> {
    let lattice = xtalforge::core::lattice_info(a, b, c, alpha, beta, gamma).map_err(map_error)?;
    let dict = pyo3::types::PyDict::new(py);
    dict.set_item("parameters", lattice.parameters().to_vec())?;
    dict.set_item("volume", lattice.volume())?;
    dict.set_item(
        "matrix",
        lattice
            .matrix()
            .iter()
            .map(|row| row.to_vec())
            .collect::<Vec<_>>(),
    )?;
    dict.set_item("right_handed", lattice.is_right_handed())?;
    Ok(dict)
}

/// Minimum-image fractional vector and distance.
#[pyfunction]
#[pyo3(signature = (a, b, c, alpha, beta, gamma, delta))]
#[allow(clippy::too_many_arguments)] // one Python argument per cell parameter
fn lattice_minimum_image<'py>(
    py: Python<'py>,
    a: f64,
    b: f64,
    c: f64,
    alpha: f64,
    beta: f64,
    gamma: f64,
    delta: Vec<f64>,
) -> PyResult<Bound<'py, pyo3::types::PyDict>> {
    let [dx, dy, dz] = delta.as_slice() else {
        return Err(InputError::new_err((
            "invalid_input",
            "delta must contain three fractional coordinates",
        )));
    };
    let delta = [*dx, *dy, *dz];
    let (vector, distance) =
        xtalforge::core::lattice_minimum_image(a, b, c, alpha, beta, gamma, delta)
            .map_err(map_error)?;
    let dict = pyo3::types::PyDict::new(py);
    dict.set_item("delta", vector.to_vec())?;
    dict.set_item("distance", distance)?;
    Ok(dict)
}

/// Soft and hard bounds for a profile and budget.
#[pyfunction]
fn limits_bounds(profile: &str, budget: &str) -> PyResult<(u64, Option<u64>)> {
    xtalforge::core::limits_bounds(profile, budget).map_err(map_error)
}

/// Generic limit check for a profile, budget and estimate.
#[pyfunction]
fn limits_check<'py>(
    py: Python<'py>,
    profile: &str,
    budget: &str,
    estimate: u64,
) -> PyResult<Bound<'py, pyo3::types::PyDict>> {
    let verdict = xtalforge::core::limits_check(profile, budget, estimate).map_err(map_error)?;
    let refused = verdict.refused();
    let warned = verdict.warned();
    let dict = pyo3::types::PyDict::new(py);
    dict.set_item("level", verdict.level.name())?;
    dict.set_item("estimate", verdict.estimate)?;
    dict.set_item("sentence", &verdict.sentence)?;
    dict.set_item("refused", refused)?;
    dict.set_item("warned", warned)?;
    Ok(dict)
}

/// Summary of a structure JSON document.
#[pyfunction]
fn structure_info<'py>(py: Python<'py>, text: &str) -> PyResult<Bound<'py, pyo3::types::PyDict>> {
    let summary = xtalforge::core::structure_info(text).map_err(map_error)?;
    let dict = pyo3::types::PyDict::new(py);
    dict.set_item("name", summary.name)?;
    dict.set_item("formula", summary.formula)?;
    dict.set_item("hall", summary.hall)?;
    dict.set_item("n_sites", summary.n_sites)?;
    dict.set_item("volume", summary.volume)?;
    dict.set_item("density", summary.density)?;
    dict.set_item("total_charge", summary.total_charge)?;
    dict.set_item("is_p1", summary.is_p1)?;
    Ok(dict)
}

/// Reference data for an element by symbol or atomic number.
#[pyfunction]
fn element_info<'py>(
    py: Python<'py>,
    input: &Bound<'_, PyAny>,
) -> PyResult<Bound<'py, pyo3::types::PyDict>> {
    let text = if let Ok(z) = input.extract::<u8>() {
        z.to_string()
    } else {
        input.extract::<String>()?
    };
    let element = xtalforge::core::element_info(&text).map_err(map_error)?;
    let dict = pyo3::types::PyDict::new(py);
    dict.set_item("z", element.z)?;
    dict.set_item("symbol", element.symbol)?;
    dict.set_item("name", element.name)?;
    dict.set_item("mass", element.mass)?;
    dict.set_item("covalent_radius", element.covalent_radius)?;
    dict.set_item("vdw_radius", element.vdw_radius)?;
    dict.set_item("color", element.color.to_vec())?;
    dict.set_item("is_metal", element.is_metal)?;
    dict.set_item("valence", element.valence)?;
    dict.set_item("electronegativity", element.electronegativity)?;
    dict.set_item("pyykko_single", element.pyykko_single)?;
    dict.set_item("pyykko_double", element.pyykko_double)?;
    dict.set_item("pyykko_aromatic", element.pyykko_aromatic)?;
    Ok(dict)
}

#[pymodule]
fn _xtalforge(module: &Bound<'_, PyModule>) -> PyResult<()> {
    module.add("__version__", xtalforge::core::version())?;
    module.add("XtalforgeError", module.py().get_type::<XtalforgeError>())?;
    module.add("InputError", module.py().get_type::<InputError>())?;
    module.add_function(wrap_pyfunction!(version, module)?)?;
    module.add_function(wrap_pyfunction!(cell_volume, module)?)?;
    module.add_function(wrap_pyfunction!(element_info, module)?)?;
    module.add_function(wrap_pyfunction!(lattice_info, module)?)?;
    module.add_function(wrap_pyfunction!(lattice_minimum_image, module)?)?;
    module.add_function(wrap_pyfunction!(limits_bounds, module)?)?;
    module.add_function(wrap_pyfunction!(limits_check, module)?)?;
    module.add_function(wrap_pyfunction!(structure_info, module)?)?;
    module.add_function(wrap_pyfunction!(build_graphene, module)?)?;
    module.add_function(wrap_pyfunction!(build_nanotube, module)?)?;
    module.add_function(wrap_pyfunction!(build_c60, module)?)?;
    module.add_function(wrap_pyfunction!(angle, module)?)?;
    module.add_function(wrap_pyfunction!(bonds, module)?)?;
    module.add_function(wrap_pyfunction!(distance, module)?)?;
    module.add_function(wrap_pyfunction!(groups, module)?)?;
    module.add_function(wrap_pyfunction!(neighbors, module)?)?;
    module.add_function(wrap_pyfunction!(read_cif, module)?)?;
    module.add_function(wrap_pyfunction!(read_structure, module)?)?;
    module.add_function(wrap_pyfunction!(rings, module)?)?;
    module.add_function(wrap_pyfunction!(ring_paths, module)?)?;
    module.add_function(wrap_pyfunction!(slab, module)?)?;
    module.add_function(wrap_pyfunction!(supercell, module)?)?;
    module.add_function(wrap_pyfunction!(torsion, module)?)?;
    module.add_function(wrap_pyfunction!(pxrd, module)?)?;
    module.add_function(wrap_pyfunction!(detect_symmetry, module)?)?;
    module.add_function(wrap_pyfunction!(asymmetrize, module)?)?;
    module.add_function(wrap_pyfunction!(to_p1, module)?)?;
    module.add_function(wrap_pyfunction!(set_space_group, module)?)?;
    module.add_function(wrap_pyfunction!(subgroups, module)?)?;
    module.add_function(wrap_pyfunction!(lower_to_subgroup, module)?)?;
    module.add_function(wrap_pyfunction!(standardize, module)?)?;
    module.add_function(wrap_pyfunction!(primitive_cell, module)?)?;
    module.add_function(wrap_pyfunction!(reduce_structure, module)?)?;
    module.add_function(wrap_pyfunction!(wyckoff_assignments, module)?)?;
    module.add_function(wrap_pyfunction!(brillouin_zone, module)?)?;
    module.add_function(wrap_pyfunction!(tolerance_ladder, module)?)?;
    module.add_function(wrap_pyfunction!(diagnose, module)?)?;
    module.add_function(wrap_pyfunction!(detect_layer_symmetry, module)?)?;
    module.add_function(wrap_pyfunction!(detect_magnetic_symmetry, module)?)?;
    module.add_function(wrap_pyfunction!(change_setting, module)?)?;
    module.add_function(wrap_pyfunction!(band_path, module)?)?;
    module.add_function(wrap_pyfunction!(scene_summary, module)?)?;
    module.add_function(wrap_pyfunction!(scene_data, module)?)?;
    module.add_function(wrap_pyfunction!(mesh_stl, module)?)?;
    module.add_function(wrap_pyfunction!(molecular_point_group, module)?)?;
    module.add_class::<Lattice>()?;
    module.add_class::<Site>()?;
    module.add_class::<Structure>()?;
    module.add_class::<P1Cell>()?;
    module.add_class::<Document>()?;
    module.add_function(wrap_pyfunction!(param_schema, module)?)?;
    module.add_function(wrap_pyfunction!(project_pack, module)?)?;
    module.add_function(wrap_pyfunction!(project_read, module)?)?;
    module.add_function(wrap_pyfunction!(workspace_open, module)?)?;
    module.add_function(wrap_pyfunction!(workspace_write_run, module)?)?;
    module.add_function(wrap_pyfunction!(workspace_runs, module)?)?;
    module.add_function(wrap_pyfunction!(workspace_read_run, module)?)?;
    module.add_function(wrap_pyfunction!(workspace_archive, module)?)?;
    module.add_function(wrap_pyfunction!(workspace_read_archive, module)?)?;
    module.add_function(wrap_pyfunction!(ff_energy, module)?)?;
    module.add_function(wrap_pyfunction!(ff_optimize, module)?)?;
    module.add_function(wrap_pyfunction!(ff_charges, module)?)?;
    module.add_function(wrap_pyfunction!(ff_vibrations, module)?)?;
    module.add_function(wrap_pyfunction!(ff_ir, module)?)?;
    module.add_function(wrap_pyfunction!(ff_thermo, module)?)?;
    module.add_function(wrap_pyfunction!(ff_types, module)?)?;
    module.add_function(wrap_pyfunction!(ff_prepare, module)?)?;
    module.add_function(wrap_pyfunction!(porosity, module)?)?;
    module.add_function(wrap_pyfunction!(zeo_parse, module)?)?;
    module.add_function(wrap_pyfunction!(zeo_run, module)?)?;
    module.add_function(wrap_pyfunction!(net_invariants, module)?)?;
    module.add_function(wrap_pyfunction!(net_overlay, module)?)?;
    module.add_function(wrap_pyfunction!(rcsr_names, module)?)?;
    module.add_function(wrap_pyfunction!(add_hydrogens, module)?)?;
    module.add_function(wrap_pyfunction!(mof_topologies, module)?)?;
    module.add_function(wrap_pyfunction!(mof_blocks, module)?)?;
    module.add_function(wrap_pyfunction!(build_mof, module)?)?;
    module.add_function(wrap_pyfunction!(build_c20, module)?)?;
    module.add_function(wrap_pyfunction!(twist, module)?)?;
    module.add_function(wrap_pyfunction!(ribbon, module)?)?;
    module.add_function(wrap_pyfunction!(build_ztc, module)?)?;
    module.add_function(wrap_pyfunction!(build_carbon, module)?)?;
    module.add_function(wrap_pyfunction!(build_schwarzite, module)?)?;
    module.add_function(wrap_pyfunction!(build_polymer, module)?)?;
    module.add_function(wrap_pyfunction!(powder_automatic, module)?)?;
    module.add_function(wrap_pyfunction!(powder_index, module)?)?;
    module.add_function(wrap_pyfunction!(powder_pawley, module)?)?;
    module.add_function(wrap_pyfunction!(powder_lebail, module)?)?;
    module.add_function(wrap_pyfunction!(powder_refine, module)?)?;
    module.add_function(wrap_pyfunction!(powder_cells, module)?)?;
    module.add_function(wrap_pyfunction!(powder_supported_groups, module)?)?;
    module.add_function(wrap_pyfunction!(rietx, module)?)?;
    module.add_function(wrap_pyfunction!(match_structures, module)?)?;
    module.add_function(wrap_pyfunction!(periodic_rmsd, module)?)?;
    module.add_function(wrap_pyfunction!(kabsch_alignment, module)?)?;
    module.add_function(wrap_pyfunction!(group_structures, module)?)?;
    module.add_function(wrap_pyfunction!(deform, module)?)?;
    module.add_function(wrap_pyfunction!(rattle, module)?)?;
    module.add_function(wrap_pyfunction!(edit_structure, module)?)?;
    module.add_function(wrap_pyfunction!(interpolate_structures, module)?)?;
    module.add_function(wrap_pyfunction!(interpolate, module)?)?;
    module.add_function(wrap_pyfunction!(optimal_supercell, module)?)?;
    module.add_function(wrap_pyfunction!(voronoi, module)?)?;
    module.add_function(wrap_pyfunction!(coordination, module)?)?;
    module.add_function(wrap_pyfunction!(bond_valence_sums, module)?)?;
    module.add_function(wrap_pyfunction!(oxidation_states, module)?)?;
    module.add_function(wrap_pyfunction!(radial_distribution, module)?)?;
    module.add_function(wrap_pyfunction!(structure_fingerprint, module)?)?;
    module.add_function(wrap_pyfunction!(fingerprint_similarity, module)?)?;
    module.add_function(wrap_pyfunction!(unique_miller_indices, module)?)?;
    module.add_function(wrap_pyfunction!(slab_terminations, module)?)?;
    module.add_function(wrap_pyfunction!(wulff, module)?)?;
    module.add_function(wrap_pyfunction!(adsorption_sites, module)?)?;
    module.add_function(wrap_pyfunction!(place_molecule, module)?)?;
    module.add_function(wrap_pyfunction!(vacancies, module)?)?;
    module.add_function(wrap_pyfunction!(substitutions, module)?)?;
    module.add_function(wrap_pyfunction!(interstitials, module)?)?;
    module.add_function(wrap_pyfunction!(orderings, module)?)?;
    module.add_function(wrap_pyfunction!(special_quasirandom, module)?)?;
    module.add_function(wrap_pyfunction!(lattice_matches, module)?)?;
    module.add_function(wrap_pyfunction!(stack_slabs, module)?)?;
    module.add_function(wrap_pyfunction!(elastic_tensor, module)?)?;
    module.add_function(wrap_pyfunction!(eos_fit, module)?)?;
    module.add_function(wrap_pyfunction!(convex_hull, module)?)?;
    module.add_function(wrap_pyfunction!(mof_check, module)?)?;
    module.add_function(wrap_pyfunction!(mof_solvent_remove, module)?)?;
    module.add_function(wrap_pyfunction!(mof_disorder, module)?)?;
    module.add_function(wrap_pyfunction!(mof_charge_balance, module)?)?;
    module.add_function(wrap_pyfunction!(mof_deconstruct, module)?)?;
    module.add_function(wrap_pyfunction!(mof_mofid, module)?)?;
    module.add_function(wrap_pyfunction!(mof_open_metal_sites, module)?)?;
    module.add_function(wrap_pyfunction!(mof_descriptors, module)?)?;
    module.add_function(wrap_pyfunction!(mof_pores, module)?)?;
    module.add_function(wrap_pyfunction!(mof_guests, module)?)?;
    module.add_function(wrap_pyfunction!(mof_forcefield, module)?)?;
    module.add_function(wrap_pyfunction!(mof_widom, module)?)?;
    module.add_function(wrap_pyfunction!(mof_gcmc, module)?)?;
    module.add_function(wrap_pyfunction!(mof_energy_grid, module)?)?;
    module.add_function(wrap_pyfunction!(mof_raspa, module)?)?;
    module.add_function(wrap_pyfunction!(mof_functionalize, module)?)?;
    module.add_function(wrap_pyfunction!(mof_defects, module)?)?;
    module.add_function(wrap_pyfunction!(mof_interpenetrate, module)?)?;
    module.add_function(wrap_pyfunction!(mof_stack_layers, module)?)?;
    module.add_function(wrap_pyfunction!(mof_mixed, module)?)?;
    module.add_function(wrap_pyfunction!(mof_read_database, module)?)?;
    module.add_function(wrap_pyfunction!(mof_screen, module)?)?;
    module.add_function(wrap_pyfunction!(mof_deduplicate, module)?)?;
    agent::register(module)?;
    plugins::register(module)?;
    ml::register(module)?;
    Ok(())
}
