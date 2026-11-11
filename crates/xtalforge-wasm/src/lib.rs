//! WebAssembly bindings for `xtalforge`.
//!
//! The module exports the same hello-world surface as the Rust, CLI and
//! Python packages. The wasm-bindgen attributes are compiled only for
//! `wasm32-unknown-unknown`; a host build of this crate is intentionally
//! empty so `cargo test --workspace` stays green.

#[cfg(target_arch = "wasm32")]
use wasm_bindgen::prelude::*;

#[cfg(target_arch = "wasm32")]
mod agent_export;
#[cfg(target_arch = "wasm32")]
mod analysis_export;
#[cfg(target_arch = "wasm32")]
mod common;
#[cfg(target_arch = "wasm32")]
mod forcefield;
#[cfg(target_arch = "wasm32")]
mod mof_export;
#[cfg(target_arch = "wasm32")]
mod phase1;
#[cfg(target_arch = "wasm32")]
mod phase3;
#[cfg(target_arch = "wasm32")]
mod phase6_export;
#[cfg(target_arch = "wasm32")]
mod powder_export;
#[cfg(target_arch = "wasm32")]
mod scene_export;
#[cfg(target_arch = "wasm32")]
mod structure_text;
#[cfg(target_arch = "wasm32")]
mod symmetry;
#[cfg(target_arch = "wasm32")]
mod symmetry_analysis;

#[cfg(target_arch = "wasm32")]
use common::{array3, array6, set, to_js_error};

#[cfg(target_arch = "wasm32")]
pub use agent_export::{AgentSession, agent_capabilities, agent_run, agent_verbs};
#[cfg(target_arch = "wasm32")]
pub use analysis_export::{hydrogens_info, net_invariants_info, porosity_info, zeo_parse_info};
#[cfg(target_arch = "wasm32")]
pub use forcefield::{
    ff_charges_info, ff_energy, ff_ir_info, ff_optimize_info, ff_optimize_text, ff_prepare_info,
    ff_thermo_info, ff_types_info, ff_vibrations_info,
};
#[cfg(target_arch = "wasm32")]
pub use mof_export::{build_mof_info, mof_blocks_info, mof_topologies_info};
#[cfg(target_arch = "wasm32")]
pub use phase1::{
    bonds, build_c60_info, build_graphene_info, build_nanotube_info, convert_structure, groups,
    mesh_stl, neighbors, pxrd_info, read_cif, read_structure, rings, scene_summary, slab_info,
    slab_text, supercell_info, supercell_text,
};
#[cfg(target_arch = "wasm32")]
pub use phase3::{
    WasmDocument, param_schema, project_info, project_pack, workspace_archive,
    workspace_read_archive,
};
#[cfg(target_arch = "wasm32")]
pub use phase6_export::{
    build_c20_info, build_carbon_recipe_info, polymer_info, ribbon_info, schwarzite_info,
    twist_info, ztc_info,
};
#[cfg(target_arch = "wasm32")]
pub use powder_export::{
    powder_automatic_info, powder_cells_info, powder_index_info, powder_pawley_info,
    powder_refine_info, rietx_inventory_info,
};
#[cfg(target_arch = "wasm32")]
pub use scene_export::scene_data;
#[cfg(target_arch = "wasm32")]
pub use structure_text::{
    asymmetrize_text, build_c20_text, build_c60_text, build_carbon_text, build_graphene_text,
    build_mof_text, build_nanotube_text, build_polymer_text, build_ribbon_text,
    build_schwarzite_text, build_twist_text, build_ztc_text, change_setting_text, hydrogens_text,
    lower_to_subgroup_text, primitive_text, reduce_text, set_space_group_text, standardize_text,
    to_p1_text,
};
#[cfg(target_arch = "wasm32")]
pub use symmetry::{
    angle, asymmetrize_info, detect_symmetry, distance, lower_to_subgroup_info,
    set_space_group_info, subgroups_info, to_p1_info, torsion,
};
#[cfg(target_arch = "wasm32")]
pub use symmetry_analysis::{
    brillouin_zone, change_setting_info, diagnose_info, kpath_info, layer_group_info,
    magnetic_group_info, pointgroup_info, primitive_info, reduce_info, standardize_info,
    tolerance_ladder_info, wyckoff_info,
};

/// Installs the panic hook in debug builds only.
#[cfg(all(target_arch = "wasm32", debug_assertions))]
#[wasm_bindgen(start)]
pub fn install_panic_hook() {
    console_error_panic_hook::set_once();
}

/// Version of the workspace.
#[cfg(target_arch = "wasm32")]
#[wasm_bindgen]
pub fn version() -> String {
    xtalforge::core::version().to_owned()
}

/// Volume of a cell in cubic Angstrom.
///
/// # Errors
///
/// Returns a JavaScript `Error` whose `code` property is the stable error
/// code from the Rust core.
#[cfg(target_arch = "wasm32")]
#[wasm_bindgen]
pub fn cell_volume(
    a: f64,
    b: f64,
    c: f64,
    alpha: f64,
    beta: f64,
    gamma: f64,
) -> Result<f64, JsValue> {
    xtalforge::core::cell_volume(a, b, c, alpha, beta, gamma).map_err(to_js_error)
}

/// Unit-cell information from the six parameters.
///
/// # Errors
///
/// Returns a JavaScript `Error` with a stable `code` property for an invalid
/// cell.
#[cfg(target_arch = "wasm32")]
#[wasm_bindgen]
pub fn lattice_info(
    a: f64,
    b: f64,
    c: f64,
    alpha: f64,
    beta: f64,
    gamma: f64,
) -> Result<JsValue, JsValue> {
    let lattice =
        xtalforge::core::lattice_info(a, b, c, alpha, beta, gamma).map_err(to_js_error)?;
    let object = js_sys::Object::new();
    set(&object, "parameters", array6(lattice.parameters()).into());
    set(&object, "volume", JsValue::from_f64(lattice.volume()));
    let matrix = js_sys::Array::new();
    for row in lattice.matrix() {
        matrix.push(&array3(row).into());
    }
    set(&object, "matrix", matrix.into());
    set(
        &object,
        "right_handed",
        JsValue::from_bool(lattice.is_right_handed()),
    );
    Ok(object.into())
}

/// Minimum-image fractional vector and distance.
///
/// # Errors
///
/// Returns a JavaScript `Error` with a stable `code` property for an invalid
/// cell or a delta that is not three numbers.
#[cfg(target_arch = "wasm32")]
#[wasm_bindgen]
pub fn lattice_minimum_image(
    a: f64,
    b: f64,
    c: f64,
    alpha: f64,
    beta: f64,
    gamma: f64,
    delta: Vec<f64>,
) -> Result<JsValue, JsValue> {
    let [dx, dy, dz] = delta.as_slice() else {
        return Err(crate::common::invalid_input(
            "delta must contain three fractional coordinates",
        ));
    };
    let (vector, distance) =
        xtalforge::core::lattice_minimum_image(a, b, c, alpha, beta, gamma, [*dx, *dy, *dz])
            .map_err(to_js_error)?;
    let object = js_sys::Object::new();
    set(&object, "delta", array3(vector).into());
    set(&object, "distance", JsValue::from_f64(distance));
    Ok(object.into())
}

/// Soft and hard bounds for a profile and budget.
///
/// # Errors
///
/// Returns a JavaScript `Error` for an unknown budget.
#[cfg(target_arch = "wasm32")]
#[wasm_bindgen]
pub fn limits_bounds(profile: &str, budget: &str) -> Result<js_sys::Array, JsValue> {
    let (soft, hard) = xtalforge::core::limits_bounds(profile, budget).map_err(to_js_error)?;
    let array = js_sys::Array::new();
    array.push(&JsValue::from_f64(soft as f64));
    array.push(&hard.map_or(JsValue::NULL, |value| JsValue::from_f64(value as f64)));
    Ok(array)
}

/// Generic limit check for a profile, budget and estimate.
///
/// # Errors
///
/// Returns a JavaScript `Error` for an unknown budget.
#[cfg(target_arch = "wasm32")]
#[wasm_bindgen]
pub fn limits_check(profile: &str, budget: &str, estimate: f64) -> Result<JsValue, JsValue> {
    let estimate = if estimate.is_finite() && estimate > 0.0 {
        estimate.round() as u64
    } else {
        0
    };
    let verdict = xtalforge::core::limits_check(profile, budget, estimate).map_err(to_js_error)?;
    let object = js_sys::Object::new();
    set(&object, "level", JsValue::from_str(verdict.level.name()));
    set(
        &object,
        "estimate",
        JsValue::from_f64(verdict.estimate as f64),
    );
    set(&object, "sentence", JsValue::from_str(&verdict.sentence));
    set(&object, "refused", JsValue::from_bool(verdict.refused()));
    set(&object, "warned", JsValue::from_bool(verdict.warned()));
    Ok(object.into())
}

/// Summary of a structure JSON document.
///
/// # Errors
///
/// Returns a JavaScript `Error` with a stable `code` property for invalid
/// input.
#[cfg(target_arch = "wasm32")]
#[wasm_bindgen]
pub fn structure_info(text: &str) -> Result<JsValue, JsValue> {
    let summary = xtalforge::core::structure_info(text).map_err(to_js_error)?;
    let structure = json_structure(text).map_err(to_js_error)?;
    let object = js_sys::Object::new();
    set(&object, "name", JsValue::from_str(&summary.name));
    set(&object, "formula", JsValue::from_str(&summary.formula));
    set(&object, "hall", JsValue::from_str(&summary.hall));
    set(
        &object,
        "n_sites",
        JsValue::from_f64(summary.n_sites as f64),
    );
    set(
        &object,
        "n_atoms",
        JsValue::from_f64(structure.p1().map_err(to_js_error)?.n_atoms() as f64),
    );
    set(&object, "volume", JsValue::from_f64(summary.volume));
    set(&object, "density", JsValue::from_f64(summary.density));
    set(
        &object,
        "total_charge",
        JsValue::from_f64(summary.total_charge),
    );
    set(&object, "is_p1", JsValue::from_bool(summary.is_p1));
    Ok(object.into())
}

/// Parses the bare or versioned structure JSON that `structure_info` accepts.
#[cfg(target_arch = "wasm32")]
fn json_structure(text: &str) -> Result<xtalforge::core::Structure, xtalforge::core::Error> {
    if text.contains("\"schema_version\"") {
        Ok(xtalforge::core::StructureDocument::from_json(text)?.structure)
    } else {
        serde_json::from_str(text).map_err(|error| {
            xtalforge::core::Error::InvalidInput(format!("invalid structure JSON: {error}"))
        })
    }
}

/// Reference data for an element by symbol or atomic number.
///
/// # Errors
///
/// Returns a JavaScript `Error` with a stable `code` property when the input
/// is not a known element.
#[cfg(target_arch = "wasm32")]
#[wasm_bindgen]
pub fn element_info(input: &str) -> Result<JsValue, JsValue> {
    let element = xtalforge::core::element_info(input).map_err(to_js_error)?;
    let object = js_sys::Object::new();
    set(&object, "z", JsValue::from_f64(f64::from(element.z)));
    set(&object, "symbol", JsValue::from_str(&element.symbol));
    set(&object, "name", JsValue::from_str(&element.name));
    set(&object, "mass", JsValue::from_f64(element.mass));
    set(
        &object,
        "covalent_radius",
        JsValue::from_f64(element.covalent_radius),
    );
    set(&object, "vdw_radius", JsValue::from_f64(element.vdw_radius));
    let color = js_sys::Array::new();
    for channel in element.color {
        color.push(&JsValue::from_f64(f64::from(channel)));
    }
    set(&object, "color", color.into());
    set(&object, "is_metal", JsValue::from_bool(element.is_metal));
    set(
        &object,
        "valence",
        element
            .valence
            .map_or(JsValue::NULL, |value| JsValue::from_f64(f64::from(value))),
    );
    set(
        &object,
        "electronegativity",
        element
            .electronegativity
            .map_or(JsValue::NULL, JsValue::from_f64),
    );
    set(
        &object,
        "pyykko_single",
        element
            .pyykko_single
            .map_or(JsValue::NULL, JsValue::from_f64),
    );
    set(
        &object,
        "pyykko_double",
        element
            .pyykko_double
            .map_or(JsValue::NULL, JsValue::from_f64),
    );
    set(
        &object,
        "pyykko_aromatic",
        element
            .pyykko_aromatic
            .map_or(JsValue::NULL, JsValue::from_f64),
    );
    Ok(object.into())
}
