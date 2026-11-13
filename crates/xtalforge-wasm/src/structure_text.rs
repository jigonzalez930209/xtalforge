//! Text-returning variants of the analysis and builder exports.
//!
//! The summary objects are enough for a form, but a 3D playground must feed
//! the result of a transform or a builder back into `read_structure`. These
//! wrappers write the same [`xtalforge::core::Structure`] the summary exports
//! build as CIF, so the browser has no second implementation of any operation.

use wasm_bindgen::prelude::*;

use xtalforge::core::Structure;

use crate::common::{parse_any, to_js_error};

fn cif(structure: &Structure) -> Result<String, JsValue> {
    xtalforge::io::write_cif(structure).map_err(to_js_error)
}

#[cfg(feature = "full")]
fn mof_error(error: xtalforge::mof::Error) -> JsValue {
    let value = js_sys::Object::new();
    let _ = js_sys::Reflect::set(
        &value,
        &JsValue::from_str("code"),
        &JsValue::from_str(error.code()),
    );
    let _ = js_sys::Reflect::set(
        &value,
        &JsValue::from_str("message"),
        &JsValue::from_str(&error.message()),
    );
    value.into()
}

/// Detects the symmetry and returns the asymmetric unit as CIF.
///
/// # Errors
///
/// Returns a JavaScript `Error` for invalid text or a failed search.
#[cfg(feature = "symmetry")]
#[wasm_bindgen]
pub fn asymmetrize_text(text: &str, symprec: f64) -> Result<String, JsValue> {
    let structure = parse_any(text).map_err(to_js_error)?;
    let reduced = xtalforge::symmetry::asymmetrize(
        &structure,
        xtalforge::symmetry::SymmetryOptions {
            symprec,
            ..Default::default()
        },
    )
    .map_err(to_js_error)?;
    cif(&reduced)
}

/// Expands a structure to P1 and returns it as CIF.
///
/// # Errors
///
/// Returns a JavaScript `Error` for invalid text or an empty structure.
#[cfg(feature = "symmetry")]
#[wasm_bindgen]
pub fn to_p1_text(text: &str) -> Result<String, JsValue> {
    let structure = parse_any(text).map_err(to_js_error)?;
    let expanded = xtalforge::symmetry::to_p1(&structure).map_err(to_js_error)?;
    cif(&expanded)
}

/// Applies a space group and returns the result as CIF.
///
/// # Errors
///
/// Returns a JavaScript `Error` for invalid text, an unknown symbol or mode.
#[cfg(feature = "symmetry")]
#[wasm_bindgen]
pub fn set_space_group_text(text: &str, symbol: &str, mode: &str) -> Result<String, JsValue> {
    let structure = parse_any(text).map_err(to_js_error)?;
    let group = xtalforge::symmetry::SpaceGroup::parse(symbol).map_err(to_js_error)?;
    let mode = match mode {
        "reinterpret" => xtalforge::symmetry::SetMode::Reinterpret,
        "impose" => xtalforge::symmetry::SetMode::Impose,
        other => {
            return Err(to_js_error(xtalforge::core::Error::InvalidInput(format!(
                "unknown set mode {other:?}; use reinterpret or impose"
            ))));
        }
    };
    let report =
        xtalforge::symmetry::set_space_group(&structure, &group, mode).map_err(to_js_error)?;
    cif(&report.structure)
}

/// Standardizes a structure and returns it as CIF.
///
/// # Errors
///
/// Returns a JavaScript `Error` for invalid text or a failed search.
#[cfg(feature = "symmetry")]
#[wasm_bindgen]
pub fn standardize_text(text: &str, symprec: f64) -> Result<String, JsValue> {
    let structure = parse_any(text).map_err(to_js_error)?;
    let result = xtalforge::symmetry::standardize(
        &structure,
        xtalforge::symmetry::SymmetryOptions {
            symprec,
            ..Default::default()
        },
    )
    .map_err(to_js_error)?;
    cif(&result.structure)
}

/// Builds the primitive cell and returns it as CIF.
///
/// # Errors
///
/// Returns a JavaScript `Error` for invalid text or a failed search.
#[cfg(feature = "symmetry")]
#[wasm_bindgen]
pub fn primitive_text(text: &str, symprec: f64) -> Result<String, JsValue> {
    let structure = parse_any(text).map_err(to_js_error)?;
    let result = xtalforge::symmetry::primitive_cell(
        &structure,
        xtalforge::symmetry::SymmetryOptions {
            symprec,
            ..Default::default()
        },
    )
    .map_err(to_js_error)?;
    cif(&result.structure)
}

/// Reduces a structure and returns it as CIF.
///
/// # Errors
///
/// Returns a JavaScript `Error` for invalid text or an unknown algorithm.
#[cfg(feature = "symmetry")]
#[wasm_bindgen]
pub fn reduce_text(text: &str, algorithm: &str) -> Result<String, JsValue> {
    let structure = parse_any(text).map_err(to_js_error)?;
    let algorithm = match algorithm {
        "niggli" => xtalforge::symmetry::Reduction::Niggli,
        "delaunay" => xtalforge::symmetry::Reduction::Delaunay,
        "minkowski" => xtalforge::symmetry::Reduction::Minkowski,
        other => {
            return Err(to_js_error(xtalforge::core::Error::InvalidInput(format!(
                "unknown reduction algorithm {other:?}; use niggli, delaunay or minkowski"
            ))));
        }
    };
    let reduced =
        xtalforge::symmetry::reduce_structure(&structure, algorithm).map_err(to_js_error)?;
    cif(&reduced.structure)
}

/// Changes the space-group setting and returns the result as CIF.
///
/// # Errors
///
/// Returns a JavaScript `Error` for invalid text or a Hall number out of range.
#[cfg(feature = "symmetry")]
#[wasm_bindgen]
pub fn change_setting_text(text: &str, hall_number: u32) -> Result<String, JsValue> {
    let structure = parse_any(text).map_err(to_js_error)?;
    let changed =
        xtalforge::symmetry::change_setting(&structure, hall_number).map_err(to_js_error)?;
    cif(&changed.structure)
}

/// Lowers a structure to a subgroup and returns it as CIF.
///
/// # Errors
///
/// Returns a JavaScript `Error` for invalid text or an unknown subgroup.
#[cfg(feature = "symmetry")]
#[wasm_bindgen]
pub fn lower_to_subgroup_text(
    text: &str,
    subgroup_number: u32,
    tolerance: f64,
) -> Result<String, JsValue> {
    let structure = parse_any(text).map_err(to_js_error)?;
    let lowered = xtalforge::symmetry::lower_to_subgroup(&structure, subgroup_number, tolerance)
        .map_err(to_js_error)?;
    cif(&lowered)
}

/// Builds a graphene sheet and returns it as CIF.
///
/// # Errors
///
/// Returns a JavaScript `Error` for invalid repeats or spacing.
#[cfg(feature = "build")]
#[wasm_bindgen]
pub fn build_graphene_text(
    nx: u32,
    ny: u32,
    layers: u32,
    interlayer: f64,
    vacuum: f64,
) -> Result<String, JsValue> {
    let structure =
        xtalforge::build::graphene(nx, ny, layers, interlayer, vacuum).map_err(to_js_error)?;
    cif(&structure)
}

/// Builds a nanotube and returns it as CIF.
///
/// # Errors
///
/// Returns a JavaScript `Error` for a bad chiral vector or length.
#[cfg(feature = "build")]
#[wasm_bindgen]
pub fn build_nanotube_text(n: i32, m: i32, length: u32, vacuum: f64) -> Result<String, JsValue> {
    let structure = xtalforge::build::nanotube(n, m, length, vacuum).map_err(to_js_error)?;
    cif(&structure)
}

/// Builds C60 and returns it as CIF.
///
/// # Errors
///
/// Returns a JavaScript `Error` when the embedded data is invalid.
#[cfg(feature = "build")]
#[wasm_bindgen]
pub fn build_c60_text() -> Result<String, JsValue> {
    let structure = xtalforge::build::fullerene(60).map_err(to_js_error)?;
    cif(&structure)
}

/// Builds C20 and returns it as CIF.
///
/// # Errors
///
/// Returns a JavaScript `Error` when the embedded data is invalid.
#[cfg(feature = "build")]
#[wasm_bindgen]
pub fn build_c20_text() -> Result<String, JsValue> {
    let structure = xtalforge::build::face_spiral(20).map_err(to_js_error)?;
    cif(&structure)
}

/// Builds a MOF and returns it as CIF.
///
/// # Errors
///
/// Returns a JavaScript `Error` for an unknown topology or block.
#[cfg(feature = "full")]
#[wasm_bindgen]
pub fn build_mof_text(topology: &str, node: &str, linker: &str) -> Result<String, JsValue> {
    let request = xtalforge::mof::MofRequest::new(topology, node, linker);
    let built = xtalforge::mof::build_mof(&request).map_err(mof_error)?;
    cif(&built.structure)
}

/// Builds a honeycomb nanoribbon and returns it as CIF.
///
/// # Errors
///
/// Returns a JavaScript `Error` for invalid dimensions.
#[cfg(feature = "build")]
#[wasm_bindgen]
pub fn build_ribbon_text(width: u32, length: u32, vacuum: f64) -> Result<String, JsValue> {
    let structure = xtalforge::build::nanoribbon(width, length, vacuum).map_err(to_js_error)?;
    cif(&structure)
}

/// Builds two twisted graphene layers and returns them as CIF.
///
/// # Errors
///
/// Returns a JavaScript `Error` for invalid repeats or angle.
#[cfg(feature = "build")]
#[wasm_bindgen]
pub fn build_twist_text(n: u32, m: u32, theta: f64, layers: u32) -> Result<String, JsValue> {
    let structure =
        xtalforge::build::commensurate_twist(n, m, theta, layers).map_err(to_js_error)?;
    cif(&structure)
}

/// Builds a gyroid ZTC carbon model and returns it as CIF.
///
/// # Errors
///
/// Returns a JavaScript `Error` for invalid options.
#[cfg(feature = "build")]
#[wasm_bindgen]
pub fn build_ztc_text(samples: usize, cell: f64, isovalue: f64) -> Result<String, JsValue> {
    let structure = xtalforge::build::build_ztc(&xtalforge::build::ZtcOptions {
        samples,
        cell,
        isovalue,
    })
    .map_err(to_js_error)?;
    cif(&structure)
}

/// Builds a Schwarz P carbon model and returns it as CIF.
///
/// # Errors
///
/// Returns a JavaScript `Error` for invalid options.
#[cfg(feature = "build")]
#[wasm_bindgen]
pub fn build_schwarzite_text(samples: usize, cell: f64, isovalue: f64) -> Result<String, JsValue> {
    let structure = xtalforge::build::build_schwarzite(&xtalforge::build::ZtcOptions {
        samples,
        cell,
        isovalue,
    })
    .map_err(to_js_error)?;
    cif(&structure)
}

/// Builds a polymer chain and returns it as CIF.
///
/// # Errors
///
/// Returns a JavaScript `Error` for an invalid SMILES or density.
#[cfg(feature = "build")]
#[wasm_bindgen]
pub fn build_polymer_text(smiles: &str, monomers: u32, density: f64) -> Result<String, JsValue> {
    let build = xtalforge::build::build_polymer(
        smiles,
        monomers,
        &xtalforge::build::PolymerOptions { density },
    )
    .map_err(to_js_error)?;
    cif(&build.structure)
}

/// Builds a disordered carbon from a net recipe and returns it as CIF.
///
/// The other recipe fields keep their defaults; `build_carbon_recipe_info`
/// exposes the full recipe and its statistics.
///
/// # Errors
///
/// Returns a JavaScript `Error` for an unknown net or invalid options.
#[cfg(feature = "build")]
#[wasm_bindgen]
pub fn build_carbon_text(
    net: &str,
    repeat: u32,
    density: f64,
    seed: u32,
) -> Result<String, JsValue> {
    let recipe = xtalforge::build::CarbonRecipe {
        net: net.to_owned(),
        repeat: [repeat as usize; 3],
        density,
        seed: u64::from(seed),
        ..xtalforge::build::CarbonRecipe::default()
    };
    let build = xtalforge::build::build_carbon(&recipe).map_err(to_js_error)?;
    cif(&build.structure)
}

/// Adds hydrogens at the missing valences and returns the result as CIF.
///
/// # Errors
///
/// Returns a JavaScript `Error` for invalid text or a bad distance.
#[cfg(feature = "build")]
#[wasm_bindgen]
pub fn hydrogens_text(text: &str, distance: f64) -> Result<String, JsValue> {
    let structure = parse_any(text).map_err(to_js_error)?;
    let options = xtalforge::build::HydrogenOptions {
        distance,
        ..Default::default()
    };
    let filled = xtalforge::build::add_hydrogens(&structure, &options).map_err(to_js_error)?;
    cif(&filled)
}
