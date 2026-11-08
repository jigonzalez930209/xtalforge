//! String-based exports for the Phase 1 operations (Task 1.17.3).
//!
//! Every function takes structure text, so the browser never needs a
//! filesystem: CIF, POSCAR and XYZ are detected for the functions without an
//! explicit format hint, and `read_structure` accepts one. Results are plain
//! JavaScript objects and arrays; errors are `Error` values whose `code`
//! property carries the stable Rust error code.

use wasm_bindgen::prelude::*;

use crate::common::{
    array_i32, array_usize, array3, fill_summary, parse_any, set, summary_object, to_js_error,
};

/// Reads a structure from text and returns its summary object.
///
/// `format_hint` accepts `cif`, `poscar`, `vasp` or `xyz`, with or without a
/// leading dot. When the hint is empty the format is detected from the text.
///
/// The object has the fields of `structure_info`: `name`, `formula`, `hall`,
/// `n_sites`, `volume`, `density`, `total_charge` and `is_p1`.
///
/// # Errors
///
/// Returns a JavaScript `Error` with a stable `code` property for an unknown
/// format or invalid structure text.
#[wasm_bindgen]
pub fn read_structure(text: &str, format_hint: &str) -> Result<JsValue, JsValue> {
    let structure = if format_hint.trim().is_empty() {
        parse_any(text)
    } else {
        xtalforge::io::read_structure(text, format_hint)
    }
    .map_err(to_js_error)?;
    summary_object(&structure)
}

/// Reads a CIF and returns its summary object.
///
/// This is the forced-CIF variant of [`read_structure`].
///
/// # Errors
///
/// Returns a JavaScript `Error` with a stable `code` property for invalid
/// CIF text.
#[wasm_bindgen]
pub fn read_cif(text: &str) -> Result<JsValue, JsValue> {
    read_structure(text, "cif")
}

/// Builds a supercell and returns its summary object.
///
/// The object carries the summary fields plus `repeats`, the `[na, nb, nc]`
/// three-number array. The supercell is always P1.
///
/// # Errors
///
/// Returns a JavaScript `Error` for invalid text, zero repeats or a cell over
/// the hard size limit.
#[wasm_bindgen]
pub fn supercell_info(text: &str, na: u32, nb: u32, nc: u32) -> Result<JsValue, JsValue> {
    let structure = parse_any(text).map_err(to_js_error)?;
    let supercell = structure.supercell(na, nb, nc).map_err(to_js_error)?;
    let object = js_sys::Object::new();
    fill_summary(&object, &supercell)?;
    set(
        &object,
        "repeats",
        array3([f64::from(na), f64::from(nb), f64::from(nc)]).into(),
    );
    Ok(object.into())
}

/// Cuts a slab and returns its summary object.
///
/// The object carries the summary fields plus `hkl`, `layers`, `thickness`,
/// `vacuum` and `cut`. The slab is always P1.
///
/// # Errors
///
/// Returns a JavaScript `Error` for invalid text, a zero layer count, a
/// negative vacuum or a plane that is degenerate in the cell.
#[wasm_bindgen]
pub fn slab_info(
    text: &str,
    h: i32,
    k: i32,
    l: i32,
    layers: u32,
    vacuum: f64,
) -> Result<JsValue, JsValue> {
    let structure = parse_any(text).map_err(to_js_error)?;
    let slab = xtalforge::core::make_slab(&structure, [h, k, l], layers, vacuum, 0.0)
        .map_err(to_js_error)?;
    let object = js_sys::Object::new();
    fill_summary(&object, &slab.structure)?;
    set(&object, "hkl", array_i32(slab.hkl).into());
    set(&object, "layers", JsValue::from_f64(f64::from(slab.layers)));
    set(&object, "thickness", JsValue::from_f64(slab.thickness));
    set(&object, "vacuum", JsValue::from_f64(slab.vacuum));
    set(&object, "cut", JsValue::from_f64(slab.cut as f64));
    Ok(object.into())
}

/// Perceives the bonds of a structure.
///
/// Each element has `i`, `j`, `image` (a three-integer array), `distance`,
/// `order` and `explicit`.
///
/// # Errors
///
/// Returns a JavaScript `Error` for invalid text or a failed bond perception.
#[wasm_bindgen]
pub fn bonds(text: &str) -> Result<js_sys::Array, JsValue> {
    let structure = parse_any(text).map_err(to_js_error)?;
    let graph = structure.bond_graph().map_err(to_js_error)?;
    let array = js_sys::Array::new();
    for bond in &graph.bonds {
        let object = js_sys::Object::new();
        set(&object, "i", JsValue::from_f64(bond.i as f64));
        set(&object, "j", JsValue::from_f64(bond.j as f64));
        set(&object, "image", array_i32(bond.image).into());
        set(&object, "distance", JsValue::from_f64(bond.distance));
        set(&object, "order", JsValue::from_f64(bond.order));
        set(&object, "explicit", JsValue::from_bool(bond.explicit));
        array.push(&object.into());
    }
    Ok(array)
}

/// Finds chordless rings up to `max_size` atoms.
///
/// Each element has `size` and `atoms`, the sorted P1 atom indices.
///
/// # Errors
///
/// Returns a JavaScript `Error` for invalid text or a size limit below three.
#[wasm_bindgen]
pub fn rings(text: &str, max_size: u32) -> Result<js_sys::Array, JsValue> {
    let structure = parse_any(text).map_err(to_js_error)?;
    let found = xtalforge::core::rings(&structure, max_size as usize).map_err(to_js_error)?;
    let array = js_sys::Array::new();
    for ring in found {
        let object = js_sys::Object::new();
        set(&object, "size", JsValue::from_f64(ring.size() as f64));
        set(&object, "atoms", array_usize(&ring.atoms).into());
        array.push(&object.into());
    }
    Ok(array)
}

/// Finds the functional groups of a structure.
///
/// Each element has `kind` (the stable group name) and `atoms`, the sorted P1
/// atom indices.
///
/// # Errors
///
/// Returns a JavaScript `Error` for invalid text or failed bond perception.
#[wasm_bindgen]
pub fn groups(text: &str) -> Result<js_sys::Array, JsValue> {
    let structure = parse_any(text).map_err(to_js_error)?;
    let found = xtalforge::core::find_groups(&structure).map_err(to_js_error)?;
    let array = js_sys::Array::new();
    for group in found {
        let object = js_sys::Object::new();
        set(&object, "kind", JsValue::from_str(group.kind.name()));
        set(&object, "atoms", array_usize(&group.atoms).into());
        array.push(&object.into());
    }
    Ok(array)
}

/// Finds the periodic neighbors within `cutoff` Angstrom.
///
/// Each element has `i`, `j`, `image` (the lattice translation of `j`) and
/// `distance`.
///
/// # Errors
///
/// Returns a JavaScript `Error` for invalid text or a non-positive cutoff.
#[wasm_bindgen]
pub fn neighbors(text: &str, cutoff: f64) -> Result<js_sys::Array, JsValue> {
    let structure = parse_any(text).map_err(to_js_error)?;
    let found = xtalforge::core::neighbors(&structure, cutoff).map_err(to_js_error)?;
    let array = js_sys::Array::new();
    for neighbor in found {
        let object = js_sys::Object::new();
        set(&object, "i", JsValue::from_f64(neighbor.i as f64));
        set(&object, "j", JsValue::from_f64(neighbor.j as f64));
        set(&object, "image", array_i32(neighbor.image).into());
        set(&object, "distance", JsValue::from_f64(neighbor.distance));
        array.push(&object.into());
    }
    Ok(array)
}

/// Converts structure text to another format and returns the file contents.
///
/// `format_hint` is the input format (`cif`, `poscar` or `xyz`, or empty to
/// detect it) and `format` is the output format (`cif`, `poscar` or `xyz`).
///
/// # Errors
///
/// Returns a JavaScript `Error` for invalid text or an unknown format.
#[wasm_bindgen]
pub fn convert_structure(text: &str, format_hint: &str, format: &str) -> Result<String, JsValue> {
    let structure = if format_hint.trim().is_empty() {
        parse_any(text)
    } else {
        xtalforge::io::read_structure(text, format_hint)
    }
    .map_err(to_js_error)?;
    xtalforge::io::write_structure(&structure, format).map_err(to_js_error)
}

/// Builds a supercell and returns it as CIF text.
///
/// # Errors
///
/// Returns a JavaScript `Error` for invalid text, zero repeats, a refused
/// size or a writer failure.
#[wasm_bindgen]
pub fn supercell_text(text: &str, na: u32, nb: u32, nc: u32) -> Result<String, JsValue> {
    let structure = parse_any(text).map_err(to_js_error)?;
    let built = structure.supercell(na, nb, nc).map_err(to_js_error)?;
    xtalforge::io::write_structure(&built, "cif").map_err(to_js_error)
}

/// Cuts a slab and returns it as CIF text.
///
/// # Errors
///
/// Returns a JavaScript `Error` for invalid text, a bad plane or a writer
/// failure.
#[wasm_bindgen]
pub fn slab_text(
    text: &str,
    h: i32,
    k: i32,
    l: i32,
    layers: u32,
    vacuum: f64,
) -> Result<String, JsValue> {
    let structure = parse_any(text).map_err(to_js_error)?;
    let slab = xtalforge::core::make_slab(&structure, [h, k, l], layers, vacuum, 0.0)
        .map_err(to_js_error)?;
    xtalforge::io::write_structure(&slab.structure, "cif").map_err(to_js_error)
}

/// Minimum-image distance between two P1 atoms, in Angstrom.
///
/// # Errors
///
/// Returns a JavaScript `Error` for invalid text or an out-of-range index.
#[wasm_bindgen]
pub fn distance(text: &str, i: u32, j: u32) -> Result<f64, JsValue> {
    let structure = parse_any(text).map_err(to_js_error)?;
    xtalforge::core::distance(&structure, i as usize, j as usize).map_err(to_js_error)
}

/// Angle `i-j-k` in degrees, using minimum images.
///
/// # Errors
///
/// Returns a JavaScript `Error` for invalid text, an out-of-range index or a
/// degenerate angle.
#[wasm_bindgen]
pub fn angle(text: &str, i: u32, j: u32, k: u32) -> Result<f64, JsValue> {
    let structure = parse_any(text).map_err(to_js_error)?;
    xtalforge::core::angle(&structure, i as usize, j as usize, k as usize).map_err(to_js_error)
}

/// Torsion `i-j-k-l` in degrees, using minimum images.
///
/// # Errors
///
/// Returns a JavaScript `Error` for invalid text, an out-of-range index or a
/// collinear segment.
#[wasm_bindgen]
pub fn torsion(text: &str, i: u32, j: u32, k: u32, l: u32) -> Result<f64, JsValue> {
    let structure = parse_any(text).map_err(to_js_error)?;
    xtalforge::core::torsion(&structure, i as usize, j as usize, k as usize, l as usize)
        .map_err(to_js_error)
}
