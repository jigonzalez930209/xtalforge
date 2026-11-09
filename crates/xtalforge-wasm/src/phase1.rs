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
        let path = js_sys::Array::new();
        for step in &ring.path {
            let pair = js_sys::Array::new();
            pair.push(&JsValue::from_f64(step.atom as f64));
            pair.push(&array_i32(step.image).into());
            path.push(&pair.into());
        }
        set(&object, "path", path.into());
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

/// Summarizes a scene for the given draw style.
///
/// # Errors
///
/// Returns a JavaScript `Error` for invalid text or an unknown style.
#[wasm_bindgen]
pub fn scene_summary(text: &str, style: &str, scale: f64) -> Result<JsValue, JsValue> {
    let structure = parse_any(text).map_err(to_js_error)?;
    let style = view_style(style)?;
    let settings = xtalforge::view::ViewSettings {
        style,
        atom_scale: scale,
        ..Default::default()
    };
    let scene = xtalforge::view::build_scene(&structure, &settings).map_err(to_js_error)?;
    let object = js_sys::Object::new();
    set(
        &object,
        "atoms",
        JsValue::from_f64(scene.atoms.len() as f64),
    );
    set(
        &object,
        "bonds",
        JsValue::from_f64(scene.bonds.len() as f64),
    );
    set(
        &object,
        "cell_edges",
        JsValue::from_f64(scene.cell_edges.len() as f64),
    );
    Ok(object.into())
}

/// Builds an ASCII STL string.
///
/// # Errors
///
/// Returns a JavaScript `Error` for invalid text or an unknown style.
#[wasm_bindgen]
pub fn mesh_stl(text: &str, style: &str, scale: f64) -> Result<String, JsValue> {
    let structure = parse_any(text).map_err(to_js_error)?;
    let style = view_style(style)?;
    let settings = xtalforge::view::ViewSettings {
        style,
        atom_scale: scale,
        ..Default::default()
    };
    let scene = xtalforge::view::build_scene(&structure, &settings).map_err(to_js_error)?;
    Ok(xtalforge::view::write_stl(
        &xtalforge::view::build_mesh(&scene),
        structure.name(),
    ))
}

fn view_style(name: &str) -> Result<xtalforge::view::DrawStyle, JsValue> {
    match name {
        "ball-stick" | "ball_stick" => Ok(xtalforge::view::DrawStyle::BallStick),
        "space-fill" | "space_fill" => Ok(xtalforge::view::DrawStyle::SpaceFill),
        "wireframe" => Ok(xtalforge::view::DrawStyle::Wireframe),
        "sticks" => Ok(xtalforge::view::DrawStyle::Sticks),
        other => Err(JsValue::from_str(&format!(
            "invalid_input: unknown draw style {other:?}"
        ))),
    }
}

/// Calculates a PXRD pattern.
///
/// # Errors
///
/// Returns a JavaScript `Error` for invalid text, radiation or options.
#[wasm_bindgen]
pub fn pxrd_info(
    text: &str,
    radiation: &str,
    max_two_theta: f64,
    step: f64,
) -> Result<JsValue, JsValue> {
    use xtalforge::analysis::{PxrdOptions, Radiation, calculate_pxrd};
    let structure = parse_any(text).map_err(to_js_error)?;
    let radiation = match radiation {
        "cu-ka" | "cuka" => Radiation::CuKa,
        "cu-ka1" | "cuka1" => Radiation::CuKa1,
        "mo-ka" | "moka" => Radiation::MoKa,
        "co-ka" | "coka" => Radiation::CoKa,
        "fe-ka" | "feka" => Radiation::FeKa,
        "cr-ka" | "crka" => Radiation::CrKa,
        "ag-ka" | "agka" => Radiation::AgKa,
        other => {
            return Err(JsValue::from_str(&format!(
                "invalid_input: unknown radiation {other:?}"
            )));
        }
    };
    let pattern = calculate_pxrd(
        &structure,
        &PxrdOptions {
            radiation,
            two_theta_range: [5.0, max_two_theta],
            step,
            ..PxrdOptions::default()
        },
    )
    .map_err(to_js_error)?;
    let object = js_sys::Object::new();
    let reflections = js_sys::Array::new();
    for reflection in &pattern.reflections {
        let item = js_sys::Object::new();
        set(&item, "hkl", array_i32(reflection.hkl).into());
        set(&item, "d", JsValue::from_f64(reflection.d));
        set(&item, "two_theta", JsValue::from_f64(reflection.two_theta));
        set(
            &item,
            "multiplicity",
            JsValue::from_f64(reflection.multiplicity as f64),
        );
        set(&item, "intensity", JsValue::from_f64(reflection.intensity));
        reflections.push(&item.into());
    }
    set(&object, "reflections", reflections.into());
    let points = js_sys::Array::new();
    for point in &pattern.points {
        let pair = js_sys::Array::new();
        pair.push(&JsValue::from_f64(point.two_theta));
        pair.push(&JsValue::from_f64(point.intensity));
        points.push(&pair.into());
    }
    set(&object, "points", points.into());
    Ok(object.into())
}

/// Builds a graphene sheet and returns its summary.
///
/// # Errors
///
/// Returns a JavaScript `Error` for invalid repeats or spacing.
#[wasm_bindgen]
pub fn build_graphene_info(
    nx: u32,
    ny: u32,
    layers: u32,
    interlayer: f64,
    vacuum: f64,
) -> Result<JsValue, JsValue> {
    let structure =
        xtalforge::build::graphene(nx, ny, layers, interlayer, vacuum).map_err(to_js_error)?;
    summary_object(&structure)
}

/// Builds a nanotube and returns its summary.
///
/// # Errors
///
/// Returns a JavaScript `Error` for a bad chiral vector or length.
#[wasm_bindgen]
pub fn build_nanotube_info(n: i32, m: i32, length: u32, vacuum: f64) -> Result<JsValue, JsValue> {
    let structure = xtalforge::build::nanotube(n, m, length, vacuum).map_err(to_js_error)?;
    summary_object(&structure)
}

/// Builds C60 and returns its summary.
///
/// # Errors
///
/// Returns a JavaScript `Error` when the embedded data is invalid.
#[wasm_bindgen]
pub fn build_c60_info() -> Result<JsValue, JsValue> {
    let structure = xtalforge::build::fullerene(60).map_err(to_js_error)?;
    summary_object(&structure)
}
