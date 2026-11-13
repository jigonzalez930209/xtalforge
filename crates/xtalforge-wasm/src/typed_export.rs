//! Typed array exports for zero-copy rendering and plotting (Task 12.3.3).

use wasm_bindgen::prelude::*;

use crate::common::{parse_any, set, to_js_error};

/// Scene geometry as typed arrays.
///
/// The object carries `positions` (Float64Array, xyz), `radii`
/// (Float32Array), `colors` (Uint8Array, rgb), `bond_positions`
/// (Float64Array, two xyz per bond), `bond_colors` (Uint8Array) and
/// `cell_edges` (Float64Array, two xyz per edge).
///
/// # Errors
///
/// Returns a JavaScript `Error` for invalid text.
#[wasm_bindgen]
pub fn scene_buffers(text: &str, format_hint: &str) -> Result<JsValue, JsValue> {
    let structure = if format_hint.trim().is_empty() {
        parse_any(text)
    } else {
        xtalforge::io::read_structure(text, format_hint)
    }
    .map_err(to_js_error)?;
    let scene = xtalforge::view::build_scene(&structure, &xtalforge::view::ViewSettings::default())
        .map_err(to_js_error)?;
    let mut positions = Vec::with_capacity(scene.atoms.len() * 3);
    let mut radii = Vec::with_capacity(scene.atoms.len());
    let mut colors = Vec::with_capacity(scene.atoms.len() * 3);
    for atom in &scene.atoms {
        positions.extend_from_slice(&atom.position);
        radii.push(atom.radius as f32);
        colors.extend_from_slice(&atom.color);
    }
    let mut bond_positions = Vec::with_capacity(scene.bonds.len() * 6);
    let mut bond_colors = Vec::with_capacity(scene.bonds.len() * 3);
    for bond in &scene.bonds {
        bond_positions.extend_from_slice(&bond.start);
        bond_positions.extend_from_slice(&bond.end);
        bond_colors.extend_from_slice(&bond.color);
    }
    let mut cell_edges = Vec::with_capacity(scene.cell_edges.len() * 6);
    for (start, end) in &scene.cell_edges {
        cell_edges.extend_from_slice(start);
        cell_edges.extend_from_slice(end);
    }
    let object = js_sys::Object::new();
    set(
        &object,
        "positions",
        js_sys::Float64Array::from(positions.as_slice()).into(),
    );
    set(
        &object,
        "radii",
        js_sys::Float32Array::from(radii.as_slice()).into(),
    );
    set(
        &object,
        "colors",
        js_sys::Uint8Array::from(colors.as_slice()).into(),
    );
    set(
        &object,
        "bond_positions",
        js_sys::Float64Array::from(bond_positions.as_slice()).into(),
    );
    set(
        &object,
        "bond_colors",
        js_sys::Uint8Array::from(bond_colors.as_slice()).into(),
    );
    set(
        &object,
        "cell_edges",
        js_sys::Float64Array::from(cell_edges.as_slice()).into(),
    );
    Ok(object.into())
}

/// PXRD `two_theta` and `intensity` as Float64Arrays.
///
/// # Errors
///
/// Returns a JavaScript `Error` for invalid text.
#[wasm_bindgen]
pub fn pxrd_buffers(text: &str, format_hint: &str) -> Result<JsValue, JsValue> {
    let structure = if format_hint.trim().is_empty() {
        parse_any(text)
    } else {
        xtalforge::io::read_structure(text, format_hint)
    }
    .map_err(to_js_error)?;
    let pattern = xtalforge::analysis::calculate_pxrd(
        &structure,
        &xtalforge::analysis::PxrdOptions::default(),
    )
    .map_err(to_js_error)?;
    let mut two_theta = Vec::with_capacity(pattern.points.len());
    let mut intensity = Vec::with_capacity(pattern.points.len());
    for point in &pattern.points {
        two_theta.push(point.two_theta);
        intensity.push(point.intensity);
    }
    let object = js_sys::Object::new();
    set(
        &object,
        "two_theta",
        js_sys::Float64Array::from(two_theta.as_slice()).into(),
    );
    set(
        &object,
        "intensity",
        js_sys::Float64Array::from(intensity.as_slice()).into(),
    );
    Ok(object.into())
}

/// Mesh triangles as vertex, normal and index buffers.
///
/// # Errors
///
/// Returns a JavaScript `Error` for invalid text.
#[wasm_bindgen]
pub fn mesh_buffers(text: &str, format_hint: &str) -> Result<JsValue, JsValue> {
    let structure = if format_hint.trim().is_empty() {
        parse_any(text)
    } else {
        xtalforge::io::read_structure(text, format_hint)
    }
    .map_err(to_js_error)?;
    let scene = xtalforge::view::build_scene(&structure, &xtalforge::view::ViewSettings::default())
        .map_err(to_js_error)?;
    let triangles = xtalforge::view::build_mesh(&scene);
    let mut vertices = Vec::with_capacity(triangles.len() * 9);
    let mut normals = Vec::with_capacity(triangles.len() * 9);
    let mut indices = Vec::with_capacity(triangles.len() * 3);
    for (index, triangle) in triangles.iter().enumerate() {
        for vertex in triangle.vertices {
            vertices.extend_from_slice(&vertex);
            normals.extend_from_slice(&triangle.normal);
        }
        let base = (index * 3) as u32;
        indices.extend_from_slice(&[base, base + 1, base + 2]);
    }
    let object = js_sys::Object::new();
    set(
        &object,
        "vertices",
        js_sys::Float64Array::from(vertices.as_slice()).into(),
    );
    set(
        &object,
        "normals",
        js_sys::Float64Array::from(normals.as_slice()).into(),
    );
    set(
        &object,
        "indices",
        js_sys::Uint32Array::from(indices.as_slice()).into(),
    );
    Ok(object.into())
}
