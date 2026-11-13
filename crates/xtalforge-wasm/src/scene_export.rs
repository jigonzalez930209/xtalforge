//! Full scene geometry for browser renderers (Task 3.14.4).
//!
//! `scene_summary` reports only counts because the first demos needed no more.
//! A 3D viewer needs the geometry itself: atom positions, radii and colours,
//! bond endpoints and the cell edges. This module exports exactly the
//! [`xtalforge::view::SceneData`] the Rust, CLI and Python surfaces build, so a
//! browser renders the same frame without a second geometry implementation.

use wasm_bindgen::prelude::*;

use xtalforge::core::Error;
use xtalforge::view::{DrawStyle, ViewSettings, build_mesh, build_scene, write_gltf, write_obj};

use crate::common::{array3, parse_any, set, to_js_error};

fn draw_style(name: &str) -> Result<DrawStyle, JsValue> {
    match name {
        "ball-stick" | "ball_stick" => Ok(DrawStyle::BallStick),
        "space-fill" | "space_fill" => Ok(DrawStyle::SpaceFill),
        "wireframe" => Ok(DrawStyle::Wireframe),
        "sticks" => Ok(DrawStyle::Sticks),
        other => Err(to_js_error(Error::InvalidInput(format!(
            "unknown draw style {other:?}; use ball-stick, space-fill, wireframe or sticks"
        )))),
    }
}

/// Builds the full scene geometry for a structure.
///
/// The object has three arrays:
///
/// - `atoms`: `element`, `position` (Cartesian Angstrom), `radius`, `color`
///   (0-255 RGB), `index` in the P1 cell and `hidden`.
/// - `bonds`: `i`, `j`, `start`, `end`, `radius`, `color` and `order`.
/// - `cell_edges`: pairs of Cartesian endpoints.
///
/// `settings` repeats the sanitized style, atom scale, bond radius, cell and
/// axes flags and background colour used to build the frame.
///
/// # Errors
///
/// Returns a JavaScript `Error` with a stable `code` property for invalid
/// structure text or an unknown draw style.
#[wasm_bindgen]
pub fn scene_data(text: &str, style: &str, scale: f64) -> Result<JsValue, JsValue> {
    let structure = parse_any(text).map_err(to_js_error)?;
    let settings = ViewSettings {
        style: draw_style(style)?,
        atom_scale: scale,
        ..ViewSettings::default()
    };
    let scene = build_scene(&structure, &settings).map_err(to_js_error)?;

    let atoms = js_sys::Array::new();
    for atom in &scene.atoms {
        let object = js_sys::Object::new();
        set(&object, "element", JsValue::from_str(&atom.element));
        set(&object, "position", array3(atom.position).into());
        set(&object, "radius", JsValue::from_f64(atom.radius));
        let color = js_sys::Array::new();
        for channel in atom.color {
            color.push(&JsValue::from_f64(f64::from(channel)));
        }
        set(&object, "color", color.into());
        set(&object, "index", JsValue::from_f64(atom.index as f64));
        set(&object, "hidden", JsValue::from_bool(atom.hidden));
        atoms.push(&object.into());
    }

    let bonds = js_sys::Array::new();
    for bond in &scene.bonds {
        let object = js_sys::Object::new();
        set(&object, "i", JsValue::from_f64(bond.i as f64));
        set(&object, "j", JsValue::from_f64(bond.j as f64));
        set(&object, "start", array3(bond.start).into());
        set(&object, "end", array3(bond.end).into());
        set(&object, "radius", JsValue::from_f64(bond.radius));
        let color = js_sys::Array::new();
        for channel in bond.color {
            color.push(&JsValue::from_f64(f64::from(channel)));
        }
        set(&object, "color", color.into());
        set(&object, "order", JsValue::from_f64(bond.order));
        bonds.push(&object.into());
    }

    let cell_edges = js_sys::Array::new();
    for (start, end) in &scene.cell_edges {
        let pair = js_sys::Array::new();
        pair.push(&array3(*start).into());
        pair.push(&array3(*end).into());
        cell_edges.push(&pair.into());
    }

    let used = scene.settings;
    let style_name = match used.style {
        DrawStyle::BallStick => "ball-stick",
        DrawStyle::SpaceFill => "space-fill",
        DrawStyle::Wireframe => "wireframe",
        DrawStyle::Sticks => "sticks",
    };
    let settings = js_sys::Object::new();
    set(&settings, "style", JsValue::from_str(style_name));
    set(&settings, "atom_scale", JsValue::from_f64(used.atom_scale));
    set(
        &settings,
        "bond_radius",
        JsValue::from_f64(used.bond_radius),
    );
    set(&settings, "show_cell", JsValue::from_bool(used.show_cell));
    set(&settings, "show_axes", JsValue::from_bool(used.show_axes));
    let background = js_sys::Array::new();
    for channel in used.background {
        background.push(&JsValue::from_f64(f64::from(channel)));
    }
    set(&settings, "background", background.into());

    let object = js_sys::Object::new();
    set(&object, "atoms", atoms.into());
    set(&object, "bonds", bonds.into());
    set(&object, "cell_edges", cell_edges.into());
    set(&object, "settings", settings.into());
    Ok(object.into())
}

fn mesh_text(text: &str, style: &str, scale: f64, gltf: bool) -> Result<String, JsValue> {
    let structure = parse_any(text).map_err(to_js_error)?;
    let settings = ViewSettings {
        style: draw_style(style)?,
        atom_scale: scale,
        ..ViewSettings::default()
    };
    let scene = build_scene(&structure, &settings).map_err(to_js_error)?;
    let mesh = build_mesh(&scene);
    if gltf {
        write_gltf(&mesh, structure.name()).map_err(to_js_error)
    } else {
        Ok(write_obj(&mesh, structure.name()))
    }
}

/// Builds a Wavefront OBJ mesh.
#[wasm_bindgen]
pub fn mesh_obj(text: &str, style: &str, scale: f64) -> Result<String, JsValue> {
    mesh_text(text, style, scale, false)
}

/// Builds a glTF 2.0 JSON mesh.
#[wasm_bindgen]
pub fn mesh_gltf(text: &str, style: &str, scale: f64) -> Result<String, JsValue> {
    mesh_text(text, style, scale, true)
}
