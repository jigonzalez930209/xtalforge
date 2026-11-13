//! View and STL bindings (Task 3.16.2).

use pyo3::prelude::*;
use pyo3::types::{PyDict, PyList};

use xtalforge::view::{
    DrawStyle, ViewSettings, build_mesh, build_scene, write_gltf, write_obj, write_stl,
};

use crate::structures::load_structure;

fn style(name: &str) -> PyResult<DrawStyle> {
    match name {
        "ball-stick" | "ball_stick" => Ok(DrawStyle::BallStick),
        "space-fill" | "space_fill" => Ok(DrawStyle::SpaceFill),
        "wireframe" => Ok(DrawStyle::Wireframe),
        "sticks" => Ok(DrawStyle::Sticks),
        other => Err(crate::invalid_input(format!(
            "unknown draw style {other:?}; use ball-stick, space-fill, wireframe or sticks"
        ))),
    }
}

/// Builds the scene summary for a structure.
#[pyfunction]
#[pyo3(signature = (text, style = "ball-stick".to_string(), scale = 1.0, format_hint = None))]
pub(crate) fn scene_summary<'py>(
    py: Python<'py>,
    text: &str,
    style: String,
    scale: f64,
    format_hint: Option<&str>,
) -> PyResult<Bound<'py, PyDict>> {
    let structure = load_structure(text, format_hint)?;
    let settings = ViewSettings {
        style: self::style(&style)?,
        atom_scale: scale,
        ..ViewSettings::default()
    };
    let scene = build_scene(&structure, &settings).map_err(crate::map_error)?;
    let dict = PyDict::new(py);
    dict.set_item("atoms", scene.atoms.len())?;
    dict.set_item("bonds", scene.bonds.len())?;
    dict.set_item("cell_edges", scene.cell_edges.len())?;
    dict.set_item(
        "elements",
        scene
            .atoms
            .iter()
            .map(|atom| atom.element.clone())
            .collect::<Vec<_>>(),
    )?;
    dict.set_item(
        "radii",
        scene
            .atoms
            .iter()
            .map(|atom| atom.radius)
            .collect::<Vec<_>>(),
    )?;
    Ok(dict)
}

/// Builds the full scene geometry for a structure.
///
/// Returns the atoms (element, position, radius, colour and index), the bonds
/// (indices, endpoints, radius, colour and order) and the cell edge segments,
/// the same data the Rust `build_scene` and the WASM `scene_data` produce.
#[pyfunction]
#[pyo3(signature = (text, style = "ball-stick".to_string(), scale = 1.0, format_hint = None))]
pub(crate) fn scene_data<'py>(
    py: Python<'py>,
    text: &str,
    style: String,
    scale: f64,
    format_hint: Option<&str>,
) -> PyResult<Bound<'py, PyDict>> {
    let structure = load_structure(text, format_hint)?;
    let settings = ViewSettings {
        style: self::style(&style)?,
        atom_scale: scale,
        ..ViewSettings::default()
    };
    let scene = build_scene(&structure, &settings).map_err(crate::map_error)?;

    let atoms = PyList::empty(py);
    for atom in &scene.atoms {
        let item = PyDict::new(py);
        item.set_item("element", &atom.element)?;
        item.set_item("position", atom.position)?;
        item.set_item("radius", atom.radius)?;
        item.set_item("color", atom.color)?;
        item.set_item("index", atom.index)?;
        item.set_item("hidden", atom.hidden)?;
        atoms.append(item)?;
    }

    let bonds = PyList::empty(py);
    for bond in &scene.bonds {
        let item = PyDict::new(py);
        item.set_item("i", bond.i)?;
        item.set_item("j", bond.j)?;
        item.set_item("start", bond.start)?;
        item.set_item("end", bond.end)?;
        item.set_item("radius", bond.radius)?;
        item.set_item("color", bond.color)?;
        item.set_item("order", bond.order)?;
        bonds.append(item)?;
    }

    let cell_edges = PyList::empty(py);
    for (start, end) in &scene.cell_edges {
        cell_edges.append((start, end))?;
    }

    let settings = PyDict::new(py);
    let style_name = match scene.settings.style {
        DrawStyle::BallStick => "ball-stick",
        DrawStyle::SpaceFill => "space-fill",
        DrawStyle::Wireframe => "wireframe",
        DrawStyle::Sticks => "sticks",
    };
    settings.set_item("style", style_name)?;
    settings.set_item("atom_scale", scene.settings.atom_scale)?;
    settings.set_item("bond_radius", scene.settings.bond_radius)?;
    settings.set_item("show_cell", scene.settings.show_cell)?;
    settings.set_item("show_axes", scene.settings.show_axes)?;
    settings.set_item("background", scene.settings.background)?;

    let dict = PyDict::new(py);
    dict.set_item("atoms", atoms)?;
    dict.set_item("bonds", bonds)?;
    dict.set_item("cell_edges", cell_edges)?;
    dict.set_item("settings", settings)?;
    Ok(dict)
}

/// Builds an ASCII STL string.
#[pyfunction]
#[pyo3(signature = (text, style = "ball-stick".to_string(), scale = 1.0, format_hint = None))]
pub(crate) fn mesh_stl(
    text: &str,
    style: String,
    scale: f64,
    format_hint: Option<&str>,
) -> PyResult<String> {
    let structure = load_structure(text, format_hint)?;
    let settings = ViewSettings {
        style: self::style(&style)?,
        atom_scale: scale,
        ..ViewSettings::default()
    };
    let scene = build_scene(&structure, &settings).map_err(crate::map_error)?;
    Ok(write_stl(&build_mesh(&scene), structure.name()))
}

/// Builds a Wavefront OBJ string.
#[pyfunction]
#[pyo3(signature = (text, style = "ball-stick".to_string(), scale = 1.0, format_hint = None))]
pub(crate) fn mesh_obj(
    text: &str,
    style: String,
    scale: f64,
    format_hint: Option<&str>,
) -> PyResult<String> {
    let structure = load_structure(text, format_hint)?;
    let settings = ViewSettings {
        style: self::style(&style)?,
        atom_scale: scale,
        ..ViewSettings::default()
    };
    let scene = build_scene(&structure, &settings).map_err(crate::map_error)?;
    Ok(write_obj(&build_mesh(&scene), structure.name()))
}

/// Builds a glTF 2.0 JSON string.
#[pyfunction]
#[pyo3(signature = (text, style = "ball-stick".to_string(), scale = 1.0, format_hint = None))]
pub(crate) fn mesh_gltf(
    text: &str,
    style: String,
    scale: f64,
    format_hint: Option<&str>,
) -> PyResult<String> {
    let structure = load_structure(text, format_hint)?;
    let settings = ViewSettings {
        style: self::style(&style)?,
        atom_scale: scale,
        ..ViewSettings::default()
    };
    let scene = build_scene(&structure, &settings).map_err(crate::map_error)?;
    write_gltf(&build_mesh(&scene), structure.name()).map_err(crate::map_error)
}
