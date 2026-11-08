//! View and STL bindings (Task 3.16.2).

use pyo3::prelude::*;
use pyo3::types::PyDict;

use xtalforge::view::{DrawStyle, ViewSettings, build_mesh, build_scene, write_stl};

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
