//! The `mesh` command: printable STL from a structure.

use xtalforge::view::{
    DrawStyle, ViewSettings, build_mesh, build_scene, write_gltf, write_obj, write_stl,
};

use crate::format::MeshFormat;
use crate::{Cli, CliError, input, output};

/// Builds an STL mesh and writes it to `out` or standard output.
///
/// # Errors
///
/// Returns [`CliError`] for a read, scene, mesh or write failure.
pub(crate) fn run(
    cli: &Cli,
    path: &str,
    style: &str,
    scale: f64,
    scene_json: bool,
    mesh_format: MeshFormat,
    out: Option<&str>,
) -> Result<(), CliError> {
    let structure = input::read_structure_path(path).map_err(CliError::core)?;
    let draw_style = match style {
        "ball-stick" | "ball_stick" => DrawStyle::BallStick,
        "space-fill" | "space_fill" => DrawStyle::SpaceFill,
        "wireframe" => DrawStyle::Wireframe,
        "sticks" => DrawStyle::Sticks,
        other => {
            return Err(CliError::input(format!(
                "unknown draw style {other:?}; use ball-stick, space-fill, wireframe or sticks"
            )));
        }
    };
    let settings = ViewSettings {
        style: draw_style,
        atom_scale: scale,
        ..ViewSettings::default()
    };
    let scene = build_scene(&structure, &settings).map_err(CliError::core)?;
    if scene_json {
        let json = serde_json::to_string_pretty(&scene).map_err(|error| {
            CliError::input(format!("could not serialize the scene geometry: {error}"))
        })?;
        if let Some(destination) = out.filter(|value| *value != "-") {
            output::write_destination(destination, &json)?;
            if !cli.quiet {
                println!(
                    "wrote {} atoms, {} bonds and {} cell edges to {destination}",
                    scene.atoms.len(),
                    scene.bonds.len(),
                    scene.cell_edges.len()
                );
            }
        } else {
            println!("{json}");
        }
        return Ok(());
    }
    let triangles = build_mesh(&scene);
    let (format, text) = match mesh_format {
        MeshFormat::Stl => ("stl", write_stl(&triangles, structure.name())),
        MeshFormat::Obj => ("obj", write_obj(&triangles, structure.name())),
        MeshFormat::Gltf => (
            "gltf",
            write_gltf(&triangles, structure.name()).map_err(CliError::core)?,
        ),
    };
    if let Some(destination) = out.filter(|value| *value != "-") {
        output::write_destination(destination, &text)?;
        if cli.json {
            println!(
                "{}",
                serde_json::json!({
                    "out": destination,
                    "format": format,
                    "triangles": triangles.len(),
                    "atoms": scene.atoms.len(),
                    "bonds": scene.bonds.len(),
                })
            );
        } else if !cli.quiet {
            println!(
                "wrote {} triangles as {format} ({} atoms, {} bonds) to {destination}",
                triangles.len(),
                scene.atoms.len(),
                scene.bonds.len()
            );
        }
        return Ok(());
    }
    if cli.json {
        println!(
            "{}",
            serde_json::json!({
                "format": format,
                "triangles": triangles.len(),
                "atoms": scene.atoms.len(),
                "bonds": scene.bonds.len(),
                "content": text,
            })
        );
    } else {
        print!("{text}");
    }
    Ok(())
}
