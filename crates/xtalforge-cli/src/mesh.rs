//! The `mesh` command: printable STL from a structure.

use xtalforge::view::{DrawStyle, ViewSettings, build_mesh, build_scene, write_stl};

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
    let triangles = build_mesh(&scene);
    let stl = write_stl(&triangles, structure.name());
    if let Some(destination) = out.filter(|value| *value != "-") {
        output::write_destination(destination, &stl)?;
        if cli.json {
            println!(
                "{}",
                serde_json::json!({
                    "out": destination,
                    "triangles": triangles.len(),
                    "atoms": scene.atoms.len(),
                    "bonds": scene.bonds.len(),
                })
            );
        } else if !cli.quiet {
            println!(
                "wrote {} triangles ({} atoms, {} bonds) to {destination}",
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
                "triangles": triangles.len(),
                "atoms": scene.atoms.len(),
                "bonds": scene.bonds.len(),
                "stl": stl,
            })
        );
    } else {
        print!("{stl}");
    }
    Ok(())
}
