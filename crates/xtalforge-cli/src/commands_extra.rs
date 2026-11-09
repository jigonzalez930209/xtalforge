//! The `formats`, `params`, `project`, `trajectory` and `script` commands.

use serde_json::json;
use xtalforge::commands::{Document, MoveSites, Origin, SetElement};

use crate::format::StructureFormat;
use crate::{Cli, CliError, input, output};

/// Lists the format registry with its capabilities.
pub(crate) fn formats(cli: &Cli) -> Result<(), CliError> {
    let entries: Vec<serde_json::Value> = xtalforge::io::registry::FORMATS
        .iter()
        .map(|format| {
            json!({
                "name": format.name,
                "extensions": format.extensions,
                "read": format.read.is_some(),
                "write": format.write.is_some(),
                "cell": format.capabilities.cell,
                "symmetry": format.capabilities.symmetry,
                "occupancy": format.capabilities.occupancy,
                "charges": format.capabilities.charges,
                "bonds": format.capabilities.bonds,
                "labels": format.capabilities.labels,
                "frames": format.capabilities.frames,
            })
        })
        .collect();
    if cli.json {
        println!("{}", json!(entries));
    } else {
        for format in xtalforge::io::registry::FORMATS {
            println!(
                "{:<12} {:<24} read={} write={}",
                format.name,
                format.extensions.join(" "),
                format.read.is_some(),
                format.write.is_some()
            );
        }
    }
    Ok(())
}

/// Lists parameter presets or prints one schema.
pub(crate) fn params(cli: &Cli, name: Option<&str>) -> Result<(), CliError> {
    match name {
        Some(name) => {
            let set = xtalforge::core::preset(name).map_err(CliError::core)?;
            let schema = set.schema();
            println!(
                "{}",
                serde_json::to_string_pretty(&schema).map_err(json_error)?
            );
        }
        None => {
            if cli.json {
                println!("{}", json!(xtalforge::core::PRESET_NAMES));
            } else {
                for name in xtalforge::core::PRESET_NAMES {
                    println!("{name}");
                }
            }
        }
    }
    Ok(())
}

/// Packs a structure into a `.xtalproj` archive.
pub(crate) fn project_pack(cli: &Cli, path: &str, out: &str) -> Result<(), CliError> {
    let structure = input::read_structure_path(path).map_err(CliError::core)?;
    let graph = structure.bond_graph().map_err(CliError::core)?;
    let bonds = json!(
        graph
            .bonds
            .iter()
            .map(|bond| json!({
                "i": bond.i,
                "j": bond.j,
                "image": bond.image,
                "distance": bond.distance,
                "order": bond.order,
                "explicit": bond.explicit,
            }))
            .collect::<Vec<_>>()
    );
    let sites = json!(
        structure
            .sites()
            .iter()
            .map(|site| json!({
                "label": site.label(),
                "element": site.element(),
                "occupancy": site.occupancy(),
                "charge": site.charge(),
                "props": site.props(),
            }))
            .collect::<Vec<_>>()
    );
    let parts = xtalforge::io::project::ProjectParts::from_structure(
        &structure,
        None,
        Some(bonds),
        Some(sites),
        None,
        None,
    )
    .map_err(CliError::core)?;
    let bytes = xtalforge::io::project::write_project(&parts).map_err(CliError::core)?;
    output::write_bytes(out, &bytes)?;
    if cli.json {
        println!("{}", json!({"out": out, "bytes": bytes.len()}));
    } else if !cli.quiet {
        println!("wrote project {out} ({} bytes)", bytes.len());
    }
    Ok(())
}

/// Prints the parts of a project archive.
pub(crate) fn project_info(cli: &Cli, path: &str) -> Result<(), CliError> {
    let bytes = std::fs::read(path)
        .map_err(|error| CliError::new("io_error", format!("could not read {path}: {error}")))?;
    let parts = xtalforge::io::project::read_project(&bytes).map_err(CliError::core)?;
    let structure = parts.structure().map_err(CliError::core)?;
    let payload = json!({
        "format_version": xtalforge::io::project::FORMAT_VERSION,
        "name": structure.name(),
        "formula": structure.formula(),
        "n_sites": structure.sites().len(),
        "cell": parts.cell,
        "bonds": parts.bonds.is_some(),
        "sites": parts.sites.is_some(),
        "view": parts.view.is_some(),
        "session": parts.session.is_some(),
    });
    if cli.json {
        println!("{payload}");
    } else {
        println!(
            "{} ({}), {} sites; parts: bonds={} sites={} view={} session={}",
            structure.name(),
            structure.formula(),
            structure.sites().len(),
            parts.bonds.is_some(),
            parts.sites.is_some(),
            parts.view.is_some(),
            parts.session.is_some()
        );
    }
    Ok(())
}

/// Extracts a project into a directory.
pub(crate) fn project_extract(cli: &Cli, path: &str, directory: &str) -> Result<(), CliError> {
    let bytes = std::fs::read(path)
        .map_err(|error| CliError::new("io_error", format!("could not read {path}: {error}")))?;
    let parts = xtalforge::io::project::read_project(&bytes).map_err(CliError::core)?;
    let root = std::path::Path::new(directory);
    std::fs::create_dir_all(root).map_err(|error| {
        CliError::new("io_error", format!("could not create {directory}: {error}"))
    })?;
    output::write_bytes(
        &root.join("structure.cif").to_string_lossy(),
        parts.structure_cif.as_bytes(),
    )?;
    for (name, part) in [
        ("cell.json", parts.cell),
        ("bonds.json", parts.bonds),
        ("sites.json", parts.sites),
        ("view.json", parts.view),
        ("session.json", parts.session),
    ] {
        if let Some(value) = part {
            let text = serde_json::to_string_pretty(&value).map_err(json_error)?;
            output::write_bytes(&root.join(name).to_string_lossy(), text.as_bytes())?;
        }
    }
    if cli.json {
        println!("{}", json!({"out": directory}));
    } else if !cli.quiet {
        println!("extracted project into {directory}");
    }
    Ok(())
}

/// Prints trajectory information or one frame as a structure.
pub(crate) fn trajectory(
    cli: &Cli,
    path: &str,
    frame: Option<usize>,
    out: Option<&str>,
    format: Option<StructureFormat>,
) -> Result<(), CliError> {
    let text = input::read_source(path).map_err(CliError::core)?;
    if let Ok(index) = xtalforge::io::XyzIndex::new(&text) {
        return match frame {
            Some(at) => {
                let structure = index.structure(at).map_err(CliError::core)?;
                output::emit_structure(cli, &structure, out, format)
            }
            None => {
                let payload = json!({
                    "frames": index.len(),
                    "elements": index.elements(),
                    "lattice": index.lattice(),
                    "energies": index.energies(),
                });
                if cli.json {
                    println!("{payload}");
                } else {
                    println!("{} frames, {} atoms", index.len(), index.elements().len());
                }
                Ok(())
            }
        };
    }
    let trajectory = xtalforge::io::read_trajectory(&text).map_err(CliError::core)?;
    match frame {
        Some(at) => {
            let structure = trajectory.structure(at).map_err(CliError::core)?;
            output::emit_structure(cli, &structure, out, format)
        }
        None => {
            let payload = json!({
                "frames": trajectory.frames.len(),
                "elements": trajectory.elements,
                "lattice": trajectory.lattice,
                "energies": trajectory.frames.iter().map(|frame| frame.energy).collect::<Vec<_>>(),
            });
            if cli.json {
                println!("{payload}");
            } else {
                println!(
                    "{} frames, {} atoms",
                    trajectory.frames.len(),
                    trajectory.elements.len()
                );
            }
            Ok(())
        }
    }
}

/// Applies a JSON command list to a structure.
///
/// The script is a JSON object with `structure` (inline structure text) and
/// `commands`; each command is `{"type": "move", "sites": [i], "delta":
/// [x, y, z]}` or `{"type": "set_element", "site": i, "element": "N"}`.
pub(crate) fn script(cli: &Cli, path: &str) -> Result<(), CliError> {
    let text = input::read_source(path).map_err(CliError::core)?;
    let script: serde_json::Value = serde_json::from_str(&text)
        .map_err(|error| CliError::input(format!("invalid script: {error}")))?;
    let structure_text = script
        .get("structure")
        .and_then(serde_json::Value::as_str)
        .ok_or_else(|| CliError::input("script has no structure"))?;
    let structure = xtalforge::io::read_structure(structure_text, "").map_err(CliError::core)?;
    let mut document = Document::new(structure);
    let commands = script
        .get("commands")
        .and_then(serde_json::Value::as_array)
        .ok_or_else(|| CliError::input("script has no commands"))?;
    for (index, command) in commands.iter().enumerate() {
        let kind = command["type"]
            .as_str()
            .ok_or_else(|| CliError::input(format!("command {index} has no type")))?;
        match kind {
            "move" => {
                let sites: Vec<usize> = serde_json::from_value(command["sites"].clone())
                    .map_err(|error| CliError::input(format!("command {index}: {error}")))?;
                let delta: [f64; 3] = serde_json::from_value(command["delta"].clone())
                    .map_err(|error| CliError::input(format!("command {index}: {error}")))?;
                let targets: Vec<(usize, [f64; 3])> = sites
                    .iter()
                    .map(|site| {
                        let atom = document.structure().sites().get(*site).ok_or_else(|| {
                            CliError::input(format!("command {index}: site {site} is out of range"))
                        })?;
                        let frac = atom.frac();
                        Ok((
                            *site,
                            [frac[0] + delta[0], frac[1] + delta[1], frac[2] + delta[2]],
                        ))
                    })
                    .collect::<Result<Vec<_>, CliError>>()?;
                document
                    .apply(Box::new(MoveSites::new(targets)), Origin::Script)
                    .map_err(CliError::core)?;
            }
            "set_element" => {
                let site: usize = serde_json::from_value(command["site"].clone())
                    .map_err(|error| CliError::input(format!("command {index}: {error}")))?;
                let element = command["element"]
                    .as_str()
                    .ok_or_else(|| CliError::input(format!("command {index} has no element")))?;
                let set = SetElement::new(vec![site], element).map_err(CliError::core)?;
                document
                    .apply(Box::new(set), Origin::Script)
                    .map_err(CliError::core)?;
            }
            other => {
                return Err(CliError::input(format!(
                    "command {index}: unknown command {other:?}"
                )));
            }
        }
    }
    output::emit_structure(cli, document.structure(), None, None)
}

fn json_error(error: serde_json::Error) -> CliError {
    CliError::new("json_error", format!("could not serialize: {error}"))
}
