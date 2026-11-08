//! The `convert` command.

use xtalforge::io::registry::{find, losses};
use xtalforge::io::write_structure;

use crate::format::StructureFormat;
use crate::{Cli, CliError, input, output};

/// Converts a structure between CIF, POSCAR and XYZ.
///
/// # Errors
///
/// Returns [`CliError`] for an unreadable input, an unknown format, a writer
/// failure or an unwritable destination.
pub(crate) fn run(
    cli: &Cli,
    path: &str,
    format: StructureFormat,
    out: Option<&str>,
) -> Result<(), CliError> {
    let structure = input::read_structure_path(path).map_err(CliError::core)?;
    let source = std::path::Path::new(path)
        .extension()
        .and_then(std::ffi::OsStr::to_str)
        .and_then(find);
    let target = find(format.name());
    if let (Some(source), Some(target)) = (source, target) {
        let warnings = losses(source.capabilities, target.capabilities);
        if !warnings.is_empty() {
            if cli.strict {
                return Err(CliError::input(format!(
                    "conversion to {} loses {}",
                    target.name,
                    warnings.join(", ")
                )));
            }
            if !cli.quiet {
                eprintln!(
                    "warning: conversion to {} loses {}",
                    target.name,
                    warnings.join(", ")
                );
            }
        }
    }
    let text = write_structure(&structure, format.name()).map_err(CliError::core)?;
    let destination = out.filter(|value| *value != "-");
    if let Some(destination) = destination {
        output::write_destination(destination, &text)?;
        if cli.json {
            let summary = structure.summary().map_err(CliError::core)?;
            println!(
                "{}",
                serde_json::json!({
                    "format": format.name(),
                    "out": destination,
                    "bytes": text.len(),
                    "summary": summary,
                })
            );
        } else if !cli.quiet {
            println!(
                "wrote {} ({} bytes) to {destination}",
                format.name(),
                text.len()
            );
        }
        return Ok(());
    }
    if cli.json {
        println!(
            "{}",
            serde_json::json!({
                "format": format.name(),
                "bytes": text.len(),
                "content": text,
            })
        );
    } else {
        print!("{text}");
    }
    Ok(())
}
