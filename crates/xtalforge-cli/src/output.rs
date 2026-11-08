//! Writes structures and command output to files or standard output.

use xtalforge::core::Structure;
use xtalforge::io::write_structure;

use crate::Cli;
use crate::CliError;
use crate::format::StructureFormat;

/// Chooses the output format from an explicit flag, the path extension or CIF.
pub(crate) fn resolve_format(
    out: Option<&str>,
    explicit: Option<StructureFormat>,
) -> StructureFormat {
    explicit
        .or_else(|| out.and_then(StructureFormat::from_path))
        .unwrap_or(StructureFormat::Cif)
}

/// Writes a structure to a path or standard output.
///
/// In JSON mode a structure written to standard output is a schema document;
/// when a file is written, the object carries its summary and size instead.
///
/// # Errors
///
/// Returns [`CliError`] for a writer, filesystem or serialization failure.
pub(crate) fn emit_structure(
    cli: &Cli,
    structure: &Structure,
    out: Option<&str>,
    explicit: Option<StructureFormat>,
) -> Result<(), CliError> {
    let format = resolve_format(out, explicit);
    let text = write_structure(structure, format.name()).map_err(CliError::core)?;
    let destination = out.filter(|path| *path != "-");
    if let Some(path) = destination {
        write_destination(path, &text)?;
        if cli.json {
            let summary = structure.summary().map_err(CliError::core)?;
            println!(
                "{}",
                serde_json::json!({
                    "format": format.name(),
                    "out": path,
                    "bytes": text.len(),
                    "summary": summary,
                })
            );
        } else if !cli.quiet {
            println!(
                "wrote {} ({} sites, {} bytes) to {path}",
                format.name(),
                structure.sites().len(),
                text.len()
            );
        }
        return Ok(());
    }
    if cli.json {
        let document = structure.to_document();
        println!("{}", document.to_json().map_err(CliError::core)?);
    } else {
        print!("{text}");
    }
    Ok(())
}

/// Writes text to a filesystem path.
///
/// # Errors
///
/// Returns [`CliError`] with code `io_error` when the file cannot be written.
pub(crate) fn write_destination(path: &str, text: &str) -> Result<(), CliError> {
    std::fs::write(path, text)
        .map_err(|error| CliError::new("io_error", format!("could not write {path}: {error}")))
}
