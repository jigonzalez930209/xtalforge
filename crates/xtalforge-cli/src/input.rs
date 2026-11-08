//! Reads structures from files, standard input or text with a format hint.

use xtalforge::core::{Error, Structure, StructureDocument};
use xtalforge::io::{decode_text, read_cif, read_poscar, read_xyz};

/// A format reader with a uniform signature.
type Reader = fn(&str) -> Result<Structure, Error>;

/// Reads a path, or standard input when `path` is `-`, as a structure.
///
/// # Errors
///
/// Returns [`Error::InvalidInput`] when the file cannot be read, the extension
/// is unknown, or the selected reader rejects the content.
pub(crate) fn read_structure_path(path: &str) -> Result<Structure, Error> {
    let text = read_source(path)?;
    structure_from_text(&text, path)
}

/// Reads decoded text as a structure, using the path extension as a hint.
///
/// # Errors
///
/// Returns the reader errors or [`Error::InvalidInput`] for an unknown format.
pub(crate) fn structure_from_text(text: &str, path: &str) -> Result<Structure, Error> {
    if path == "-" {
        return detect(text);
    }
    match extension(path) {
        None => detect(text),
        Some(extension) => read_format(text, &extension),
    }
}

/// Reads a path, or standard input when `path` is `-`, as decoded text.
///
/// # Errors
///
/// Returns [`Error::InvalidInput`] when the file or stream cannot be read.
pub(crate) fn read_source(path: &str) -> Result<String, Error> {
    if path == "-" {
        let mut buffer = String::new();
        std::io::Read::read_to_string(&mut std::io::stdin(), &mut buffer).map_err(|error| {
            Error::InvalidInput(format!("could not read standard input: {error}"))
        })?;
        return Ok(buffer);
    }
    let bytes = std::fs::read(path)
        .map_err(|error| Error::InvalidInput(format!("could not read {path}: {error}")))?;
    Ok(decode_text(&bytes))
}

/// Parses a bare structure or a versioned `StructureDocument` envelope.
///
/// # Errors
///
/// Returns [`Error::InvalidInput`] for invalid JSON or an unsupported schema.
pub(crate) fn structure_from_json(text: &str) -> Result<Structure, Error> {
    if text.contains("\"schema_version\"") {
        StructureDocument::from_json(text).map(|document| document.structure)
    } else {
        serde_json::from_str(text)
            .map_err(|error| Error::InvalidInput(format!("invalid structure JSON: {error}")))
    }
}

fn extension(path: &str) -> Option<String> {
    std::path::Path::new(path)
        .extension()
        .and_then(std::ffi::OsStr::to_str)
        .map(str::to_ascii_lowercase)
}

fn read_format(text: &str, format: &str) -> Result<Structure, Error> {
    match format {
        "json" => structure_from_json(text),
        "cif" => read_cif(text),
        "poscar" | "contcar" | "vasp" => read_poscar(text),
        "xyz" => read_xyz(text),
        other => Err(Error::InvalidInput(format!(
            "unknown structure format: {other}; expected CIF, POSCAR or XYZ"
        ))),
    }
}

fn detect(text: &str) -> Result<Structure, Error> {
    let readers: [Reader; 3] = [read_cif, read_poscar, read_xyz];
    let mut last = None;
    for reader in readers {
        match reader(text) {
            Ok(structure) => return Ok(structure),
            Err(error) => last = Some(error),
        }
    }
    Err(Error::InvalidInput(format!(
        "could not detect the structure format; expected CIF, POSCAR or XYZ ({})",
        last.map_or_else(|| "no reader was tried".to_owned(), |error| error.message())
    )))
}
