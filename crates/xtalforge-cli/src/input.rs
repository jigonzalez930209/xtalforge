//! Reads structures from files, standard input or text with a format hint.

use std::sync::OnceLock;

use xtalforge::core::{Error, Structure, StructureDocument};

static SPACEGROUP_OVERRIDE: OnceLock<Option<String>> = OnceLock::new();

/// Installs the `--spacegroup` override for every structure read afterwards.
pub(crate) fn set_spacegroup_override(symbol: Option<String>) {
    let _ = SPACEGROUP_OVERRIDE.set(symbol);
}

fn spacegroup_override() -> Option<&'static str> {
    SPACEGROUP_OVERRIDE.get().and_then(|value| value.as_deref())
}
use xtalforge::io::decode_text;

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
    let structure = if path == "-" {
        detect(text)?
    } else {
        match format_hint(path) {
            None => detect(text)?,
            Some(extension) => read_format(text, &extension)?,
        }
    };
    apply_spacegroup_override(structure)
}

fn apply_spacegroup_override(structure: Structure) -> Result<Structure, Error> {
    let Some(symbol) = spacegroup_override() else {
        return Ok(structure);
    };
    let group = xtalforge::symmetry::SpaceGroup::parse(symbol)?;
    let report = xtalforge::symmetry::set_space_group(
        &structure,
        &group,
        xtalforge::symmetry::SetMode::Reinterpret,
    )?;
    Ok(report.structure)
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

fn format_hint(path: &str) -> Option<String> {
    let lower = path.to_ascii_lowercase();
    if lower.ends_with("geometry.in") {
        return Some("fhi-aims".to_owned());
    }
    if lower.ends_with("vasprun.xml") {
        return Some("vasprun".to_owned());
    }
    extension(path)
}

fn extension(path: &str) -> Option<String> {
    std::path::Path::new(path)
        .extension()
        .and_then(std::ffi::OsStr::to_str)
        .map(str::to_ascii_lowercase)
}

fn read_format(text: &str, format: &str) -> Result<Structure, Error> {
    if format == "json" {
        return structure_from_json(text);
    }
    xtalforge::io::read_structure(text, format)
}

fn detect(text: &str) -> Result<Structure, Error> {
    if text.trim_start().starts_with('{') {
        return structure_from_json(text);
    }
    xtalforge::io::read_structure(text, "")
}
