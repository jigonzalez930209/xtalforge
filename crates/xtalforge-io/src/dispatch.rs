//! Minimal format dispatch for Phase 1 (Task 1.16.3).

use xtalforge_core::{Error, Structure};

use crate::{read_cif, read_poscar, read_xyz, write_cif, write_poscar, write_xyz};

/// Reads a structure from text with a format hint.
///
/// The hint may be a format name (`cif`, `poscar`, `xyz`) or a file extension.
///
/// # Errors
///
/// Returns [`Error::InvalidInput`] for an unknown format or a reader error.
pub fn read_structure(text: &str, format_hint: &str) -> Result<Structure, Error> {
    let format = normalized(format_hint);
    match format.as_str() {
        "cif" => read_cif(text),
        "poscar" | "contcar" | "vasp" => read_poscar(text),
        "xyz" => read_xyz(text),
        other => Err(Error::InvalidInput(format!(
            "unknown structure format: {other}"
        ))),
    }
}

/// Writes a structure in a named format.
///
/// # Errors
///
/// Returns [`Error::InvalidInput`] for an unknown format or a writer error.
pub fn write_structure(structure: &Structure, format: &str) -> Result<String, Error> {
    let format = normalized(format);
    match format.as_str() {
        "cif" => write_cif(structure),
        "poscar" | "contcar" | "vasp" => write_poscar(structure),
        "xyz" => write_xyz(structure),
        other => Err(Error::InvalidInput(format!(
            "unknown structure format: {other}"
        ))),
    }
}

fn normalized(format: &str) -> String {
    let lower = format.trim().to_ascii_lowercase();
    match lower.as_str() {
        ".cif" => "cif".to_owned(),
        ".vasp" | ".poscar" | ".contcar" => "poscar".to_owned(),
        ".xyz" => "xyz".to_owned(),
        _ => lower,
    }
}
