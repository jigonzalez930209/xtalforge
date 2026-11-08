//! Format dispatch through the registry (Task 3.1.2).

use xtalforge_core::{Error, Structure};

use crate::registry;

/// Reads a structure from text with a format hint or content sniffing.
///
/// # Errors
///
/// Returns [`Error::InvalidInput`] for an unknown format or a reader error.
pub fn read_structure(text: &str, format_hint: &str) -> Result<Structure, Error> {
    let hint = format_hint.trim();
    if !hint.is_empty() {
        if let Some(format) = registry::find(hint) {
            return registry::read_with(text, format.name);
        }
        return Err(Error::InvalidInput(format!(
            "unknown structure format: {hint}"
        )));
    }
    if let Some(format) = registry::sniff(text) {
        return registry::read_with(text, format.name);
    }
    Err(Error::InvalidInput(
        "could not detect the structure format; pass the format hint".to_owned(),
    ))
}

/// Writes a structure in a named format.
///
/// # Errors
///
/// Returns [`Error::InvalidInput`] for an unknown format or a writer error.
pub fn write_structure(structure: &Structure, format: &str) -> Result<String, Error> {
    registry::write_with(structure, format)
}
