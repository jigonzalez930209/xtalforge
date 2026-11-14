//! Helpers shared by the engine writers.

use serde_json::Value;
use xtalforge_core::{P1Atom, Site, Structure, StructureBuilder};

use crate::Error;

/// Returns the P1 expansion with labels and occupancies preserved.
pub(crate) fn p1_snapshot(structure: &Structure) -> Result<Structure, Error> {
    let p1 = structure.p1()?;
    let mut builder = StructureBuilder::new("p1", p1.lattice);
    for atom in &p1.atoms {
        let site = Site::new(atom.label.clone(), &atom.element, atom.frac)?
            .with_occupancy(atom.occupancy)?;
        builder = builder.site(site);
    }
    Ok(builder.build()?)
}

/// Unique symbols in first-appearance order.
pub(crate) fn species(atoms: &[P1Atom]) -> Vec<String> {
    let mut symbols: Vec<String> = Vec::with_capacity(atoms.len());
    for atom in atoms {
        if !symbols.iter().any(|symbol| symbol == &atom.element) {
            symbols.push(atom.element.clone());
        }
    }
    symbols
}

/// Replaces everything outside `[A-Za-z0-9_]` with `_`, falling back to a
/// fixed name when nothing remains.
pub(crate) fn identifier(text: &str, fallback: &str) -> String {
    let cleaned: String = text
        .chars()
        .map(|character| {
            if character.is_ascii_alphanumeric() || character == '_' {
                character
            } else {
                '_'
            }
        })
        .collect();
    if cleaned.is_empty() {
        fallback.to_owned()
    } else {
        cleaned
    }
}

/// Rejects whitespace, control characters and shell punctuation in names that
/// are written as bare engine tokens.
pub(crate) fn token<'a>(value: &'a str, what: &str) -> Result<&'a str, Error> {
    let valid = !value.is_empty()
        && value.chars().all(|character| {
            character.is_ascii_alphanumeric()
                || matches!(character, '_' | '-' | '.' | '+' | '(' | ')')
        });
    if valid {
        Ok(value)
    } else {
        Err(Error::Job(format!("invalid {what} token {value:?}")))
    }
}

/// Validates a relative path written as a bare engine token.
pub(crate) fn relative_path<'a>(value: &'a str, what: &str) -> Result<&'a str, Error> {
    let valid = !value.is_empty()
        && !value.starts_with('/')
        && !value.contains("..")
        && !value.contains('\\')
        && value.chars().all(|character| {
            character.is_ascii_alphanumeric()
                || matches!(character, '_' | '-' | '.' | '/' | '+' | '(' | ')')
        });
    if valid {
        Ok(value)
    } else {
        Err(Error::Job(format!("invalid {what} path {value:?}")))
    }
}

/// Returns a site property, then the matching structure metadata value.
///
/// Metadata may be a scalar shared by every atom or an array indexed by the
/// P1 atom order.
pub(crate) fn site_or_metadata<'a>(
    structure: &'a Structure,
    atom: &P1Atom,
    index: usize,
    key: &str,
) -> Option<&'a Value> {
    let site = structure.sites().get(atom.site_idx)?;
    if let Some(value) = site.props().get(key) {
        return Some(value);
    }
    let value = structure.metadata().get(key)?;
    match value {
        Value::Array(values) => values.get(index),
        scalar => Some(scalar),
    }
}

/// Reads a finite number from a JSON value.
pub(crate) fn number(value: &Value) -> Option<f64> {
    value.as_f64().filter(|value| value.is_finite())
}

/// Maps an orbital name to the VASP `LDAUL` code.
pub(crate) fn orbital_code(name: &str) -> Result<i64, Error> {
    match name.trim().to_ascii_lowercase().as_str() {
        "s" => Ok(0),
        "p" => Ok(1),
        "d" => Ok(2),
        "f" => Ok(3),
        other => Err(Error::Job(format!(
            "unknown LDAU orbital {other:?}; use s, p, d or f"
        ))),
    }
}

/// Maps a `xtalforge_ff` error into a preset-level error.
pub(crate) fn force_field_error(error: xtalforge_ff::Error) -> Error {
    Error::Job(format!("force field: {error}"))
}
