//! Shared helpers for the WebAssembly object surface.

use wasm_bindgen::prelude::*;

use xtalforge::core::{Error, Structure};

/// Sets a property on a JavaScript object, ignoring the (unreachable) failure.
pub(crate) fn set(object: &js_sys::Object, key: &str, value: JsValue) {
    let _ = js_sys::Reflect::set(object, &JsValue::from_str(key), &value);
}

/// Three numbers as a JavaScript array.
pub(crate) fn array3(values: [f64; 3]) -> js_sys::Array {
    let array = js_sys::Array::new();
    for value in values {
        array.push(&JsValue::from_f64(value));
    }
    array
}

/// Six numbers as a JavaScript array.
pub(crate) fn array6(values: [f64; 6]) -> js_sys::Array {
    let array = js_sys::Array::new();
    for value in values {
        array.push(&JsValue::from_f64(value));
    }
    array
}

/// Three integers as a JavaScript array.
pub(crate) fn array_i32(values: [i32; 3]) -> js_sys::Array {
    let array = js_sys::Array::new();
    for value in values {
        array.push(&JsValue::from_f64(f64::from(value)));
    }
    array
}

/// Atom indices as a JavaScript array.
pub(crate) fn array_usize(values: &[usize]) -> js_sys::Array {
    let array = js_sys::Array::new();
    for value in values {
        array.push(&JsValue::from_f64(*value as f64));
    }
    array
}

/// Adds the fields of a structure summary to an existing object.
pub(crate) fn fill_summary(object: &js_sys::Object, structure: &Structure) -> Result<(), JsValue> {
    let summary = structure.summary().map_err(to_js_error)?;
    set(object, "name", JsValue::from_str(&summary.name));
    set(object, "formula", JsValue::from_str(&summary.formula));
    set(object, "hall", JsValue::from_str(&summary.hall));
    set(object, "n_sites", JsValue::from_f64(summary.n_sites as f64));
    set(
        object,
        "n_atoms",
        JsValue::from_f64(structure.p1().map_err(to_js_error)?.n_atoms() as f64),
    );
    set(object, "volume", JsValue::from_f64(summary.volume));
    set(object, "density", JsValue::from_f64(summary.density));
    set(
        object,
        "total_charge",
        JsValue::from_f64(summary.total_charge),
    );
    set(object, "is_p1", JsValue::from_bool(summary.is_p1));
    let warnings = js_sys::Array::new();
    for warning in &summary.warnings {
        let item = js_sys::Object::new();
        set(&item, "code", JsValue::from_str(&warning.code));
        set(&item, "message", JsValue::from_str(&warning.message));
        let atoms = js_sys::Array::new();
        for atom in &warning.atoms {
            atoms.push(&JsValue::from_f64(*atom as f64));
        }
        set(&item, "atoms", atoms.into());
        let fields = js_sys::Array::new();
        for field in &warning.fields {
            fields.push(&JsValue::from_str(field));
        }
        set(&item, "fields", fields.into());
        if let Some(hint) = &warning.hint {
            set(&item, "hint", JsValue::from_str(hint));
        }
        warnings.push(&item.into());
    }
    set(object, "warnings", warnings.into());
    Ok(())
}

/// A summary object for a structure.
pub(crate) fn summary_object(structure: &Structure) -> Result<JsValue, JsValue> {
    let object = js_sys::Object::new();
    fill_summary(&object, structure)?;
    Ok(object.into())
}

/// Parses a structure from text, detecting CIF, POSCAR or XYZ.
///
/// The no-hint exports (`bonds`, `rings`, `groups`, `neighbors`, the measure
/// functions) use this entry point; `read_structure` takes an explicit hint.
pub(crate) fn parse_any(text: &str) -> Result<Structure, Error> {
    let format = detect_format(text)?;
    xtalforge::io::read_structure(text, format)
}

/// Detects the format of a structure text.
fn detect_format(text: &str) -> Result<&'static str, Error> {
    let trimmed = text.trim_start_matches('\u{feff}').trim_start();
    if trimmed.starts_with("data_")
        || text.contains("_atom_site_")
        || text.contains("_cell_length_")
    {
        return Ok("cif");
    }
    let lines: Vec<&str> = trimmed.lines().collect();
    if lines
        .first()
        .is_some_and(|line| line.trim().parse::<usize>().is_ok())
        && looks_like_xyz(&lines)
    {
        return Ok("xyz");
    }
    if looks_like_poscar(text) {
        return Ok("poscar");
    }
    Err(Error::InvalidInput(
        "could not detect the structure format; pass the format hint to read_structure".to_owned(),
    ))
}

/// Whether the first atom line of an XYZ file is well formed.
fn looks_like_xyz(lines: &[&str]) -> bool {
    let Some(line) = lines.get(2) else {
        return false;
    };
    let mut fields = line.split_whitespace();
    let Some(element) = fields.next() else {
        return false;
    };
    element.chars().next().is_some_and(char::is_alphabetic) && fields.count() >= 3
}

/// Whether the text starts like a POSCAR with a species line.
fn looks_like_poscar(text: &str) -> bool {
    let lines: Vec<&str> = text.lines().collect();
    if lines.len() < 8 {
        return false;
    }
    if lines
        .get(1)
        .is_none_or(|line| line.trim().parse::<f64>().is_err())
    {
        return false;
    }
    let lattice_ok = lines.get(2..5).is_some_and(|rows| {
        rows.iter().all(|row| {
            row.split_whitespace()
                .filter(|value| value.parse::<f64>().is_ok())
                .count()
                >= 3
        })
    });
    let species_ok = lines.get(5).is_some_and(|line| {
        line.split_whitespace()
            .next()
            .is_some_and(|species| species.chars().any(char::is_alphabetic))
    });
    let counts_ok = lines.get(6).is_some_and(|line| {
        line.split_whitespace()
            .next()
            .is_some_and(|count| count.parse::<usize>().is_ok())
    });
    lattice_ok && species_ok && counts_ok
}

/// A JavaScript `Error` with the `invalid_input` code, like a core
/// [`Error::InvalidInput`].
pub(crate) fn invalid_input(message: &str) -> JsValue {
    coded_error("invalid_input", &format!("invalid input: {message}"))
}

/// A JavaScript `Error` with an explicit stable `code`.
pub(crate) fn coded_error(code: &str, message: &str) -> JsValue {
    let js_error = js_sys::Error::new(message);
    let _ = js_sys::Reflect::set(
        &js_error,
        &JsValue::from_str("code"),
        &JsValue::from_str(code),
    );
    js_error.into()
}

/// Converts a core error into a JavaScript `Error` with a stable `code`.
pub(crate) fn to_js_error(error: Error) -> JsValue {
    let js_error = js_sys::Error::new(&error.message());
    let _ = js_sys::Reflect::set(
        &js_error,
        &JsValue::from_str("code"),
        &JsValue::from_str(error.code()),
    );
    let _ = js_sys::Reflect::set(
        &js_error,
        &JsValue::from_str("category"),
        &JsValue::from_str(error.category()),
    );
    js_error.into()
}
