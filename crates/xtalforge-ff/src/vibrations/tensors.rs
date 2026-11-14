//! Simple tensor files for Born effective charges and Raman tensors.
//!
//! Two formats are accepted, selected by the first non-space character:
//!
//! * JSON: an object `{"comment": "...", "natoms": N, "tensors": [...]}` or a
//!   bare array of `3 x 3` arrays.  `tensors` may also be spelled `blocks`.
//!   `comment` and `natoms` are optional; when `natoms` is present it must
//!   equal the number of tensors.  Unknown keys are ignored so an
//!   engine-specific export can carry extra metadata.
//! * Text: one tensor per line, either nine numbers in row-major order or a
//!   one-based atom index followed by nine numbers.  Blank lines and `#`
//!   comments are ignored; the first full-line comment becomes the comment.
//!   `D` exponents (`1.0D-02`, Fortran and Quantum ESPRESSO style) are accepted.
//!
//! This is a deliberately small interchange format, not a full `OUTCAR` or
//! Quantum ESPRESSO output parser.  Task 19.4.2 asks for a parser for a simple
//! text/JSON tensor file; producing that file from an engine log stays outside
//! this crate.

use serde::Deserialize;

use crate::Error;

/// Tensors read from a text or JSON tensor file.
#[derive(Debug, Clone, PartialEq)]
pub struct TensorFile {
    /// One `3 x 3` tensor per atom or per mode, row-major.
    pub tensors: Vec<[[f64; 3]; 3]>,
    /// Free-form comment, when the file carries one.
    pub comment: Option<String>,
    /// Declared tensor count, when the file carries one.
    pub natoms: Option<usize>,
}

#[derive(Deserialize)]
struct JsonDocument {
    #[serde(default)]
    comment: Option<String>,
    #[serde(default)]
    natoms: Option<usize>,
    #[serde(default)]
    tensors: Option<Vec<[[f64; 3]; 3]>>,
    #[serde(default)]
    blocks: Option<Vec<[[f64; 3]; 3]>>,
}

/// Parses a text or JSON tensor file.
///
/// The accepted formats are documented on this module.  A file with no tensor
/// is rejected, and every component must be finite.
///
/// # Errors
///
/// Returns [`Error::Json`] for malformed JSON and [`Error::InvalidInput`] for an
/// empty file, a missing tensor list, an inconsistent declared count, a wrong
/// column count, a non-numeric component or a non-finite component.
pub fn parse_tensor_file(text: &str) -> Result<TensorFile, Error> {
    let trimmed = text.trim();
    if trimmed.is_empty() {
        return Err(Error::InvalidInput(
            "tensor file is empty; expected 3 x 3 tensors".to_owned(),
        ));
    }
    if trimmed.starts_with('{') {
        parse_json_object(trimmed)
    } else if trimmed.starts_with('[') {
        let tensors = serde_json::from_str::<Vec<[[f64; 3]; 3]>>(trimmed)
            .map_err(|error| Error::Json(error.to_string()))?;
        finish(TensorFile {
            tensors,
            comment: None,
            natoms: None,
        })
    } else {
        parse_text(text)
    }
}

fn parse_json_object(text: &str) -> Result<TensorFile, Error> {
    let document = serde_json::from_str::<JsonDocument>(text)
        .map_err(|error| Error::Json(error.to_string()))?;
    let tensors = match (document.tensors, document.blocks) {
        (Some(tensors), _) => tensors,
        (None, Some(blocks)) => blocks,
        (None, None) => {
            return Err(Error::InvalidInput(
                "JSON tensor file has no `tensors` array".to_owned(),
            ));
        }
    };
    finish(TensorFile {
        tensors,
        comment: document.comment,
        natoms: document.natoms,
    })
}

fn finish(file: TensorFile) -> Result<TensorFile, Error> {
    if let Some(declared) = file.natoms
        && declared != file.tensors.len()
    {
        return Err(Error::InvalidInput(format!(
            "tensor file declares {declared} tensors but stores {}",
            file.tensors.len()
        )));
    }
    if file.tensors.is_empty() {
        return Err(Error::InvalidInput(
            "tensor file contains no 3 x 3 tensor".to_owned(),
        ));
    }
    for (index, tensor) in file.tensors.iter().enumerate() {
        if tensor.iter().flatten().any(|value| !value.is_finite()) {
            return Err(Error::InvalidInput(format!(
                "tensor {index} has a non-finite component"
            )));
        }
    }
    Ok(file)
}

fn parse_text(text: &str) -> Result<TensorFile, Error> {
    let mut tensors = Vec::new();
    let mut comment = None;
    for (line_index, raw) in text.lines().enumerate() {
        if let Some((before, after)) = raw.split_once('#')
            && before.trim().is_empty()
            && comment.is_none()
        {
            comment = Some(after.trim().to_owned());
        }
        let content = raw.split_once('#').map_or(raw, |(before, _)| before).trim();
        if content.is_empty() {
            continue;
        }
        let tokens = content.split_whitespace().collect::<Vec<_>>();
        let start = match tokens.len() {
            9 => 0,
            10 => {
                let index = tokens
                    .first()
                    .ok_or_else(|| {
                        Error::InvalidInput(format!("line {} has no atom index", line_index + 1))
                    })?
                    .parse::<usize>()
                    .map_err(|_| {
                        Error::InvalidInput(format!(
                            "line {}: tensor index must be a positive integer",
                            line_index + 1
                        ))
                    })?;
                if index != tensors.len() + 1 {
                    return Err(Error::InvalidInput(format!(
                        "line {}: expected tensor index {}, got {index}",
                        line_index + 1,
                        tensors.len() + 1
                    )));
                }
                1
            }
            other => {
                return Err(Error::InvalidInput(format!(
                    "line {}: expected 9 numbers or an index plus 9 numbers, got {other} fields",
                    line_index + 1
                )));
            }
        };
        let mut tensor = [[0.0; 3]; 3];
        for (row, values) in tensor.iter_mut().enumerate() {
            for (column, slot) in values.iter_mut().enumerate() {
                let position = start + 3 * row + column;
                let token = tokens.get(position).ok_or_else(|| {
                    Error::InvalidInput(format!(
                        "line {} ends before tensor component {position}",
                        line_index + 1
                    ))
                })?;
                *slot = number(token, line_index + 1)?;
            }
        }
        tensors.push(tensor);
    }
    finish(TensorFile {
        tensors,
        comment,
        natoms: None,
    })
}

fn number(token: &str, line: usize) -> Result<f64, Error> {
    let normalised = token.replace(['D', 'd'], "E");
    let value = normalised
        .parse::<f64>()
        .map_err(|_| Error::InvalidInput(format!("line {line}: `{token}` is not a number")))?;
    if !value.is_finite() {
        return Err(Error::InvalidInput(format!(
            "line {line}: `{token}` is not finite"
        )));
    }
    Ok(value)
}
