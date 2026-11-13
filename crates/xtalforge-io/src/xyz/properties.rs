//! Extended-XYZ property columns (Task 13.2.6).
//!
//! `Properties=species:S:1:pos:R:3:forces:R:3` declares the per-atom columns.
//! Values become site properties with their original names; the comment keys
//! that are not `Lattice` or `Properties` become structure metadata.

#![allow(clippy::indexing_slicing)] // fixed-format parser, lengths checked

use std::collections::BTreeMap;

use serde_json::Value;
use xtalforge_core::elements::{element_by_z, parse_symbol};
use xtalforge_core::{Error, Site};

/// One declared column of the atom table.
#[derive(Debug, Clone)]
struct Field {
    key: String,
    kind: char,
    count: usize,
}

/// A parsed `Properties=` declaration.
#[derive(Debug, Clone)]
pub(crate) struct PropertySpec {
    fields: Vec<Field>,
    position: Option<usize>,
    total: usize,
}

impl PropertySpec {
    /// The default `species:S:1:pos:R:3` layout of a plain XYZ file.
    pub(crate) fn plain() -> Self {
        Self {
            fields: vec![
                Field {
                    key: "species".to_owned(),
                    kind: 'S',
                    count: 1,
                },
                Field {
                    key: "pos".to_owned(),
                    kind: 'R',
                    count: 3,
                },
            ],
            position: Some(1),
            total: 4,
        }
    }
}

/// Parses the `Properties=` token of an extended-XYZ comment line.
///
/// # Errors
///
/// Returns [`Error::InvalidInput`] for an unknown column type or a count of
/// zero.
pub(crate) fn parse_properties(comment: &str) -> Result<Option<PropertySpec>, Error> {
    let Some(start) = comment.find("Properties=") else {
        return Ok(None);
    };
    let rest = &comment[start + "Properties=".len()..];
    let end = rest
        .find(|character: char| character.is_whitespace())
        .unwrap_or(rest.len());
    let declaration = rest[..end].trim_matches('"');
    if declaration.is_empty() {
        return Ok(None);
    }
    let tokens: Vec<&str> = declaration.split(':').collect();
    if !tokens.len().is_multiple_of(3) {
        return Err(Error::InvalidInput(format!(
            "extXYZ Properties needs key:type:count triples, got {declaration:?}"
        )));
    }
    let mut fields = Vec::new();
    let mut offset = 0usize;
    let mut position = None;
    for triple in tokens.chunks(3) {
        let (Some(key), Some(kind), Some(count)) = (triple.first(), triple.get(1), triple.get(2))
        else {
            return Err(Error::InvalidInput(
                "extXYZ Properties triple is incomplete".to_owned(),
            ));
        };
        let count: usize = count
            .parse()
            .map_err(|_| Error::InvalidInput(format!("bad extXYZ property count: {count:?}")))?;
        if count == 0 {
            return Err(Error::InvalidInput(
                "an extXYZ property needs at least one column".to_owned(),
            ));
        }
        let kind = kind
            .chars()
            .next()
            .filter(|character| matches!(character.to_ascii_uppercase(), 'R' | 'I' | 'S' | 'L'))
            .map(|character| character.to_ascii_uppercase())
            .ok_or_else(|| {
                Error::InvalidInput(format!("unknown extXYZ property type: {kind:?}"))
            })?;
        if *key == "pos" {
            if count != 3 {
                return Err(Error::InvalidInput(
                    "the extXYZ pos column needs three values".to_owned(),
                ));
            }
            position = Some(offset);
        }
        fields.push(Field {
            key: (*key).to_owned(),
            kind,
            count,
        });
        offset += count;
    }
    Ok(Some(PropertySpec {
        fields,
        position,
        total: offset,
    }))
}

/// One atom parsed from an extended-XYZ line.
pub(crate) struct AtomRecord {
    /// Element symbol.
    pub element: String,
    /// Cartesian coordinates in Angstrom.
    pub cart: [f64; 3],
    /// Site label from a `label` column, when present.
    pub label: Option<String>,
    /// Non-position columns.
    pub props: BTreeMap<String, Value>,
}

/// Parses one atom line against a [`PropertySpec`].
///
/// # Errors
///
/// Returns [`Error::InvalidInput`] for too few fields, a bad number, a bad
/// logical or an unknown element.
pub(crate) fn parse_atom(fields: &[String], spec: &PropertySpec) -> Result<AtomRecord, Error> {
    if fields.len() < spec.total {
        return Err(Error::InvalidInput(format!(
            "extXYZ atom line has {} values but Properties declares {}",
            fields.len(),
            spec.total
        )));
    }
    let mut offset = 0usize;
    let mut element = String::new();
    let mut cart = [0.0; 3];
    let mut label = None;
    let mut props = BTreeMap::new();
    for field in &spec.fields {
        let slice = fields.get(offset..offset + field.count).ok_or_else(|| {
            Error::InvalidInput("extXYZ atom line is shorter than Properties".to_owned())
        })?;
        offset += field.count;
        if field.key == "species" || field.key == "element" {
            element = parse_element(&slice[0])?;
        } else if field.key == "Z" {
            let z: u8 = slice[0]
                .parse()
                .map_err(|_| Error::InvalidInput("bad extXYZ atomic number".to_owned()))?;
            element = element_by_z(z)?.symbol.clone();
        } else if field.key == "pos" {
            for (axis, value) in slice.iter().enumerate() {
                cart[axis] = parse_real(value, "extXYZ position")?;
            }
        } else {
            let value = parse_value(field.kind, slice)?;
            if field.key == "label"
                && let Some(text) = value.as_str()
            {
                label = Some(text.to_owned());
            }
            props.insert(field.key.clone(), value);
        }
    }
    if element.is_empty() {
        return Err(Error::InvalidInput(
            "extXYZ Properties declares no species column".to_owned(),
        ));
    }
    if spec.position.is_none() {
        return Err(Error::InvalidInput(
            "extXYZ Properties declares no pos column".to_owned(),
        ));
    }
    Ok(AtomRecord {
        element,
        cart,
        label,
        props,
    })
}

fn parse_element(text: &str) -> Result<String, Error> {
    let cleaned = text.trim_matches(|character| character == '"' || character == '\'');
    parse_symbol(cleaned)
}

fn parse_real(text: &str, what: &str) -> Result<f64, Error> {
    text.trim()
        .parse()
        .map_err(|_| Error::InvalidInput(format!("bad {what}: {text:?}")))
}

fn parse_value(kind: char, slice: &[String]) -> Result<Value, Error> {
    if slice.len() == 1 {
        return scalar_value(kind, &slice[0]);
    }
    let mut values = Vec::with_capacity(slice.len());
    for text in slice {
        values.push(scalar_value(kind, text)?);
    }
    Ok(Value::Array(values))
}

fn scalar_value(kind: char, text: &str) -> Result<Value, Error> {
    let text = text.trim();
    match kind {
        'R' => Ok(Value::from(parse_real(text, "extXYZ real")?)),
        'I' => text
            .parse::<i64>()
            .map(Value::from)
            .map_err(|_| Error::InvalidInput(format!("bad extXYZ integer: {text:?}"))),
        'L' => match text.to_ascii_lowercase().as_str() {
            "t" | "true" | "1" => Ok(Value::Bool(true)),
            "f" | "false" | "0" => Ok(Value::Bool(false)),
            _ => Err(Error::InvalidInput(format!("bad extXYZ logical: {text:?}"))),
        },
        _ => Ok(Value::String(
            text.trim_matches(|character| character == '"' || character == '\'')
                .to_owned(),
        )),
    }
}

/// Splits an atom line into fields, keeping quoted strings together.
pub(crate) fn split_fields(line: &str) -> Vec<String> {
    let mut fields = Vec::new();
    let mut current = String::new();
    let mut quote: Option<char> = None;
    for character in line.chars() {
        if let Some(active) = quote {
            if character == active {
                quote = None;
            } else {
                current.push(character);
            }
            continue;
        }
        match character {
            '"' | '\'' => quote = Some(character),
            value if value.is_whitespace() => {
                if !current.is_empty() {
                    fields.push(std::mem::take(&mut current));
                }
            }
            value => current.push(value),
        }
    }
    if !current.is_empty() {
        fields.push(current);
    }
    fields
}

/// Applies the parsed properties to a site.
pub(crate) fn decorate(mut site: Site, record: &AtomRecord) -> Site {
    for (key, value) in &record.props {
        site = site.with_prop(key.clone(), value.clone());
    }
    site
}
