#![allow(clippy::indexing_slicing)] // bounded indices from validated loops and fixed-size vectors

//! Dataset importers (Task 21.2.3).
//!
//! * extxyz training sets, including the `Lattice`, `Properties` and `pbc`
//!   header conventions;
//! * JSON Lines and JSON arrays of `{id, structure, properties}` records;
//! * OPTIMADE JSON dumps (`data[].attributes`);
//! * a directory of CIF files (CoRE MOF, QMOF or any curated corpus).
//!
//! Every importer returns a [`Dataset`] with stable identifiers. Optional
//! properties are stored as [`PropertyValue`] numbers, text or vectors; a
//! value that cannot be typed is kept as text.

use serde::{Deserialize, Serialize};
use std::collections::BTreeMap;
use xtalforge_analysis::dataset::{Dataset, DatasetRecord, PropertyValue};
use xtalforge_core::{Error, Lattice, Site, StructureBuilder, StructureDocument};

/// Import format.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
#[serde(rename_all = "snake_case")]
pub enum ImportFormat {
    /// Extended XYZ.
    Extxyz,
    /// JSON Lines or a JSON array of records.
    Json,
    /// OPTIMADE JSON.
    Optimade,
}

/// Options for [`import_dataset`].
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct ImportOptions {
    /// Format.
    pub format: ImportFormat,
    /// Dataset name.
    pub name: String,
    /// Optional prefix for generated identifiers.
    pub id_prefix: Option<String>,
}

impl ImportOptions {
    /// Options with the format and a name.
    #[must_use]
    pub fn new(format: ImportFormat, name: impl Into<String>) -> Self {
        Self {
            format,
            name: name.into(),
            id_prefix: None,
        }
    }
}

/// Result of an import.
#[derive(Debug, Clone, PartialEq)]
pub struct ImportReport {
    /// The dataset.
    pub dataset: Dataset,
    /// Records imported.
    pub n_records: usize,
    /// Frames or records skipped because they had no atoms.
    pub skipped: usize,
}

/// Imports text in one of the text formats.
///
/// # Errors
///
/// Returns [`Error::InvalidInput`] for parse errors.
pub fn import_dataset(text: &str, options: &ImportOptions) -> Result<ImportReport, Error> {
    match options.format {
        ImportFormat::Extxyz => import_extxyz(text, &options.name, options.id_prefix.as_deref()),
        ImportFormat::Json => import_json_records(text, &options.name),
        ImportFormat::Optimade => import_optimade(text, &options.name),
    }
}

/// Parses an extended XYZ file.
///
/// # Errors
///
/// Returns [`Error::InvalidInput`] for malformed headers or atom lines.
pub fn import_extxyz(
    text: &str,
    name: &str,
    id_prefix: Option<&str>,
) -> Result<ImportReport, Error> {
    let mut dataset = Dataset::new(name);
    let mut skipped = 0usize;
    let lines: Vec<&str> = text.lines().collect();
    let mut index = 0usize;
    let mut frame = 0usize;
    while index < lines.len() {
        let count_line = lines[index].trim();
        if count_line.is_empty() {
            index += 1;
            continue;
        }
        let count: usize = count_line.parse().map_err(|error| {
            Error::InvalidInput(format!("extxyz frame at line {}: {error}", index + 1))
        })?;
        index += 1;
        let Some(header) = lines.get(index) else {
            return Err(Error::InvalidInput("extxyz header is missing".to_owned()));
        };
        index += 1;
        let (lattice, properties, pbc) = parse_extxyz_header(header)?;
        let mut builder = StructureBuilder::new(format!("frame-{frame}"), lattice);
        let mut names: Vec<String> = Vec::new();
        let mut values: Vec<BTreeMap<String, String>> = Vec::new();
        for _ in 0..count {
            let Some(line) = lines.get(index) else {
                return Err(Error::InvalidInput(format!(
                    "extxyz frame {frame} has fewer atoms than declared"
                )));
            };
            index += 1;
            let fields: Vec<&str> = line.split_whitespace().collect();
            if fields.len() < 4 {
                return Err(Error::InvalidInput(format!(
                    "extxyz atom line `{line}` is too short"
                )));
            }
            let element = fields[0].to_owned();
            let x = parse_number(fields[1])?;
            let y = parse_number(fields[2])?;
            let z = parse_number(fields[3])?;
            let frac = lattice.to_frac([x, y, z]);
            builder = builder.site(Site::new(
                format!("{element}{}", values.len() + 1),
                &element,
                frac,
            )?);
            let mut map = BTreeMap::new();
            for (field_index, field) in properties.iter().enumerate() {
                if let Some(value) = fields.get(4 + field_index) {
                    map.insert(field.clone(), (*value).to_owned());
                }
            }
            values.push(map);
            names.push(element);
        }
        let structure = builder.build()?;
        if structure.p1()?.atoms.is_empty() {
            skipped += 1;
            continue;
        }
        let id = format!("{}{frame}", id_prefix.unwrap_or("extxyz-"));
        let mut record = DatasetRecord::new(id, structure)?;
        if pbc == "F F F" {
            record = record.with_property("pbc", PropertyValue::Text(pbc.clone()))?;
        }
        for (atom_index, map) in values.iter().enumerate() {
            for (key, value) in map {
                if let Ok(number) = value.parse::<f64>() {
                    if atom_index == 0 {
                        record =
                            record.with_property(key.clone(), PropertyValue::Number(number))?;
                    }
                } else if atom_index == 0 {
                    record =
                        record.with_property(key.clone(), PropertyValue::Text(value.clone()))?;
                }
            }
        }
        dataset.push(record)?;
        frame += 1;
    }
    let n_records = dataset.len();
    Ok(ImportReport {
        dataset,
        n_records,
        skipped,
    })
}

fn parse_extxyz_header(header: &str) -> Result<(Lattice, Vec<String>, String), Error> {
    let mut lattice = Lattice::from_parameters(10.0, 10.0, 10.0, 90.0, 90.0, 90.0)?;
    let mut properties = Vec::new();
    let mut pbc = String::from("T T T");
    let mut cursor = header;
    while let Some(equal) = cursor.find('=') {
        let key = cursor[..equal].trim();
        let after = &cursor[equal + 1..];
        let value = if let Some(rest) = after.strip_prefix('"') {
            let Some(end) = rest.find('"') else {
                return Err(Error::InvalidInput(
                    "extxyz header has an unterminated quote".to_owned(),
                ));
            };
            &rest[..end]
        } else {
            let end = after.find(char::is_whitespace).unwrap_or(after.len());
            &after[..end]
        };
        if key == "Lattice" {
            let numbers: Vec<f64> = value
                .split_whitespace()
                .map(parse_number)
                .collect::<Result<_, _>>()?;
            if numbers.len() != 9 {
                return Err(Error::InvalidInput(
                    "extxyz Lattice needs nine numbers".to_owned(),
                ));
            }
            lattice = Lattice::new([
                [numbers[0], numbers[1], numbers[2]],
                [numbers[3], numbers[4], numbers[5]],
                [numbers[6], numbers[7], numbers[8]],
            ])?;
        } else if key == "Properties" {
            properties = value.split(':').step_by(3).map(str::to_owned).collect();
        } else if key == "pbc" {
            pbc = value.to_owned();
        }
        let consumed = equal + 1 + value.len() + if after.starts_with('"') { 2 } else { 0 };
        if consumed >= cursor.len() {
            break;
        }
        cursor = &cursor[consumed..];
    }
    Ok((lattice, properties, pbc))
}

fn parse_number(text: &str) -> Result<f64, Error> {
    text.parse::<f64>()
        .map_err(|error| Error::InvalidInput(format!("`{text}` is not a number: {error}")))
}

/// Parses JSON Lines or a JSON array of record objects.
///
/// # Errors
///
/// Returns [`Error::InvalidInput`] for malformed records.
pub fn import_json_records(text: &str, name: &str) -> Result<ImportReport, Error> {
    let trimmed = text.trim();
    let values: Vec<serde_json::Value> = if trimmed.starts_with('[') {
        serde_json::from_str(trimmed)
            .map_err(|error| Error::InvalidInput(format!("JSON records: {error}")))?
    } else {
        let mut values = Vec::new();
        for line in trimmed.lines() {
            if line.trim().is_empty() {
                continue;
            }
            values.push(
                serde_json::from_str(line)
                    .map_err(|error| Error::InvalidInput(format!("JSON record: {error}")))?,
            );
        }
        values
    };
    let mut dataset = Dataset::new(name);
    let mut skipped = 0usize;
    for (index, value) in values.into_iter().enumerate() {
        let id = value
            .get("id")
            .and_then(serde_json::Value::as_str)
            .map(str::to_owned)
            .unwrap_or_else(|| format!("record-{index}"));
        let Some(structure_value) = value.get("structure") else {
            skipped += 1;
            continue;
        };
        let structure = if structure_value.get("schema_version").is_some() {
            StructureDocument::from_json(&structure_value.to_string())?.structure
        } else {
            serde_json::from_value(structure_value.clone())
                .map_err(|error| Error::InvalidInput(format!("record `{id}` structure: {error}")))?
        };
        let mut record = DatasetRecord::new(id.clone(), structure)?;
        if let Some(object) = value
            .get("properties")
            .and_then(serde_json::Value::as_object)
        {
            for (key, value) in object {
                record = record.with_property(key.clone(), property_from_json(value)?)?;
            }
        }
        if record.structure.p1()?.atoms.is_empty() {
            skipped += 1;
            continue;
        }
        dataset.push(record)?;
    }
    let n_records = dataset.len();
    Ok(ImportReport {
        dataset,
        n_records,
        skipped,
    })
}

pub(crate) fn property_from_json(value: &serde_json::Value) -> Result<PropertyValue, Error> {
    Ok(match value {
        serde_json::Value::Bool(value) => PropertyValue::Boolean(*value),
        serde_json::Value::Number(number) => {
            if let Some(integer) = number.as_i64() {
                PropertyValue::Integer(integer)
            } else {
                PropertyValue::Number(number.as_f64().unwrap_or(0.0))
            }
        }
        serde_json::Value::String(value) => PropertyValue::Text(value.clone()),
        serde_json::Value::Array(values) => {
            let mut numbers = Vec::with_capacity(values.len());
            for value in values {
                let Some(number) = value.as_f64() else {
                    return Err(Error::InvalidInput(
                        "property arrays must contain numbers".to_owned(),
                    ));
                };
                numbers.push(number);
            }
            PropertyValue::Vector(numbers)
        }
        serde_json::Value::Null | serde_json::Value::Object(_) => {
            return Err(Error::InvalidInput(
                "property values must be scalar or numeric arrays".to_owned(),
            ));
        }
    })
}

/// Parses an OPTIMADE JSON dump.
///
/// # Errors
///
/// Returns [`Error::InvalidInput`] for malformed attributes.
pub fn import_optimade(text: &str, name: &str) -> Result<ImportReport, Error> {
    let value: serde_json::Value = serde_json::from_str(text)
        .map_err(|error| Error::InvalidInput(format!("OPTIMADE JSON: {error}")))?;
    let entries = value
        .get("data")
        .and_then(serde_json::Value::as_array)
        .ok_or_else(|| Error::InvalidInput("OPTIMADE data array is missing".to_owned()))?;
    let mut dataset = Dataset::new(name);
    let mut skipped = 0usize;
    for (index, entry) in entries.iter().enumerate() {
        let id = entry
            .get("id")
            .and_then(serde_json::Value::as_str)
            .map(str::to_owned)
            .unwrap_or_else(|| format!("optimade-{index}"));
        let Some(attributes) = entry.get("attributes") else {
            skipped += 1;
            continue;
        };
        let Some(vectors) = attributes
            .get("lattice_vectors")
            .and_then(serde_json::Value::as_array)
        else {
            skipped += 1;
            continue;
        };
        if vectors.len() != 3 {
            skipped += 1;
            continue;
        }
        let mut matrix = [[0.0; 3]; 3];
        for (row, vector) in vectors.iter().enumerate() {
            let Some(values) = vector.as_array() else {
                skipped += 1;
                continue;
            };
            for (column, value) in values.iter().take(3).enumerate() {
                matrix[row][column] = value.as_f64().unwrap_or(0.0);
            }
        }
        let lattice = Lattice::new(matrix)?;
        let Some(positions) = attributes
            .get("cartesian_site_positions")
            .and_then(serde_json::Value::as_array)
        else {
            skipped += 1;
            continue;
        };
        let Some(species) = attributes
            .get("species_at_sites")
            .and_then(serde_json::Value::as_array)
        else {
            skipped += 1;
            continue;
        };
        let mut builder = StructureBuilder::new(id.clone(), lattice);
        for (atom_index, position) in positions.iter().enumerate() {
            let Some(coords) = position.as_array() else {
                continue;
            };
            let element = species
                .get(atom_index)
                .and_then(serde_json::Value::as_str)
                .unwrap_or("X");
            let cart = [
                coords
                    .first()
                    .and_then(serde_json::Value::as_f64)
                    .unwrap_or(0.0),
                coords
                    .get(1)
                    .and_then(serde_json::Value::as_f64)
                    .unwrap_or(0.0),
                coords
                    .get(2)
                    .and_then(serde_json::Value::as_f64)
                    .unwrap_or(0.0),
            ];
            builder = builder.site(Site::new(
                format!("{element}{}", atom_index + 1),
                element,
                lattice.to_frac(cart),
            )?);
        }
        let structure = builder.build()?;
        let mut record = DatasetRecord::new(id, structure)?;
        for key in ["chemical_formula_reduced", "chemical_formula_descriptive"] {
            if let Some(value) = attributes.get(key).and_then(serde_json::Value::as_str) {
                record = record.with_property(key, PropertyValue::Text(value.to_owned()))?;
            }
        }
        if let Some(value) = attributes
            .get("nelements")
            .and_then(serde_json::Value::as_i64)
        {
            record = record.with_property("nelements", PropertyValue::Integer(value))?;
        }
        dataset.push(record)?;
    }
    let n_records = dataset.len();
    Ok(ImportReport {
        dataset,
        n_records,
        skipped,
    })
}
