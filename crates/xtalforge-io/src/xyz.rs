//! XYZ and extended XYZ reader and writer (Tasks 1.16.2 and 13.2.6).
//!
//! Plain XYZ keeps its historical behaviour: a padded P1 box and a label made
//! from the element and index. Extended XYZ adds `Lattice=`, a `Properties=`
//! declaration and arbitrary `key=value` comment pairs. Property columns are
//! kept as site properties and comment pairs as structure metadata, so a
//! write/read round trip preserves both.

#![allow(clippy::indexing_slicing)] // fixed-format parser, lengths checked

use std::collections::{BTreeMap, BTreeSet};

use serde_json::Value;
use xtalforge_core::{
    Error, Lattice, PROVENANCE_METADATA_KEY, Provenance, Site, Structure, StructureBuilder,
};

use crate::molecule::{MoleculeAtom, molecule_box};
use crate::provenance::{XYZ_PROVENANCE_KEY, provenance_from_xyz_info, xyz_provenance_pair};
use crate::trajectory::parse_info;

mod properties;

use properties::{AtomRecord, PropertySpec, parse_properties};

/// Reads plain or extended XYZ.
///
/// `Lattice="..."` makes the file periodic; without it the atoms become a
/// molecule in a padded P1 box.
///
/// # Errors
///
/// Returns [`Error::InvalidInput`] for a bad atom count, a malformed line, an
/// invalid lattice or a property value that does not match its declared type.
pub fn read_xyz(text: &str) -> Result<Structure, Error> {
    let mut lines = text.lines();
    let count: usize = lines
        .next()
        .and_then(|line| line.trim().parse().ok())
        .ok_or_else(|| Error::InvalidInput("XYZ atom count is missing".to_owned()))?;
    let comment = lines.next().unwrap_or("").trim();
    let lattice = parse_lattice(comment)?;
    let spec = parse_properties(comment)?.unwrap_or_else(PropertySpec::plain);
    let info = parse_info(comment);
    let mut records = Vec::with_capacity(count);
    for _ in 0..count {
        let line = lines
            .next()
            .ok_or_else(|| Error::InvalidInput("XYZ ended before all atoms".to_owned()))?;
        let fields = properties::split_fields(line);
        records.push(properties::parse_atom(&fields, &spec)?);
    }
    let name = if comment.is_empty() { "xyz" } else { comment };
    match lattice {
        Some(lattice) => {
            let mut builder = StructureBuilder::new(name, lattice);
            for (index, record) in records.iter().enumerate() {
                let label = record
                    .label
                    .clone()
                    .unwrap_or_else(|| format!("{}{}", record.element, index + 1));
                let site = Site::new(label, &record.element, lattice.to_frac(record.cart))?;
                builder = builder.site(properties::decorate(site, record));
            }
            let builder = apply_info(builder, &info)?;
            builder.build()
        }
        None => molecule_from_records(name, &records, &info),
    }
}

/// Copies comment pairs into metadata and attaches an embedded provenance
/// record.
fn apply_info(
    mut builder: StructureBuilder,
    info: &BTreeMap<String, String>,
) -> Result<StructureBuilder, Error> {
    for (key, value) in info {
        if key == "Lattice" || key == "Properties" || key == XYZ_PROVENANCE_KEY {
            continue;
        }
        builder = builder.metadata(key.clone(), Value::String(value.clone()));
    }
    if let Some(record) = provenance_from_xyz_info(info)? {
        builder = builder.metadata(PROVENANCE_METADATA_KEY, Value::String(record.to_json()?));
    }
    Ok(builder)
}

fn molecule_from_records(
    name: &str,
    records: &[AtomRecord],
    info: &BTreeMap<String, String>,
) -> Result<Structure, Error> {
    let atoms: Vec<MoleculeAtom> = records
        .iter()
        .enumerate()
        .map(|(index, record)| MoleculeAtom {
            label: record
                .label
                .clone()
                .unwrap_or_else(|| format!("{}{}", record.element, index + 1)),
            element: record.element.clone(),
            cart: record.cart,
        })
        .collect();
    let base = molecule_box(name, &atoms, &[])?;
    let mut builder = StructureBuilder::new(name, *base.lattice());
    for (site, record) in base.sites().iter().zip(records) {
        builder = builder.site(properties::decorate(site.clone(), record));
    }
    let builder = apply_info(builder, info)?;
    builder.build()
}

/// Writes an extended XYZ with `Lattice`, `Properties` and scalar metadata.
///
/// Per-atom properties shared by every site become extra `Properties` columns.
/// Scalar structure metadata becomes comment pairs. A structure carrying a
/// [`Provenance`] record writes it as the base64 `xtalforge_provenance` pair.
/// A structure with no extra properties writes the historical
/// `species:S:1:pos:R:3` header.
///
/// # Errors
///
/// Returns an error if the structure cannot be expanded, or if two sites name
/// the same property with incompatible types.
pub fn write_xyz(structure: &Structure) -> Result<String, Error> {
    let p1 = structure.p1()?;
    let columns = property_columns(structure)?;
    let lattice = structure.lattice().matrix();
    let lattice_text = lattice
        .iter()
        .flat_map(|row| row.iter())
        .map(|value| format!("{value:.10}"))
        .collect::<Vec<_>>()
        .join(" ");
    let mut header = format!("Lattice=\"{lattice_text}\"");
    let mut declaration = String::from("species:S:1:pos:R:3");
    for column in &columns {
        declaration.push_str(&format!(":{}:{}:{}", column.key, column.kind, column.count));
    }
    header.push_str(&format!(" Properties={declaration}"));
    if let Some(record) = Provenance::from_structure(structure)? {
        header.push(' ');
        header.push_str(&xyz_provenance_pair(&record)?);
    }
    for (key, value) in structure.metadata() {
        if key.contains("warning")
            || !is_scalar(value)
            || key == "extxyz_info"
            || key == PROVENANCE_METADATA_KEY
            || key == XYZ_PROVENANCE_KEY
        {
            continue;
        }
        let Some(text) = scalar_text(value) else {
            continue;
        };
        header.push(' ');
        header.push_str(key);
        header.push('=');
        if text.contains(char::is_whitespace) || text.contains('"') {
            if text.contains('"') {
                continue;
            }
            header.push('"');
            header.push_str(&text);
            header.push('"');
        } else {
            header.push_str(&text);
        }
    }
    let mut out = format!("{}\n{header}\n", p1.atoms.len());
    for atom in &p1.atoms {
        let cart = structure.lattice().to_cart(atom.frac);
        out.push_str(&format!(
            "{} {:.10} {:.10} {:.10}",
            atom.element, cart[0], cart[1], cart[2]
        ));
        if let Some(site) = structure.sites().get(atom.site_idx) {
            for column in &columns {
                out.push(' ');
                out.push_str(&column_value(site.props().get(&column.key), column));
            }
        }
        out.push('\n');
    }
    Ok(out)
}

struct Column {
    key: String,
    kind: char,
    count: usize,
}

fn property_columns(structure: &Structure) -> Result<Vec<Column>, Error> {
    let mut keys: BTreeSet<String> = BTreeSet::new();
    for site in structure.sites() {
        keys.extend(site.props().keys().cloned());
    }
    let mut columns = Vec::new();
    for key in keys {
        let mut shape: Option<(char, usize)> = None;
        let mut complete = true;
        for site in structure.sites() {
            let Some(value) = site.props().get(&key) else {
                complete = false;
                break;
            };
            let Some(current) = column_shape(value) else {
                complete = false;
                break;
            };
            match shape {
                None => shape = Some(current),
                Some(expected) if expected == current => {}
                Some(_) => {
                    return Err(Error::InvalidInput(format!(
                        "property {key:?} has incompatible types across sites"
                    )));
                }
            }
        }
        if !complete {
            continue;
        }
        if let Some((kind, count)) = shape {
            columns.push(Column { key, kind, count });
        }
    }
    Ok(columns)
}

fn column_shape(value: &Value) -> Option<(char, usize)> {
    match value {
        Value::Number(number) if number.is_i64() || number.is_u64() => Some(('I', 1)),
        Value::Number(_) => Some(('R', 1)),
        Value::Bool(_) => Some(('L', 1)),
        Value::String(_) => Some(('S', 1)),
        Value::Array(values) => {
            let first = values.first()?;
            let (kind, _) = column_shape_scalar(first)?;
            for value in values {
                if column_shape_scalar(value)?.0 != kind {
                    return None;
                }
            }
            Some((kind, values.len()))
        }
        Value::Object(_) | Value::Null => None,
    }
}

fn column_shape_scalar(value: &Value) -> Option<(char, usize)> {
    match value {
        Value::Number(number) if number.is_i64() || number.is_u64() => Some(('I', 1)),
        Value::Number(_) => Some(('R', 1)),
        Value::Bool(_) => Some(('L', 1)),
        Value::String(_) => Some(('S', 1)),
        _ => None,
    }
}

fn column_value(value: Option<&Value>, column: &Column) -> String {
    match value {
        Some(Value::Array(values)) => values
            .iter()
            .map(scalar_token)
            .collect::<Vec<_>>()
            .join(" "),
        Some(value) => scalar_token(value),
        None => {
            let placeholder = match column.kind {
                'R' => "0.0",
                'I' => "0",
                'L' => "F",
                _ => "?",
            };
            vec![placeholder; column.count].join(" ")
        }
    }
}

fn scalar_token(value: &Value) -> String {
    match value {
        Value::String(text) => format!("\"{}\"", text.replace('"', "'")),
        other => other.to_string(),
    }
}

fn is_scalar(value: &Value) -> bool {
    matches!(value, Value::String(_) | Value::Number(_) | Value::Bool(_))
}

fn scalar_text(value: &Value) -> Option<String> {
    match value {
        Value::String(text) => Some(text.clone()),
        Value::Number(number) => Some(number.to_string()),
        Value::Bool(flag) => Some(if *flag { "T" } else { "F" }.to_owned()),
        _ => None,
    }
}

fn parse_lattice(comment: &str) -> Result<Option<Lattice>, Error> {
    let Some(start) = comment.find("Lattice=\"") else {
        return Ok(None);
    };
    let rest = &comment[start + 9..];
    let Some(end) = rest.find('"') else {
        return Err(Error::InvalidInput(
            "XYZ Lattice attribute is not closed".to_owned(),
        ));
    };
    let values: Vec<f64> = rest[..end]
        .split_whitespace()
        .filter_map(|value| value.parse().ok())
        .collect();
    if values.len() != 9 {
        return Err(Error::InvalidInput(
            "XYZ Lattice needs nine numbers".to_owned(),
        ));
    }
    Ok(Some(Lattice::new([
        [values[0], values[1], values[2]],
        [values[3], values[4], values[5]],
        [values[6], values[7], values[8]],
    ])?))
}
