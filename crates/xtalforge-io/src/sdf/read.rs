//! SDF and MDL MOL reader (Task 13.2.5).
//!
//! V2000 is read in full. V3000 is read as a documented subset: the `COUNTS`,
//! `ATOM` and `BOND` blocks, atom-map numbers and `CHG=` attributes. Files
//! without a cell become a molecule in a padded P1 box. Data fields after
//! `M  END` are kept in the `sdf_data` metadata object.

#![allow(clippy::indexing_slicing)] // fixed-format parser, lengths checked

use std::collections::BTreeMap;

use serde_json::{Map, Value};
use xtalforge_core::elements::parse_symbol;
use xtalforge_core::{Error, ExplicitBond, Structure};

use crate::molecule::{MoleculeAtom, molecule_box, parse_number};

struct Record {
    title: String,
    atoms: Vec<MoleculeAtom>,
    charges: Vec<f64>,
    bonds: Vec<(usize, usize, f64)>,
    data: BTreeMap<String, String>,
}

impl Record {
    fn empty() -> Self {
        Self {
            title: String::from("sdf"),
            atoms: Vec::new(),
            charges: Vec::new(),
            bonds: Vec::new(),
            data: BTreeMap::new(),
        }
    }
}

/// Reads the first SDF or MOL record.
///
/// # Errors
///
/// Returns [`Error::InvalidInput`] for a truncated counts line, a malformed
/// atom or bond line, or an unsupported V3000 field.
pub fn read_sdf(text: &str) -> Result<Structure, Error> {
    let records = parse_records(text)?;
    let record = records
        .into_iter()
        .next()
        .ok_or_else(|| Error::InvalidInput("SDF has no record".to_owned()))?;
    to_structure(record)
}

/// Reads every SDF record.
///
/// # Errors
///
/// Returns the error of the first record that cannot be parsed.
pub fn read_sdf_all(text: &str) -> Result<Vec<Structure>, Error> {
    parse_records(text)?.into_iter().map(to_structure).collect()
}

fn parse_records(text: &str) -> Result<Vec<Record>, Error> {
    let mut records = Vec::new();
    let mut chunk: Vec<&str> = Vec::new();
    for line in text.lines() {
        if line.trim() == "$$$$" {
            records.push(parse_record(&chunk.join("\n"))?);
            chunk.clear();
        } else {
            chunk.push(line);
        }
    }
    if chunk.iter().any(|line| !line.trim().is_empty()) {
        records.push(parse_record(&chunk.join("\n"))?);
    }
    Ok(records)
}

fn parse_record(text: &str) -> Result<Record, Error> {
    let lines: Vec<&str> = text.lines().collect();
    if lines.len() < 4 {
        return Err(Error::InvalidInput(
            "a MOL record has three header lines and a counts line".to_owned(),
        ));
    }
    let mut record = Record::empty();
    record.title = lines[0].trim().to_owned();
    let counts = lines[3];
    let v3000 = counts.contains("V3000");
    let (n_atoms, n_bonds) = if v3000 {
        find_v3000_counts(&lines)?
    } else {
        let atoms = counts.get(0..3).unwrap_or("").trim();
        let bonds = counts.get(3..6).unwrap_or("").trim();
        let atoms = atoms.parse::<usize>().map_err(|_| {
            Error::InvalidInput(format!("MOL counts line has no atom count: {counts:?}"))
        })?;
        let bonds = bonds.parse::<usize>().map_err(|_| {
            Error::InvalidInput(format!("MOL counts line has no bond count: {counts:?}"))
        })?;
        (atoms, bonds)
    };
    let mut index = 4;
    if v3000 {
        parse_v3000_body(&lines, &mut record)?;
    } else {
        for atom_index in 0..n_atoms {
            let line = lines
                .get(index)
                .ok_or_else(|| Error::InvalidInput("MOL ended before all atom lines".to_owned()))?;
            index += 1;
            record.atoms.push(parse_v2000_atom(line, atom_index)?);
            record.charges.push(parse_v2000_charge(line, atom_index)?);
        }
        for _ in 0..n_bonds {
            let line = lines
                .get(index)
                .ok_or_else(|| Error::InvalidInput("MOL ended before all bond lines".to_owned()))?;
            index += 1;
            if let Some(bond) = parse_v2000_bond(line)? {
                record.bonds.push(bond);
            }
        }
        while let Some(line) = lines.get(index) {
            index += 1;
            let trimmed = line.trim();
            if trimmed == "M  END" {
                break;
            }
            if let Some(entries) = trimmed.strip_prefix("M  CHG") {
                apply_charges(entries, &mut record.charges)?;
            }
        }
    }
    record.data = parse_data_fields(&lines[index.min(lines.len())..]);
    Ok(record)
}

fn parse_v2000_atom(line: &str, index: usize) -> Result<MoleculeAtom, Error> {
    let symbol = line
        .get(31..34)
        .map(str::trim)
        .filter(|value| !value.is_empty())
        .or_else(|| line.split_whitespace().nth(3))
        .ok_or_else(|| Error::InvalidInput(format!("MOL atom {} has no element", index + 1)))?;
    let element = parse_symbol(symbol)?;
    let cart = [
        parse_number(line.get(0..10).unwrap_or(""), "MOL x")?,
        parse_number(line.get(10..20).unwrap_or(""), "MOL y")?,
        parse_number(line.get(20..30).unwrap_or(""), "MOL z")?,
    ];
    Ok(MoleculeAtom {
        label: format!("{element}{}", index + 1),
        element,
        cart,
    })
}

fn parse_v2000_charge(line: &str, index: usize) -> Result<f64, Error> {
    let code = line
        .get(36..39)
        .map(str::trim)
        .filter(|value| !value.is_empty())
        .and_then(|value| value.parse::<i32>().ok())
        .unwrap_or(0);
    let charge = match code {
        1 => 3.0,
        2 => 2.0,
        3 => 1.0,
        5 => -1.0,
        6 => -2.0,
        7 => -3.0,
        0 | 4 => 0.0,
        other => {
            return Err(Error::InvalidInput(format!(
                "MOL atom {} has an unknown charge code {other}",
                index + 1
            )));
        }
    };
    Ok(charge)
}

fn parse_v2000_bond(line: &str) -> Result<Option<(usize, usize, f64)>, Error> {
    let fields: Vec<&str> = line.split_whitespace().collect();
    if fields.len() < 3 {
        return Ok(None);
    }
    let first = fields[0]
        .parse::<usize>()
        .map_err(|_| Error::InvalidInput(format!("bad MOL bond line: {line:?}")))?;
    let second = fields[1]
        .parse::<usize>()
        .map_err(|_| Error::InvalidInput(format!("bad MOL bond line: {line:?}")))?;
    let kind = fields[2]
        .parse::<u32>()
        .map_err(|_| Error::InvalidInput(format!("bad MOL bond order: {line:?}")))?;
    if first == 0 || second == 0 {
        return Err(Error::InvalidInput(
            "MOL bond atom indices are one-based and start at 1".to_owned(),
        ));
    }
    Ok(Some((first - 1, second - 1, bond_order(kind))))
}

fn bond_order(kind: u32) -> f64 {
    match kind {
        2 => 2.0,
        3 => 3.0,
        4 => 1.5,
        _ => 1.0,
    }
}

fn apply_charges(entries: &str, charges: &mut [f64]) -> Result<(), Error> {
    let fields: Vec<&str> = entries.split_whitespace().collect();
    let count = fields
        .first()
        .and_then(|value| value.parse::<usize>().ok())
        .unwrap_or(fields.len() / 2);
    let entries = fields.get(1..).unwrap_or(&[]);
    for pair in entries.chunks(2).take(count) {
        let (Some(atom), Some(charge)) = (pair.first(), pair.get(1)) else {
            break;
        };
        let atom = atom
            .parse::<usize>()
            .map_err(|_| Error::InvalidInput(format!("bad M  CHG atom: {atom:?}")))?;
        let charge = charge
            .parse::<f64>()
            .map_err(|_| Error::InvalidInput(format!("bad M  CHG charge: {charge:?}")))?;
        if let Some(slot) = atom.checked_sub(1).and_then(|at| charges.get_mut(at)) {
            *slot = charge;
        }
    }
    Ok(())
}

fn find_v3000_counts(lines: &[&str]) -> Result<(usize, usize), Error> {
    for line in lines {
        let trimmed = line.trim();
        if let Some(rest) = trimmed.strip_prefix("M  V30 COUNTS") {
            let fields: Vec<&str> = rest.split_whitespace().collect();
            let atoms = fields.first().and_then(|value| value.parse().ok());
            let bonds = fields.get(1).and_then(|value| value.parse().ok());
            if let (Some(atoms), Some(bonds)) = (atoms, bonds) {
                return Ok((atoms, bonds));
            }
        }
    }
    Err(Error::InvalidInput(
        "a V3000 record needs a COUNTS line".to_owned(),
    ))
}

fn parse_v3000_body(lines: &[&str], record: &mut Record) -> Result<(), Error> {
    let mut section = "";
    for line in lines {
        let trimmed = line.trim();
        let Some(body) = trimmed.strip_prefix("M  V30 ") else {
            continue;
        };
        let body = body.trim();
        if let Some(name) = body.strip_prefix("BEGIN ") {
            section = name.trim();
            continue;
        }
        if let Some(name) = body.strip_prefix("END ") {
            if name.trim() == section {
                section = "";
            }
            continue;
        }
        match section {
            "ATOM" => {
                let fields: Vec<&str> = body.split_whitespace().collect();
                if fields.len() < 5 {
                    return Err(Error::InvalidInput(format!(
                        "short V3000 atom line: {body:?}"
                    )));
                }
                let cart = [
                    parse_number(fields[2], "V3000 x")?,
                    parse_number(fields[3], "V3000 y")?,
                    parse_number(fields[4], "V3000 z")?,
                ];
                let element = parse_symbol(fields[1])?;
                let count = record.atoms.len();
                record.atoms.push(MoleculeAtom {
                    label: format!("{element}{}", count + 1),
                    element,
                    cart,
                });
                record
                    .charges
                    .push(v3000_attribute(&fields[6..], "CHG").unwrap_or(0.0));
            }
            "BOND" => {
                let fields: Vec<&str> = body.split_whitespace().collect();
                if fields.len() < 4 {
                    return Err(Error::InvalidInput(format!(
                        "short V3000 bond line: {body:?}"
                    )));
                }
                let kind = fields[1].parse::<u32>().unwrap_or(1);
                let first = fields[2]
                    .parse::<usize>()
                    .map_err(|_| Error::InvalidInput("bad V3000 bond atom".to_owned()))?;
                let second = fields[3]
                    .parse::<usize>()
                    .map_err(|_| Error::InvalidInput("bad V3000 bond atom".to_owned()))?;
                if first == 0 || second == 0 {
                    return Err(Error::InvalidInput(
                        "V3000 bond atom indices are one-based".to_owned(),
                    ));
                }
                record.bonds.push((first - 1, second - 1, bond_order(kind)));
            }
            _ => {}
        }
    }
    Ok(())
}

fn v3000_attribute(fields: &[&str], key: &str) -> Option<f64> {
    let prefix = format!("{key}=");
    fields
        .iter()
        .find_map(|field| field.strip_prefix(prefix.as_str()))
        .and_then(|value| value.parse().ok())
}

fn parse_data_fields(lines: &[&str]) -> BTreeMap<String, String> {
    let mut data = BTreeMap::new();
    let mut index = 0;
    while index < lines.len() {
        let line = lines[index].trim();
        index += 1;
        let Some(rest) = line.strip_prefix('>') else {
            continue;
        };
        let name = rest
            .trim()
            .trim_start_matches('<')
            .trim_end_matches('>')
            .to_owned();
        let mut value = String::new();
        while let Some(line) = lines.get(index) {
            index += 1;
            if line.trim().is_empty() {
                break;
            }
            if !value.is_empty() {
                value.push('\n');
            }
            value.push_str(line);
        }
        if !name.is_empty() {
            data.insert(name, value);
        }
    }
    data
}

fn to_structure(record: Record) -> Result<Structure, Error> {
    let bonds: Vec<ExplicitBond> = record
        .bonds
        .iter()
        .map(|&(i, j, order)| ExplicitBond {
            i,
            j,
            image: [0, 0, 0],
            order,
            remove: false,
        })
        .collect();
    let mut structure = molecule_box(&record.title, &record.atoms, &bonds)?;
    let mut builder = xtalforge_core::StructureBuilder::new(structure.name(), *structure.lattice());
    for site in structure.sites() {
        let mut site = site.clone();
        if let Some(charge) = record.charges.get(site_index(&site, &record))
            && *charge != 0.0
        {
            site = site.with_charge(*charge)?;
        }
        builder = builder.site(site);
    }
    for bond in &bonds {
        builder = builder.explicit_bond(*bond);
    }
    if !record.data.is_empty() {
        let map: Map<String, Value> = record
            .data
            .iter()
            .map(|(key, value)| (key.clone(), Value::String(value.clone())))
            .collect();
        builder = builder.metadata("sdf_data", Value::Object(map));
    }
    structure = builder.build()?;
    Ok(structure)
}

fn site_index(site: &xtalforge_core::Site, record: &Record) -> usize {
    record
        .atoms
        .iter()
        .position(|atom| atom.label == site.label())
        .unwrap_or(0)
}
