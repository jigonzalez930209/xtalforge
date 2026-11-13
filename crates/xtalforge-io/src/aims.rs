//! FHI-aims `geometry.in` reader and writer (Task 13.2.2).
//!
//! Reads `lattice_vector`, `atom`, `atom_frac`, `velocity`,
//! `constrain_relaxation`, `initial_charge` and `initial_moment`. The writer
//! emits the periodic form; per-atom constraint, velocity and initial-charge
//! properties are written back when they are present.

#![allow(clippy::indexing_slicing)] // fixed-format parser, lengths checked

use serde_json::Value;
use xtalforge_core::elements::parse_symbol;
use xtalforge_core::{Error, Lattice, Site, Structure, StructureBuilder};

use crate::molecule::{MoleculeAtom, molecule_box, parse_number};

/// Reads an FHI-aims `geometry.in` file.
///
/// # Errors
///
/// Returns [`Error::InvalidInput`] for a malformed line, an incomplete lattice
/// or fractional coordinates without a lattice.
pub fn read_aims(text: &str) -> Result<Structure, Error> {
    let mut lattice_rows: Vec<[f64; 3]> = Vec::new();
    let mut records: Vec<AtomRecord> = Vec::new();
    for line in text.lines() {
        let line = line.split('#').next().unwrap_or("").trim();
        if line.is_empty() {
            continue;
        }
        let fields: Vec<&str> = line.split_whitespace().collect();
        match fields.first().copied().unwrap_or("") {
            "lattice_vector" => {
                if fields.len() < 4 {
                    return Err(Error::InvalidInput(
                        "FHI-aims lattice_vector needs three values".to_owned(),
                    ));
                }
                lattice_rows.push([
                    parse_number(fields[1], "FHI-aims lattice x")?,
                    parse_number(fields[2], "FHI-aims lattice y")?,
                    parse_number(fields[3], "FHI-aims lattice z")?,
                ]);
            }
            "atom" | "atom_frac" => {
                if fields.len() < 5 {
                    return Err(Error::InvalidInput(
                        "FHI-aims atom needs three coordinates and an element".to_owned(),
                    ));
                }
                let element = parse_symbol(fields[4])?;
                let cart = [
                    parse_number(fields[1], "FHI-aims x")?,
                    parse_number(fields[2], "FHI-aims y")?,
                    parse_number(fields[3], "FHI-aims z")?,
                ];
                records.push(AtomRecord {
                    element,
                    position: cart,
                    fractional: fields[0] == "atom_frac",
                    props: Vec::new(),
                });
            }
            "velocity" => {
                let values = numeric_row(&fields, "FHI-aims velocity")?;
                if let Some(last) = records.last_mut() {
                    last.props
                        .push(("velocity".to_owned(), Value::from(values.to_vec())));
                }
            }
            "constrain_relaxation" => {
                if let Some(last) = records.last_mut() {
                    let flag = fields
                        .get(1)
                        .is_some_and(|value| !value.starts_with(".false."));
                    last.props
                        .push(("constrain_relaxation".to_owned(), Value::Bool(flag)));
                }
            }
            "initial_charge" | "initial_moment" => {
                let value = fields
                    .get(1)
                    .and_then(|value| value.parse::<f64>().ok())
                    .ok_or_else(|| {
                        Error::InvalidInput(format!("FHI-aims {} needs a number", fields[0]))
                    })?;
                if let Some(last) = records.last_mut() {
                    last.props.push((fields[0].to_owned(), Value::from(value)));
                }
            }
            "empty" | "homogeneous_field" | "symmetry_region" => {}
            _ => {}
        }
    }
    if records.is_empty() {
        return Err(Error::InvalidInput(
            "FHI-aims geometry.in has no atoms".to_owned(),
        ));
    }
    if lattice_rows.is_empty() {
        if records.iter().any(|record| record.fractional) {
            return Err(Error::InvalidInput(
                "FHI-aims atom_frac needs three lattice_vector rows".to_owned(),
            ));
        }
        let atoms: Vec<MoleculeAtom> = records
            .iter()
            .enumerate()
            .map(|(index, record)| MoleculeAtom {
                label: format!("{}{}", record.element, index + 1),
                element: record.element.clone(),
                cart: record.position,
            })
            .collect();
        let base = molecule_box("fhi-aims", &atoms, &[])?;
        let mut builder = StructureBuilder::new("fhi-aims", *base.lattice());
        for (site, record) in base.sites().iter().zip(records.iter()) {
            let mut site = site.clone();
            for (key, value) in &record.props {
                site = site.with_prop(key.clone(), value.clone());
            }
            builder = builder.site(site);
        }
        return builder.build();
    }
    if lattice_rows.len() != 3 {
        return Err(Error::InvalidInput(format!(
            "FHI-aims needs exactly three lattice_vector rows, got {}",
            lattice_rows.len()
        )));
    }
    let lattice = Lattice::new([lattice_rows[0], lattice_rows[1], lattice_rows[2]])?;
    let mut builder = StructureBuilder::new("fhi-aims", lattice);
    for (index, record) in records.iter().enumerate() {
        let frac = if record.fractional {
            record.position
        } else {
            lattice.to_frac(record.position)
        };
        let mut site = Site::new(
            format!("{}{}", record.element, index + 1),
            &record.element,
            frac,
        )?;
        for (key, value) in &record.props {
            site = site.with_prop(key.clone(), value.clone());
        }
        builder = builder.site(site);
    }
    builder.build()
}

/// Writes a P1 structure as an FHI-aims `geometry.in`.
///
/// # Errors
///
/// Returns an error for an empty structure or a failed P1 expansion.
pub fn write_aims(structure: &Structure) -> Result<String, Error> {
    let p1 = structure.p1()?;
    if p1.atoms.is_empty() {
        return Err(Error::InvalidInput(
            "an FHI-aims geometry cannot hold a structure with no atoms".to_owned(),
        ));
    }
    let mut out = String::from("# geometry.in written by xtalforge\n");
    for row in structure.lattice().matrix() {
        out.push_str(&format!(
            "lattice_vector {:.10} {:.10} {:.10}\n",
            row[0], row[1], row[2]
        ));
    }
    for atom in &p1.atoms {
        let site = structure.sites().get(atom.site_idx);
        out.push_str(&format!(
            "atom_frac {:.10} {:.10} {:.10} {}\n",
            atom.frac[0], atom.frac[1], atom.frac[2], atom.element
        ));
        let Some(site) = site else {
            continue;
        };
        if let Some(value) = site.props().get("velocity").and_then(Value::as_array) {
            let values: Vec<String> = value
                .iter()
                .filter_map(Value::as_f64)
                .map(|value| format!("{value:.10}"))
                .collect();
            if values.len() == 3 {
                out.push_str(&format!("  velocity {}\n", values.join(" ")));
            }
        }
        if let Some(charge) = site.props().get("initial_charge").and_then(Value::as_f64) {
            out.push_str(&format!("  initial_charge {charge:.6}\n"));
        }
        if let Some(moment) = site.props().get("initial_moment").and_then(Value::as_f64) {
            out.push_str(&format!("  initial_moment {moment:.6}\n"));
        }
        if site
            .props()
            .get("constrain_relaxation")
            .and_then(Value::as_bool)
            .unwrap_or(false)
        {
            out.push_str("  constrain_relaxation .true.\n");
        }
    }
    Ok(out)
}

struct AtomRecord {
    element: String,
    position: [f64; 3],
    fractional: bool,
    props: Vec<(String, Value)>,
}

fn numeric_row(fields: &[&str], what: &str) -> Result<[f64; 3], Error> {
    if fields.len() < 4 {
        return Err(Error::InvalidInput(format!("{what} needs three values")));
    }
    Ok([
        parse_number(fields[1], what)?,
        parse_number(fields[2], what)?,
        parse_number(fields[3], what)?,
    ])
}
