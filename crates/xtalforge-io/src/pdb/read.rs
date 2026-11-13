//! PDB reader (Task 13.2.3).
//!
//! Reads the first model of an ATOM/HETATM file, the `CRYST1` cell and space
//! group, `CONECT` bonds, occupancy, B factors (as `U_iso = B / 8 pi^2`) and
//! the residue/chain hierarchy as site properties. A file without `CRYST1`
//! becomes a molecule in a padded P1 box, like the plain-XYZ path.

use std::collections::{BTreeMap, BTreeSet};

use serde_json::Value;
use xtalforge_core::elements::parse_symbol;
use xtalforge_core::{Error, ExplicitBond, Lattice, Site, Structure, StructureBuilder};

use crate::molecule::{MoleculeAtom, molecule_box, unique_label};

struct PdbAtom {
    serial: u32,
    label: String,
    element: String,
    cart: [f64; 3],
    occupancy: f64,
    b_factor: Option<f64>,
    residue: String,
    chain: String,
    residue_number: String,
    alt_loc: String,
}

/// Reads a PDB file.
///
/// # Errors
///
/// Returns [`Error::InvalidInput`] for a malformed atom record, an invalid
/// `CRYST1` cell or a structure with no atoms.
pub fn read_pdb(text: &str) -> Result<Structure, Error> {
    let mut cell: Option<[f64; 6]> = None;
    let mut space_group = String::new();
    let mut title = String::new();
    let mut atoms: Vec<PdbAtom> = Vec::new();
    let mut used = BTreeSet::new();
    let mut conect: Vec<(u32, u32)> = Vec::new();
    let mut saw_model = false;
    let mut first_model_done = false;
    let mut in_model = true;
    for line in text.lines() {
        let record = line.get(0..6).unwrap_or("").trim();
        match record {
            "HEADER" => {
                if title.is_empty() {
                    title = line.get(10..50).unwrap_or("").trim().to_owned();
                }
            }
            "TITLE" => {
                let part = line.get(10..).unwrap_or("").trim();
                if !part.is_empty() {
                    if !title.is_empty() {
                        title.push(' ');
                    }
                    title.push_str(part);
                }
            }
            "CRYST1" => {
                let (values, group) = parse_cryst1(line)?;
                cell = Some(values);
                space_group = group;
            }
            "MODEL" => {
                if saw_model {
                    first_model_done = true;
                }
                saw_model = true;
                in_model = !first_model_done;
            }
            "ENDMDL" => in_model = false,
            "ATOM" | "HETATM" => {
                if first_model_done || (saw_model && !in_model) {
                    continue;
                }
                atoms.push(parse_atom(line, &mut used)?);
            }
            "CONECT" => conect.extend(parse_conect(line)?),
            _ => {}
        }
    }
    if atoms.is_empty() {
        return Err(Error::InvalidInput(
            "PDB has no ATOM or HETATM records".to_owned(),
        ));
    }
    let serial_to_atom: BTreeMap<u32, usize> = atoms
        .iter()
        .enumerate()
        .map(|(index, atom)| (atom.serial, index))
        .collect();
    let mut bonds = Vec::new();
    for (first, second) in conect {
        let (Some(&i), Some(&j)) = (serial_to_atom.get(&first), serial_to_atom.get(&second)) else {
            continue;
        };
        if i == j {
            continue;
        }
        let bond = if i <= j {
            ExplicitBond {
                i,
                j,
                image: [0, 0, 0],
                order: 1.0,
                remove: false,
            }
        } else {
            ExplicitBond {
                i: j,
                j: i,
                image: [0, 0, 0],
                order: 1.0,
                remove: false,
            }
        };
        if !bonds.contains(&bond) {
            bonds.push(bond);
        }
    }
    let metadata = Metadata {
        title,
        space_group,
        residue_chains: atoms
            .iter()
            .map(|atom| (atom.chain.clone(), atom.residue.clone()))
            .collect(),
    };
    match cell {
        Some(parameters) => periodic(atoms, parameters, &metadata, &bonds),
        None => molecular(atoms, &metadata, &bonds),
    }
}

struct Metadata {
    title: String,
    space_group: String,
    residue_chains: BTreeSet<(String, String)>,
}

fn periodic(
    atoms: Vec<PdbAtom>,
    parameters: [f64; 6],
    metadata: &Metadata,
    bonds: &[ExplicitBond],
) -> Result<Structure, Error> {
    let lattice = Lattice::from_parameters(
        parameters[0],
        parameters[1],
        parameters[2],
        parameters[3],
        parameters[4],
        parameters[5],
    )?;
    let name = if metadata.title.is_empty() {
        "pdb"
    } else {
        metadata.title.as_str()
    };
    let (hall, operations, warning) = space_group(&metadata.space_group);
    let mut builder = StructureBuilder::new(name, lattice).hall(hall);
    for operation in operations {
        builder = builder.operation(operation);
    }
    for atom in &atoms {
        builder = builder.site(atom_site(atom, lattice.to_frac(atom.cart))?);
    }
    for bond in bonds {
        builder = builder.explicit_bond(*bond);
    }
    builder = builder.metadata(
        "residue_chains",
        Value::Array(
            metadata
                .residue_chains
                .iter()
                .map(|(chain, residue)| {
                    Value::Array(vec![
                        Value::String(chain.clone()),
                        Value::String(residue.clone()),
                    ])
                })
                .collect(),
        ),
    );
    if let Some(warning) = warning {
        builder = builder.metadata("pdb_symmetry_warning", Value::String(warning));
    }
    if !metadata.space_group.is_empty() {
        builder = builder.metadata(
            "pdb_space_group",
            Value::String(metadata.space_group.clone()),
        );
    }
    builder.build()
}

fn molecular(
    atoms: Vec<PdbAtom>,
    metadata: &Metadata,
    bonds: &[ExplicitBond],
) -> Result<Structure, Error> {
    let name = if metadata.title.is_empty() {
        "pdb"
    } else {
        metadata.title.as_str()
    };
    let points: Vec<MoleculeAtom> = atoms
        .iter()
        .map(|atom| MoleculeAtom {
            label: atom.label.clone(),
            element: atom.element.clone(),
            cart: atom.cart,
        })
        .collect();
    let structure = molecule_box(name, &points, bonds)?;
    let mut builder = StructureBuilder::new(structure.name(), *structure.lattice());
    for (site, atom) in structure.sites().iter().zip(atoms.iter()) {
        builder = builder.site(decorate(site.clone(), atom));
    }
    for bond in bonds {
        builder = builder.explicit_bond(*bond);
    }
    builder = builder.metadata(
        "residue_chains",
        Value::Array(
            metadata
                .residue_chains
                .iter()
                .map(|(chain, residue)| {
                    Value::Array(vec![
                        Value::String(chain.clone()),
                        Value::String(residue.clone()),
                    ])
                })
                .collect(),
        ),
    );
    builder.build()
}

fn space_group(text: &str) -> (String, Vec<xtalforge_core::SymOp>, Option<String>) {
    if text.trim().is_empty() {
        return ("P 1".to_owned(), Vec::new(), None);
    }
    match xtalforge_symmetry::SpaceGroup::parse(text) {
        Ok(group) => (group.hall_symbol, group.operations, None),
        Err(error) => (
            "P 1".to_owned(),
            Vec::new(),
            Some(format!(
                "CRYST1 space group {text:?} was not recognised ({error}); the file is read as P1"
            )),
        ),
    }
}

fn atom_site(atom: &PdbAtom, frac: [f64; 3]) -> Result<Site, Error> {
    Ok(decorate(Site::new(&atom.label, &atom.element, frac)?, atom))
}

fn decorate(mut site: Site, atom: &PdbAtom) -> Site {
    if atom.occupancy < 1.0
        && let Ok(updated) = site.clone().with_occupancy(atom.occupancy)
    {
        site = updated;
    }
    if let Some(b_factor) = atom.b_factor {
        let u_iso = b_factor / (8.0 * std::f64::consts::PI * std::f64::consts::PI);
        if u_iso >= 0.0
            && let Ok(updated) = site.clone().with_u_iso(u_iso)
        {
            site = updated;
        }
    }
    site = site
        .with_prop("pdb_atom_name", Value::String(atom.label.clone()))
        .with_prop("pdb_serial", Value::from(atom.serial))
        .with_prop("residue_name", Value::String(atom.residue.clone()))
        .with_prop("chain_id", Value::String(atom.chain.clone()))
        .with_prop("residue_number", Value::String(atom.residue_number.clone()));
    if !atom.alt_loc.is_empty() {
        site = site.with_prop("alt_loc", Value::String(atom.alt_loc.clone()));
    }
    site
}

fn parse_atom(line: &str, used: &mut BTreeSet<String>) -> Result<PdbAtom, Error> {
    let serial = field(line, 6, 11)
        .and_then(|text| text.parse().ok())
        .ok_or_else(|| Error::InvalidInput("PDB atom record has no serial".to_owned()))?;
    let name = field(line, 12, 16).unwrap_or("X").to_owned();
    let alt_loc = field(line, 16, 17).unwrap_or("").to_owned();
    let residue = field(line, 17, 20).unwrap_or("UNL").trim().to_owned();
    let chain = field(line, 21, 22).unwrap_or("").trim().to_owned();
    let residue_number = field(line, 22, 26).unwrap_or("1").trim().to_owned();
    let cart = [
        parse_field(line, 30, 38, "x")?,
        parse_field(line, 38, 46, "y")?,
        parse_field(line, 46, 54, "z")?,
    ];
    let occupancy = field(line, 54, 60)
        .and_then(|text| text.parse::<f64>().ok())
        .filter(|value| value.is_finite())
        .unwrap_or(1.0)
        .clamp(0.0, 1.0);
    let b_factor = field(line, 60, 66)
        .and_then(|text| text.parse::<f64>().ok())
        .filter(|value| value.is_finite() && *value >= 0.0);
    let declared = field(line, 76, 78).unwrap_or("");
    let element = element_from_name(&name, declared);
    let label = unique_label(name.trim(), used);
    Ok(PdbAtom {
        serial,
        label,
        element,
        cart,
        occupancy,
        b_factor,
        residue,
        chain,
        residue_number,
        alt_loc,
    })
}

fn element_from_name(name: &str, declared: &str) -> String {
    for candidate in [declared.trim(), name.trim()] {
        if candidate.is_empty() {
            continue;
        }
        let stripped = candidate.trim_start_matches(|character: char| character.is_ascii_digit());
        if let Ok(symbol) = parse_symbol(stripped) {
            return symbol;
        }
        let mut characters = stripped.chars();
        if let Some(first) = characters.next() {
            let two: String = std::iter::once(first).chain(characters.next()).collect();
            if let Ok(symbol) = parse_symbol(&two) {
                return symbol;
            }
            if let Ok(symbol) = parse_symbol(&first.to_string()) {
                return symbol;
            }
        }
    }
    "X".to_owned()
}

fn parse_cryst1(line: &str) -> Result<([f64; 6], String), Error> {
    let tokens: Vec<&str> = line
        .split_whitespace()
        .skip_while(|token| token.eq_ignore_ascii_case("CRYST1"))
        .collect();
    let mut values = [0.0; 6];
    let mut parsed = 0usize;
    for token in &tokens {
        if parsed == 6 {
            break;
        }
        if let Ok(value) = token.parse::<f64>() {
            if let Some(slot) = values.get_mut(parsed) {
                *slot = value;
            }
            parsed += 1;
        }
    }
    if parsed < 6 {
        for (slot, (start, end)) in [(6, 15), (15, 25), (25, 35), (35, 45), (45, 55), (55, 65)]
            .iter()
            .enumerate()
        {
            let text = field(line, *start, *end).unwrap_or("");
            if let Ok(value) = text.parse::<f64>() {
                if let Some(slot) = values.get_mut(slot) {
                    *slot = value;
                }
                parsed += 1;
            }
        }
    }
    if parsed < 6 {
        return Err(Error::InvalidInput(
            "a PDB CRYST1 record needs six cell parameters".to_owned(),
        ));
    }
    let fixed = line.get(65..70).unwrap_or("").trim();
    let group = if fixed.is_empty() || fixed.chars().any(|character| character.is_ascii_digit()) {
        let mut rest: Vec<&str> = tokens.iter().skip(6).copied().collect();
        if rest
            .last()
            .is_some_and(|token| token.parse::<u32>().is_ok())
        {
            rest.pop();
        }
        rest.join(" ")
    } else {
        fixed.to_owned()
    };
    Ok((values, group))
}

fn parse_conect(line: &str) -> Result<Vec<(u32, u32)>, Error> {
    let mut serials = Vec::new();
    for chunk in line.get(6..).unwrap_or("").as_bytes().chunks(5) {
        let text = std::str::from_utf8(chunk).unwrap_or("").trim();
        if text.is_empty() {
            continue;
        }
        if let Ok(serial) = text.parse::<u32>() {
            serials.push(serial);
        }
    }
    let Some(first) = serials.first().copied() else {
        return Err(Error::InvalidInput(
            "a PDB CONECT record needs an atom serial".to_owned(),
        ));
    };
    Ok(serials
        .iter()
        .skip(1)
        .map(|serial| (first, *serial))
        .collect())
}

fn field(line: &str, start: usize, end: usize) -> Option<&str> {
    line.get(start..end).map(str::trim)
}

fn parse_field(line: &str, start: usize, end: usize, what: &str) -> Result<f64, Error> {
    let text = field(line, start, end).unwrap_or("");
    text.parse()
        .map_err(|_| Error::InvalidInput(format!("PDB atom has a bad {what} coordinate")))
}
