//! Tripos MOL2 reader and writer (Task 13.2.5).
//!
//! Reads the `MOLECULE`, `ATOM` and `BOND` sections of one molecule. A MOL2
//! file carries no cell, so the reader places the atoms in a padded P1 box.
//! Bond orders `ar`, `am`, `du` and `un` are accepted; the writer emits
//! `1`/`2`/`3`/`ar` from the stored orders and element symbols as atom types.

#![allow(clippy::indexing_slicing)] // fixed-format parser, lengths checked

use std::collections::BTreeSet;

use serde_json::Value;
use xtalforge_core::elements::parse_symbol;
use xtalforge_core::{Error, ExplicitBond, Structure};

use crate::molecule::{MoleculeAtom, molecule_box, parse_index, parse_number, unique_label};

/// Reads the first molecule of a Tripos MOL2 file.
///
/// # Errors
///
/// Returns [`Error::InvalidInput`] for a missing section, a bad coordinate, an
/// unknown atom type or a bond that names a missing atom.
pub fn read_mol2(text: &str) -> Result<Structure, Error> {
    let mut name = String::from("mol2");
    let mut atoms: Vec<MoleculeAtom> = Vec::new();
    let mut charges: Vec<Option<f64>> = Vec::new();
    let mut bonds: Vec<(usize, usize, f64)> = Vec::new();
    let mut charge_type: Option<String> = None;
    let mut section = "";
    let mut used = BTreeSet::new();
    let mut molecule_line = 0usize;
    for line in text.lines() {
        let trimmed = line.trim();
        if let Some(header) = trimmed
            .strip_prefix("@<TRIPOS>")
            .or_else(|| trimmed.strip_prefix("@<tripos>"))
        {
            section = match header.trim().to_ascii_uppercase().as_str() {
                "MOLECULE" => "molecule",
                "ATOM" => "atom",
                "BOND" => "bond",
                _ => "other",
            };
            molecule_line = 0;
            continue;
        }
        if trimmed.is_empty() || trimmed.starts_with('#') {
            continue;
        }
        match section {
            "molecule" => {
                if molecule_line == 0 && name == "mol2" {
                    name = trimmed.to_owned();
                } else if molecule_line == 3 {
                    charge_type = Some(trimmed.to_owned());
                }
                molecule_line += 1;
            }
            "atom" => {
                let fields: Vec<&str> = trimmed.split_whitespace().collect();
                if fields.len() < 6 {
                    return Err(Error::InvalidInput(
                        "a MOL2 atom line needs id, name, x, y, z and type".to_owned(),
                    ));
                }
                let cart = [
                    parse_number(fields[2], "MOL2 x")?,
                    parse_number(fields[3], "MOL2 y")?,
                    parse_number(fields[4], "MOL2 z")?,
                ];
                let element = element_from_type(fields[5])?;
                let label = unique_label(fields[1], &mut used);
                atoms.push(MoleculeAtom {
                    label,
                    element,
                    cart,
                });
                charges.push(fields.get(8).and_then(|value| value.parse().ok()));
            }
            "bond" => {
                let fields: Vec<&str> = trimmed.split_whitespace().collect();
                if fields.len() < 4 {
                    return Err(Error::InvalidInput(
                        "a MOL2 bond line needs id and two atoms".to_owned(),
                    ));
                }
                let first = parse_index(fields[1], "MOL2 bond")?;
                let second = parse_index(fields[2], "MOL2 bond")?;
                if first == 0 || second == 0 || first > atoms.len() || second > atoms.len() {
                    return Err(Error::InvalidInput(format!(
                        "MOL2 bond names atom {first} or {second}, but there are {} atoms",
                        atoms.len()
                    )));
                }
                bonds.push((first - 1, second - 1, bond_order(fields[3])));
            }
            _ => {}
        }
    }
    if atoms.is_empty() {
        return Err(Error::InvalidInput(
            "MOL2 has no @<TRIPOS>ATOM section".to_owned(),
        ));
    }
    let explicit: Vec<ExplicitBond> = bonds
        .iter()
        .map(|&(i, j, order)| ExplicitBond {
            i,
            j,
            image: [0, 0, 0],
            order,
            remove: false,
        })
        .collect();
    let mut structure = molecule_box(&name, &atoms, &explicit)?;
    if let Some(charge_type) = charge_type {
        structure = rebuild_with_metadata(structure, "mol2_charge_type", charge_type)?;
    }
    let mut sites = Vec::with_capacity(structure.sites().len());
    for (index, site) in structure.sites().iter().enumerate() {
        let mut site = site.clone();
        if let Some(Some(charge)) = charges.get(index) {
            site = site.with_charge(*charge)?;
        }
        sites.push(site);
    }
    rebuild_with_sites(&structure, sites)
}

fn element_from_type(atom_type: &str) -> Result<String, Error> {
    let base = atom_type
        .split_once('.')
        .map_or(atom_type, |(base, _)| base)
        .trim();
    let base = match base {
        "Du" | "LP" | "Lp" | "Any" | "ANY" => "X",
        _ => base,
    };
    let symbol = parse_symbol(base)?;
    Ok(symbol)
}

fn bond_order(kind: &str) -> f64 {
    match kind.to_ascii_lowercase().as_str() {
        "2" => 2.0,
        "3" => 3.0,
        "ar" => 1.5,
        _ => 1.0,
    }
}

fn rebuild_with_metadata(
    structure: Structure,
    key: &str,
    value: String,
) -> Result<Structure, Error> {
    let mut builder = xtalforge_core::StructureBuilder::new(structure.name(), *structure.lattice())
        .hall(structure.hall().to_owned())
        .metadata(key, Value::String(value));
    for operation in structure.operations() {
        builder = builder.operation(*operation);
    }
    for site in structure.sites() {
        builder = builder.site(site.clone());
    }
    for bond in structure.explicit_bonds() {
        builder = builder.explicit_bond(*bond);
    }
    for (key, value) in structure.metadata() {
        builder = builder.metadata(key.clone(), value.clone());
    }
    builder.build()
}

fn rebuild_with_sites(
    structure: &Structure,
    sites: Vec<xtalforge_core::Site>,
) -> Result<Structure, Error> {
    let mut builder = xtalforge_core::StructureBuilder::new(structure.name(), *structure.lattice())
        .hall(structure.hall().to_owned());
    for operation in structure.operations() {
        builder = builder.operation(*operation);
    }
    for site in sites {
        builder = builder.site(site);
    }
    for bond in structure.explicit_bonds() {
        builder = builder.explicit_bond(*bond);
    }
    for (key, value) in structure.metadata() {
        builder = builder.metadata(key.clone(), value.clone());
    }
    builder.build()
}

/// Writes a structure as a single Tripos MOL2 molecule.
///
/// # Errors
///
/// Returns an error for an empty structure or a failed P1 expansion.
pub fn write_mol2(structure: &Structure) -> Result<String, Error> {
    let p1 = structure.p1()?;
    if p1.atoms.is_empty() {
        return Err(Error::InvalidInput(
            "a MOL2 file cannot hold a structure with no atoms".to_owned(),
        ));
    }
    let mut out = String::from("@<TRIPOS>MOLECULE\n");
    out.push_str(&format!("{}\n", structure.name().replace('\n', " ")));
    out.push_str(&format!(
        "{} {}\nSMALL\nUSER_CHARGES\n\n",
        p1.atoms.len(),
        structure
            .explicit_bonds()
            .iter()
            .filter(|bond| !bond.remove)
            .count()
    ));
    out.push_str("@<TRIPOS>ATOM\n");
    let mut used = BTreeSet::new();
    for (index, atom) in p1.atoms.iter().enumerate() {
        let site = structure.sites().get(atom.site_idx);
        let charge = site.and_then(|site| site.charge()).unwrap_or(0.0);
        let base = if atom.label.is_empty() {
            format!("{}{}", atom.element, index + 1)
        } else {
            atom.label.clone()
        };
        let label = unique_label(&base, &mut used);
        let [x, y, z] = structure.lattice().to_cart(atom.frac);
        out.push_str(&format!(
            "{:>7} {:<8} {:>10.4} {:>10.4} {:>10.4} {:<5} 1 MOL  {:>8.4}\n",
            index + 1,
            label,
            x,
            y,
            z,
            atom.element,
            charge
        ));
    }
    let mut site_atom: Vec<Option<usize>> = vec![None; structure.sites().len()];
    for (index, atom) in p1.atoms.iter().enumerate() {
        if let Some(slot) = site_atom.get_mut(atom.site_idx)
            && slot.is_none()
        {
            *slot = Some(index);
        }
    }
    let bonds: Vec<(usize, usize, f64)> = structure
        .explicit_bonds()
        .iter()
        .filter(|bond| !bond.remove)
        .filter_map(|bond| {
            let i = site_atom.get(bond.i).copied().flatten()?;
            let j = site_atom.get(bond.j).copied().flatten()?;
            (i != j).then_some((i, j, bond.order))
        })
        .collect();
    if !bonds.is_empty() {
        out.push_str("@<TRIPOS>BOND\n");
        for (index, (i, j, order)) in bonds.iter().enumerate() {
            out.push_str(&format!(
                "{:>6} {:>4} {:>4} {}\n",
                index + 1,
                i + 1,
                j + 1,
                bond_kind(*order)
            ));
        }
    }
    Ok(out)
}

fn bond_kind(order: f64) -> &'static str {
    if (order - 1.5).abs() < 0.2 {
        "ar"
    } else if (order - 3.0).abs() < 0.2 {
        "3"
    } else if (order - 2.0).abs() < 0.2 {
        "2"
    } else {
        "1"
    }
}
