//! CP2K input reader and writer (Task 13.2.2).
//!
//! The reader walks the `&SECTION` tree and uses `&CELL` (`ABC` plus
//! `ALPHA_BETA_GAMMA`, or the three vectors `A`, `B`, `C`) together with
//! `&COORD` (`SCALED` or Cartesian) and the `ELEMENT` entry of `&KIND`. The
//! writer emits a complete `&FORCE_EVAL` input with an angstrom cell.

#![allow(clippy::indexing_slicing)] // fixed-format parser, lengths checked

use std::collections::BTreeMap;

use serde_json::Value;
use xtalforge_core::{Error, Lattice, Site, Structure, StructureBuilder};

use crate::molecule::{MoleculeAtom, molecule_box, parse_number};

struct Section {
    name: String,
    args: Vec<String>,
    entries: Vec<(String, Vec<String>)>,
    children: Vec<Section>,
}

impl Section {
    fn find(&self, name: &str) -> Option<&Section> {
        if self.name.eq_ignore_ascii_case(name) {
            return Some(self);
        }
        self.children.iter().find_map(|child| child.find(name))
    }

    fn keyword(&self, name: &str) -> Option<&[String]> {
        self.entries
            .iter()
            .find(|(key, _)| key.eq_ignore_ascii_case(name))
            .map(|(_, values)| values.as_slice())
    }

    fn has(&self, name: &str) -> bool {
        self.entries
            .iter()
            .any(|(key, _)| key.eq_ignore_ascii_case(name))
            || self
                .children
                .iter()
                .any(|child| child.name.eq_ignore_ascii_case(name))
    }
}

/// Reads a CP2K input file.
///
/// # Errors
///
/// Returns [`Error::InvalidInput`] for an unbalanced section, a missing
/// `&COORD` or a malformed cell or coordinate.
pub fn read_cp2k(text: &str) -> Result<Structure, Error> {
    let cleaned: Vec<&str> = text
        .lines()
        .map(|line| line.split('#').next().unwrap_or("").trim())
        .filter(|line| !line.is_empty())
        .collect();
    let mut index = 0usize;
    let root = parse_section(&cleaned, &mut index, "")?;
    let coord = root
        .find("COORD")
        .ok_or_else(|| Error::InvalidInput("CP2K input has no &COORD section".to_owned()))?;
    let scaled = coord
        .args
        .iter()
        .any(|arg| arg.eq_ignore_ascii_case("SCALED"))
        || coord.has("SCALED");
    let kinds = kind_elements(&root);
    let mut atoms = Vec::new();
    for (keyword, values) in &coord.entries {
        if keyword.eq_ignore_ascii_case("SCALED")
            || keyword.eq_ignore_ascii_case("UNIT")
            || keyword.eq_ignore_ascii_case("FILE")
        {
            continue;
        }
        if values.len() < 3 {
            return Err(Error::InvalidInput(format!(
                "CP2K &COORD line {keyword:?} needs three coordinates"
            )));
        }
        let element = kinds
            .get(keyword)
            .cloned()
            .or_else(|| xtalforge_core::elements::parse_symbol(keyword).ok())
            .ok_or_else(|| {
                Error::InvalidInput(format!("CP2K &COORD uses unknown element {keyword:?}"))
            })?;
        let position = [
            parse_number(&values[0], "CP2K x")?,
            parse_number(&values[1], "CP2K y")?,
            parse_number(&values[2], "CP2K z")?,
        ];
        atoms.push((element, position));
    }
    if atoms.is_empty() {
        return Err(Error::InvalidInput(
            "CP2K &COORD section has no atoms".to_owned(),
        ));
    }
    let lattice = match cp2k_cell(&root)? {
        Some(lattice) => lattice,
        None => {
            if scaled {
                return Err(Error::InvalidInput(
                    "CP2K &COORD is SCALED but there is no &CELL".to_owned(),
                ));
            }
            let points: Vec<MoleculeAtom> = atoms
                .iter()
                .enumerate()
                .map(|(index, (element, position))| MoleculeAtom {
                    label: format!("{element}{}", index + 1),
                    element: element.clone(),
                    cart: *position,
                })
                .collect();
            let base = molecule_box("cp2k", &points, &[])?;
            *base.lattice()
        }
    };
    let name = root
        .find("GLOBAL")
        .and_then(|global| global.keyword("PROJECT"))
        .and_then(|values| values.first())
        .cloned()
        .unwrap_or_else(|| "cp2k".to_owned());
    let mut builder = StructureBuilder::new(&name, lattice);
    for (index, (element, position)) in atoms.iter().enumerate() {
        let frac = if scaled {
            *position
        } else {
            lattice.to_frac(*position)
        };
        builder = builder.site(Site::new(format!("{element}{}", index + 1), element, frac)?);
    }
    builder = builder.metadata("cp2k_project", Value::String(name));
    builder.build()
}

fn parse_section(lines: &[&str], index: &mut usize, expected: &str) -> Result<Section, Error> {
    let mut section = Section {
        name: expected.to_owned(),
        args: Vec::new(),
        entries: Vec::new(),
        children: Vec::new(),
    };
    while *index < lines.len() {
        let line = lines.get(*index).copied().unwrap_or("").trim();
        *index += 1;
        if let Some(rest) = line.strip_prefix('&') {
            let mut parts = rest.split_whitespace();
            let name = parts.next().unwrap_or("").to_owned();
            if name.eq_ignore_ascii_case("END") {
                if expected.is_empty() {
                    continue;
                }
                return Ok(section);
            }
            let child = parse_section(lines, index, &name)?;
            section.children.push(child);
            continue;
        }
        let mut parts = line.split_whitespace();
        let Some(keyword) = parts.next() else {
            continue;
        };
        let values: Vec<String> = parts.map(str::to_owned).collect();
        section.entries.push((keyword.to_owned(), values));
    }
    if expected.is_empty() {
        Ok(section)
    } else {
        Err(Error::InvalidInput(format!(
            "CP2K section &{expected} is not closed by &END"
        )))
    }
}

fn kind_elements(root: &Section) -> BTreeMap<String, String> {
    let mut out = BTreeMap::new();
    collect_kinds(root, &mut out);
    out
}

fn collect_kinds(section: &Section, out: &mut BTreeMap<String, String>) {
    for child in &section.children {
        if child.name.eq_ignore_ascii_case("KIND")
            && let Some(label) = child.args.first()
        {
            let element = child
                .keyword("ELEMENT")
                .and_then(|values| values.first())
                .cloned()
                .or_else(|| xtalforge_core::elements::parse_symbol(label).ok());
            if let Some(element) = element {
                out.insert(label.clone(), element.clone());
                out.entry(label.to_ascii_uppercase()).or_insert(element);
            }
        }
        collect_kinds(child, out);
    }
}

fn cp2k_cell(root: &Section) -> Result<Option<Lattice>, Error> {
    let Some(cell) = root.find("CELL") else {
        return Ok(None);
    };
    if let Some(rows) = cell_vectors(cell)? {
        return Lattice::new(rows).map(Some);
    }
    if let Some(abc) = cell.keyword("ABC") {
        let values = numbers(abc, "CP2K CELL ABC")?;
        if values.len() < 3 {
            return Err(Error::InvalidInput("CP2K CELL ABC needs a b c".to_owned()));
        }
        let angles = cell
            .keyword("ALPHA_BETA_GAMMA")
            .map(|values| numbers(values, "CP2K CELL angles"))
            .transpose()?
            .unwrap_or_else(|| vec![90.0, 90.0, 90.0]);
        if angles.len() < 3 {
            return Err(Error::InvalidInput(
                "CP2K ALPHA_BETA_GAMMA needs three angles".to_owned(),
            ));
        }
        return Lattice::from_parameters(
            values[0], values[1], values[2], angles[0], angles[1], angles[2],
        )
        .map(Some);
    }
    Err(Error::InvalidInput(
        "CP2K &CELL needs ABC or the A, B, C vectors".to_owned(),
    ))
}

fn cell_vectors(cell: &Section) -> Result<Option<[[f64; 3]; 3]>, Error> {
    let mut rows = [[0.0; 3]; 3];
    let mut found = false;
    for (slot, label) in ["A", "B", "C"].iter().enumerate() {
        let Some(values) = cell.keyword(label) else {
            continue;
        };
        let numbers = numbers(values, "CP2K CELL vector")?;
        if numbers.len() < 3 {
            return Err(Error::InvalidInput(format!(
                "CP2K CELL {label} needs three values"
            )));
        }
        if let Some(row) = rows.get_mut(slot) {
            *row = [numbers[0], numbers[1], numbers[2]];
        }
        found = true;
    }
    Ok(found.then_some(rows))
}

fn numbers(values: &[String], what: &str) -> Result<Vec<f64>, Error> {
    values
        .iter()
        .map(|value| parse_number(value, what))
        .collect()
}

/// Writes a P1 structure as a CP2K input file.
///
/// # Errors
///
/// Returns an error for an empty structure or a failed P1 expansion.
pub fn write_cp2k(structure: &Structure) -> Result<String, Error> {
    let p1 = structure.p1()?;
    if p1.atoms.is_empty() {
        return Err(Error::InvalidInput(
            "a CP2K input cannot hold a structure with no atoms".to_owned(),
        ));
    }
    let [a, b, c, alpha, beta, gamma] = structure.lattice().parameters();
    let project: String = structure
        .name()
        .chars()
        .map(|character| {
            if character.is_ascii_alphanumeric() {
                character
            } else {
                '_'
            }
        })
        .collect();
    let mut out = String::from("&GLOBAL\n");
    out.push_str(&format!(
        "  PROJECT {project}\n  RUN_TYPE ENERGY\n&END GLOBAL\n"
    ));
    out.push_str("&FORCE_EVAL\n  METHOD Quickstep\n  &SUBSYS\n    &CELL\n");
    out.push_str(&format!("      ABC {a:.10} {b:.10} {c:.10}\n"));
    out.push_str(&format!(
        "      ALPHA_BETA_GAMMA {alpha:.10} {beta:.10} {gamma:.10}\n"
    ));
    out.push_str("    &END CELL\n    &COORD\n      SCALED\n");
    for atom in &p1.atoms {
        out.push_str(&format!(
            "      {} {:.10} {:.10} {:.10}\n",
            atom.element, atom.frac[0], atom.frac[1], atom.frac[2]
        ));
    }
    out.push_str("    &END COORD\n");
    let mut species: Vec<String> = Vec::new();
    for atom in &p1.atoms {
        if !species.iter().any(|known| known == &atom.element) {
            species.push(atom.element.clone());
        }
    }
    for element in &species {
        out.push_str(&format!(
            "    &KIND {element}\n      ELEMENT {element}\n      BASIS_SET DZVP-MOLOPT-SR-GAUSSIAN\n      POTENTIAL GTH-PBE\n    &END KIND\n"
        ));
    }
    out.push_str("  &END SUBSYS\n  &DFT\n    BASIS_SET_FILE_NAME BASIS_MOLOPT\n");
    out.push_str("    POTENTIAL_FILE_NAME GTH_POTENTIALS\n  &END DFT\n&END FORCE_EVAL\n");
    Ok(out)
}
