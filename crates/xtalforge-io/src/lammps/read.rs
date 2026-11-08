//! LAMMPS data reader (Task 3.3.4).
//!
//! Ported from Crystal-Builder `xtal/io/lammps.py` (MIT). Atoms are read in
//! `atom_style full`; masses name the atom types and bonds become
//! explicit-bond records. Nothing is perceived here.

use std::collections::BTreeMap;

use serde_json::json;
use xtalforge_core::elements::{all, element_by_symbol};
use xtalforge_core::{Element, Error, ExplicitBond, Lattice, Site, Structure, StructureBuilder};

/// One atom row of the `Atoms` section.
#[derive(Debug, Clone, Copy)]
struct AtomRow {
    molecule: usize,
    type_id: u32,
    charge: f64,
    cart: [f64; 3],
}

/// Header fields found before the first section.
#[derive(Debug, Default)]
struct Header {
    atom_count: Option<usize>,
    origin: [f64; 3],
    matrix: [[f64; 3]; 3],
    bounds: [bool; 3],
}

/// The section the scanner is currently inside.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum Section {
    Header,
    Masses,
    Atoms,
    Bonds,
    Other,
}

/// Reads a LAMMPS data file with `atom_style full`.
///
/// # Errors
///
/// Returns [`Error::InvalidInput`] for a missing box or `Atoms` section, an
/// unsupported atom style, a malformed row or a bond index outside the file.
pub fn read_lammps(text: &str) -> Result<Structure, Error> {
    let mut header = Header::default();
    let mut section = Section::Header;
    let mut header_started = false;
    let mut atom_style: Option<String> = None;
    let mut masses: BTreeMap<u32, (f64, Option<String>)> = BTreeMap::new();
    let mut atoms: Vec<(usize, AtomRow)> = Vec::new();
    let mut bonds: Vec<(usize, usize)> = Vec::new();

    for raw in text.lines() {
        let (body, comment) = raw
            .split_once('#')
            .map_or((raw, ""), |(before, after)| (before, after.trim()));
        let words: Vec<&str> = body.split_whitespace().collect();
        let Some(first) = words.first().copied() else {
            continue;
        };
        if let Some(next) = section_of(first) {
            section = next;
            header_started = true;
            if next == Section::Atoms {
                atom_style = atom_style_of(&words, comment);
            }
            continue;
        }
        if header_started && first.bytes().all(|byte| byte.is_ascii_alphabetic()) {
            section = Section::Other;
            continue;
        }
        match section {
            Section::Header => {
                header_started |= parse_header_row(&words, &mut header)?;
            }
            Section::Masses => {
                let type_id = parse_u32(field(&words, 0, "Masses")?, "atom type")?;
                let mass = parse_f64(field(&words, 1, "Masses")?, "mass")?;
                let named = comment
                    .split_whitespace()
                    .next()
                    .and_then(|token| canonical_symbol(token).ok().flatten());
                if masses.insert(type_id, (mass, named)).is_some() {
                    return Err(Error::InvalidInput(format!(
                        "LAMMPS atom type {type_id} is declared twice in Masses"
                    )));
                }
            }
            Section::Atoms => {
                if let Some(style) = &atom_style
                    && !style.eq_ignore_ascii_case("full")
                {
                    return Err(Error::InvalidInput(format!(
                        "LAMMPS atom_style {style:?} is not supported; use full"
                    )));
                }
                if atom_style.is_none() {
                    if words.len() < 7 {
                        return Err(Error::InvalidInput(
                            "LAMMPS atom_style full needs id mol type q x y z".to_owned(),
                        ));
                    }
                    atom_style = Some("full".to_owned());
                }
                let id = parse_usize(field(&words, 0, "Atoms")?, "atom id")?;
                let row = parse_atom_row(&words)?;
                atoms.push((id, row));
            }
            Section::Bonds => {
                let (a, b) = parse_bond_row(&words)?;
                bonds.push((a, b));
            }
            Section::Other => {}
        }
    }

    if atoms.is_empty() {
        return Err(Error::InvalidInput(
            "a LAMMPS data file needs a non-empty Atoms section".to_owned(),
        ));
    }
    for (axis, name) in ["x", "y", "z"].iter().enumerate() {
        let present = header.bounds.get(axis).copied().unwrap_or(false);
        let diagonal = header
            .matrix
            .get(axis)
            .and_then(|row| row.get(axis))
            .copied()
            .unwrap_or(0.0);
        if !present || diagonal <= 0.0 {
            return Err(Error::InvalidInput(format!(
                "LAMMPS box has no positive {name}lo/{name}hi pair"
            )));
        }
    }
    let lattice = Lattice::new(header.matrix)?;
    let (symbols, warnings) = resolve_masses(&masses)?;

    let mut atoms_by_id: BTreeMap<usize, AtomRow> = BTreeMap::new();
    for (id, row) in atoms {
        if id == 0 {
            return Err(Error::InvalidInput(
                "LAMMPS atom ids are one-based".to_owned(),
            ));
        }
        if atoms_by_id.insert(id, row).is_some() {
            return Err(Error::InvalidInput(format!(
                "LAMMPS atom id {id} appears twice"
            )));
        }
    }
    if let Some(expected) = header.atom_count
        && expected != atoms_by_id.len()
    {
        return Err(Error::InvalidInput(format!(
            "LAMMPS header declares {expected} atoms but {} rows were read",
            atoms_by_id.len()
        )));
    }

    let mut index_of_id: BTreeMap<usize, usize> = BTreeMap::new();
    let mut builder = StructureBuilder::new("lammps", lattice);
    let [ox, oy, oz] = header.origin;
    for (index, (id, row)) in atoms_by_id.iter().enumerate() {
        index_of_id.insert(*id, index);
        let symbol = symbols.get(&row.type_id).ok_or_else(|| {
            Error::InvalidInput(format!(
                "LAMMPS atom type {} has no mass in the file",
                row.type_id
            ))
        })?;
        let [x, y, z] = row.cart;
        let frac = Lattice::wrap(lattice.to_frac([x - ox, y - oy, z - oz]));
        let site = Site::new(format!("{symbol}{id}"), symbol, frac)?
            .with_charge(row.charge)?
            .with_prop("lammps_molecule", json!(row.molecule));
        builder = builder.site(site);
    }
    for (a, b) in bonds {
        let i = index_of_id
            .get(&a)
            .copied()
            .ok_or_else(|| Error::InvalidInput(format!("LAMMPS bond names missing atom id {a}")))?;
        let j = index_of_id
            .get(&b)
            .copied()
            .ok_or_else(|| Error::InvalidInput(format!("LAMMPS bond names missing atom id {b}")))?;
        if i == j {
            return Err(Error::InvalidInput(format!(
                "LAMMPS bond joins atom {a} to itself"
            )));
        }
        builder = builder.explicit_bond(ExplicitBond {
            i,
            j,
            image: [0, 0, 0],
            order: 1.0,
            remove: false,
        });
    }
    if !warnings.is_empty() {
        builder = builder.metadata("warnings", json!(warnings));
    }
    builder.build()
}

fn section_of(word: &str) -> Option<Section> {
    match word.to_ascii_lowercase().as_str() {
        "masses" => Some(Section::Masses),
        "atoms" => Some(Section::Atoms),
        "bonds" => Some(Section::Bonds),
        "angles" | "dihedrals" | "impropers" | "velocities" | "ellipsoids" | "lines"
        | "triangles" | "bodies" => Some(Section::Other),
        _ => None,
    }
}

fn atom_style_of(words: &[&str], comment: &str) -> Option<String> {
    for token in comment.split_whitespace() {
        let lower = token.to_ascii_lowercase();
        if matches!(
            lower.as_str(),
            "full" | "charge" | "molecular" | "atomic" | "bond"
        ) {
            return Some(lower);
        }
    }
    words.get(1).map(|word| word.to_ascii_lowercase())
}

fn parse_header_row(words: &[&str], header: &mut Header) -> Result<bool, Error> {
    let first = field(words, 0, "header")?;
    if words
        .get(1)
        .is_some_and(|word| word.eq_ignore_ascii_case("atoms"))
    {
        header.atom_count = Some(parse_usize(first, "atom count")?);
        return Ok(true);
    }
    if words.len() >= 3
        && words
            .get(1)
            .is_some_and(|word| word.eq_ignore_ascii_case("atom"))
        && words
            .get(2)
            .is_some_and(|word| word.eq_ignore_ascii_case("types"))
    {
        return Ok(true);
    }
    if words.len() >= 4 {
        let high_text = words.get(1).copied().unwrap_or("");
        let label = words.get(2).copied().unwrap_or("");
        let high_label = words.get(3).copied().unwrap_or("");
        let axis = match label {
            "xlo" if high_label == "xhi" => Some(0),
            "ylo" if high_label == "yhi" => Some(1),
            "zlo" if high_label == "zhi" => Some(2),
            _ => None,
        };
        if let Some(axis) = axis {
            let low = parse_f64(first, "box lower bound")?;
            let high = parse_f64(high_text, "box upper bound")?;
            if let Some(slot) = header.origin.get_mut(axis) {
                *slot = low;
            }
            if let Some(row) = header.matrix.get_mut(axis)
                && let Some(diagonal) = row.get_mut(axis)
            {
                *diagonal = high - low;
            }
            if let Some(bound) = header.bounds.get_mut(axis) {
                *bound = true;
            }
            return Ok(true);
        }
    }
    if words.len() >= 6
        && words.get(3).is_some_and(|word| *word == "xy")
        && words.get(4).is_some_and(|word| *word == "xz")
        && words.get(5).is_some_and(|word| *word == "yz")
    {
        let xy = parse_f64(first, "xy tilt")?;
        let xz = parse_f64(field(words, 1, "tilt")?, "xz tilt")?;
        let yz = parse_f64(field(words, 2, "tilt")?, "yz tilt")?;
        if let Some(row) = header.matrix.get_mut(1)
            && let Some(slot) = row.get_mut(0)
        {
            *slot = xy;
        }
        if let Some(row) = header.matrix.get_mut(2) {
            if let Some(slot) = row.get_mut(0) {
                *slot = xz;
            }
            if let Some(slot) = row.get_mut(1) {
                *slot = yz;
            }
        }
        return Ok(true);
    }
    Ok(false)
}

fn parse_atom_row(words: &[&str]) -> Result<AtomRow, Error> {
    if words.len() < 7 {
        return Err(Error::InvalidInput(
            "LAMMPS atom_style full needs id mol type q x y z".to_owned(),
        ));
    }
    Ok(AtomRow {
        molecule: parse_usize(field(words, 1, "Atoms")?, "molecule id")?,
        type_id: parse_u32(field(words, 2, "Atoms")?, "atom type")?,
        charge: parse_f64(field(words, 3, "Atoms")?, "charge")?,
        cart: [
            parse_f64(field(words, 4, "Atoms")?, "x")?,
            parse_f64(field(words, 5, "Atoms")?, "y")?,
            parse_f64(field(words, 6, "Atoms")?, "z")?,
        ],
    })
}

fn parse_bond_row(words: &[&str]) -> Result<(usize, usize), Error> {
    if words.len() < 4 {
        return Err(Error::InvalidInput(
            "a LAMMPS Bonds row needs id type atom1 atom2".to_owned(),
        ));
    }
    parse_usize(field(words, 0, "Bonds")?, "bond id")?;
    parse_u32(field(words, 1, "Bonds")?, "bond type")?;
    Ok((
        parse_usize(field(words, 2, "Bonds")?, "first atom id")?,
        parse_usize(field(words, 3, "Bonds")?, "second atom id")?,
    ))
}

fn resolve_masses(
    masses: &BTreeMap<u32, (f64, Option<String>)>,
) -> Result<(BTreeMap<u32, String>, Vec<String>), Error> {
    let mut symbols = BTreeMap::new();
    let mut warnings = Vec::new();
    for (type_id, (mass, named)) in masses {
        let symbol = match named {
            Some(symbol) => element_by_symbol(symbol)?.symbol.clone(),
            None => {
                let (symbol, warning) = element_for_mass(*mass)?;
                if let Some(warning) = warning {
                    warnings.push(format!("atom type {type_id}: {warning}"));
                }
                symbol
            }
        };
        symbols.insert(*type_id, symbol);
    }
    Ok((symbols, warnings))
}

fn element_for_mass(mass: f64) -> Result<(String, Option<String>), Error> {
    if !mass.is_finite() || mass <= 0.0 {
        return Err(Error::InvalidInput(format!(
            "LAMMPS mass must be positive and finite, got {mass}"
        )));
    }
    let mut best: Option<(&'static Element, f64)> = None;
    let mut second: Option<(&'static Element, f64)> = None;
    for element in all()? {
        if element.symbol == "X" {
            continue;
        }
        let distance = (element.mass - mass).abs();
        match best {
            None => best = Some((element, distance)),
            Some((_, best_distance)) if distance < best_distance => {
                second = best;
                best = Some((element, distance));
            }
            _ => {
                if second.is_none_or(|(_, best_second)| distance < best_second) {
                    second = Some((element, distance));
                }
            }
        }
    }
    let Some((closest, closest_distance)) = best else {
        return Err(Error::InvalidInput("the element table is empty".to_owned()));
    };
    let warning = second.and_then(|(other, other_distance)| {
        (other_distance - closest_distance < 0.5).then(|| {
            format!(
                "mass {mass:.3} is ambiguous between {} ({:.3}) and {} ({:.3}); using {}",
                closest.symbol, closest.mass, other.symbol, other.mass, closest.symbol
            )
        })
    });
    Ok((closest.symbol.clone(), warning))
}

fn canonical_symbol(text: &str) -> Result<Option<String>, Error> {
    xtalforge_core::elements::canonical_symbol(text)
}

fn field<'a>(words: &'a [&'a str], index: usize, section: &str) -> Result<&'a str, Error> {
    words
        .get(index)
        .copied()
        .ok_or_else(|| Error::InvalidInput(format!("the LAMMPS {section} row is missing a field")))
}

fn parse_f64(text: &str, what: &str) -> Result<f64, Error> {
    text.parse::<f64>()
        .map_err(|_| Error::InvalidInput(format!("bad LAMMPS {what}: {text:?}")))
}

fn parse_u32(text: &str, what: &str) -> Result<u32, Error> {
    text.parse::<u32>()
        .map_err(|_| Error::InvalidInput(format!("bad LAMMPS {what}: {text:?}")))
}

fn parse_usize(text: &str, what: &str) -> Result<usize, Error> {
    text.parse::<usize>()
        .map_err(|_| Error::InvalidInput(format!("bad LAMMPS {what}: {text:?}")))
}
