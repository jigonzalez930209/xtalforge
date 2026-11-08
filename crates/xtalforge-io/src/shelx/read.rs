//! SHELX instruction-file reader (Task 3.3.5).
//!
//! Ported from Crystal-Builder `xtal/io/shelx.py` (MIT). `LATT` and `SYMM`
//! are recorded as metadata but not applied: the atoms are read in P1.

use serde_json::json;
use xtalforge_core::elements::element_by_symbol;
use xtalforge_core::{Error, Lattice, Site, Structure, StructureBuilder};

use super::expr::{resolve_field, resolve_token};

/// One raw atom record, kept until `CELL` and `SFAC` are known.
#[derive(Debug, Clone)]
struct AtomLine {
    label: String,
    element: Option<String>,
    sfac: Option<usize>,
    values: Vec<String>,
    part: i32,
    default_sof: Option<f64>,
}

/// Reads a SHELX instruction file.
///
/// # Errors
///
/// Returns [`Error::InvalidInput`] for a missing or malformed `CELL`, an atom
/// that names an undeclared `SFAC`, an unresolvable free variable or a bad
/// occupancy.
pub fn read_shelx(text: &str) -> Result<Structure, Error> {
    let mut title = String::from("shelx");
    let mut cell: Option<[f64; 6]> = None;
    let mut latt = 1i32;
    let mut symm: Vec<String> = Vec::new();
    let mut sfac: Vec<String> = Vec::new();
    let mut fvar: Vec<f64> = Vec::new();
    let mut atoms: Vec<AtomLine> = Vec::new();
    let mut part = 0i32;
    let mut part_sof: Option<f64> = None;

    for line in logical_lines(text) {
        let words: Vec<&str> = line.split_whitespace().collect();
        let Some(first) = words.first().copied() else {
            continue;
        };
        let head = first.to_ascii_uppercase();
        let keyword = head.split('_').next().unwrap_or("");
        if keyword == "HKLF" || keyword == "END" {
            break;
        }
        let rest: &[&str] = words.get(1..).unwrap_or(&[]);
        match keyword {
            "TITL" => {
                title = rest.join(" ");
                if title.trim().is_empty() {
                    title = String::from("shelx");
                }
            }
            "CELL" => cell = Some(parse_cell(rest)?),
            "LATT" => {
                let value = parse_i32(field(&words, 1, "LATT")?, "LATT")?;
                if !(1..=7).contains(&value.abs()) {
                    return Err(Error::InvalidInput(format!(
                        "SHELX LATT must be between -7 and 7, got {value}"
                    )));
                }
                latt = value;
            }
            "SYMM" if !rest.is_empty() => symm.push(rest.join(" ")),
            "SFAC" => sfac.extend(parse_sfac(rest)?),
            "FVAR" => {
                for token in rest {
                    fvar.push(parse_f64(token, "FVAR value")?);
                }
            }
            "PART" => {
                part = words
                    .get(1)
                    .map(|token| parse_i32(token, "PART number"))
                    .transpose()?
                    .unwrap_or(0);
                part_sof = words
                    .get(2)
                    .map(|token| parse_f64(token, "PART occupancy"))
                    .transpose()?;
            }
            _ if !is_instruction(keyword) => {
                if let Some(atom) = parse_atom_line(&words, part, part_sof) {
                    atoms.push(atom);
                }
            }
            _ => {}
        }
    }

    let cell = cell.ok_or_else(|| {
        Error::InvalidInput("a SHELX instruction file must have a CELL line".to_owned())
    })?;
    if atoms.is_empty() {
        return Err(Error::InvalidInput(
            "no SHELX atom records were found".to_owned(),
        ));
    }
    let [cell_a, cell_b, cell_c, cell_alpha, cell_beta, cell_gamma] = cell;
    let lattice =
        Lattice::from_parameters(cell_a, cell_b, cell_c, cell_alpha, cell_beta, cell_gamma)?;

    let mut builder = StructureBuilder::new(title, lattice);
    for atom in atoms {
        let element = match (&atom.element, atom.sfac) {
            (Some(element), _) => element.clone(),
            (None, Some(number)) => {
                let index = number.checked_sub(1).ok_or_else(|| {
                    Error::InvalidInput("SHELX SFAC numbers are one-based".to_owned())
                })?;
                sfac.get(index).cloned().ok_or_else(|| {
                    Error::InvalidInput(format!(
                        "SHELX atom {} names SFAC {number}, but only {} types are declared",
                        atom.label,
                        sfac.len()
                    ))
                })?
            }
            (None, None) => {
                return Err(Error::InvalidInput(format!(
                    "SHELX atom {} has no SFAC number",
                    atom.label
                )));
            }
        };
        let x = resolve_token(value_at(&atom.values, 0, &atom.label)?, &fvar)?;
        let y = resolve_token(value_at(&atom.values, 1, &atom.label)?, &fvar)?;
        let z = resolve_token(value_at(&atom.values, 2, &atom.label)?, &fvar)?;
        let occupancy = match atom.values.get(3) {
            Some(token) => resolve_token(token, &fvar)?,
            None => resolve_field(atom.default_sof.unwrap_or(11.0), &fvar),
        };
        if !(0.0..=1.0).contains(&occupancy) {
            return Err(Error::InvalidInput(format!(
                "SHELX atom {} has occupancy {occupancy}, outside [0, 1]",
                atom.label
            )));
        }
        let mut site =
            Site::new(atom.label.clone(), &element, [x, y, z])?.with_occupancy(occupancy)?;
        if atom.values.len() < 10
            && let Some(token) = atom.values.get(4)
        {
            let u_iso = resolve_token(token, &fvar)?;
            if u_iso.is_finite() && u_iso >= 0.0 {
                site = site.with_u_iso(u_iso)?;
            }
        }
        let mut group = atom.part;
        if atom.values.len() == 6
            && let Some(token) = atom.values.get(5)
            && let Ok(value) = token.parse::<i32>()
        {
            group = value;
        }
        if group != 0 {
            site = site.with_disorder_group(group.unsigned_abs());
        }
        builder = builder.site(site);
    }

    builder = builder.metadata("shelx_latt", json!(latt));
    if !symm.is_empty() {
        builder = builder.metadata("shelx_symm", json!(symm));
    }
    if latt.abs() != 1 || !symm.is_empty() {
        builder = builder.metadata(
            "symmetry_warning",
            json!("SHELX LATT/SYMM symmetry is recorded but not applied; atoms are read in P1"),
        );
    }
    builder.build()
}

fn logical_lines(text: &str) -> Vec<String> {
    let mut lines = Vec::new();
    let mut pending = String::new();
    for raw in text.lines() {
        let body = raw.split_once('!').map_or(raw, |(before, _)| before);
        let trimmed = body.trim_end();
        if trimmed
            .split_whitespace()
            .next()
            .is_some_and(|word| word.eq_ignore_ascii_case("REM"))
        {
            continue;
        }
        if let Some(continuation) = trimmed.strip_suffix('=') {
            pending.push_str(continuation);
            pending.push(' ');
            continue;
        }
        let mut full = std::mem::take(&mut pending);
        full.push_str(trimmed);
        lines.push(full);
    }
    if !pending.is_empty() {
        lines.push(pending);
    }
    lines
}

fn parse_atom_line(words: &[&str], part: i32, default_sof: Option<f64>) -> Option<AtomLine> {
    if words.len() < 5 {
        return None;
    }
    let first = *words.first()?;
    let second = *words.get(1)?;
    if let Ok(number) = second.parse::<usize>()
        && number > 0
    {
        let values = words
            .get(2..)?
            .iter()
            .map(|word| (*word).to_owned())
            .collect();
        return Some(AtomLine {
            label: first.to_owned(),
            element: None,
            sfac: Some(number),
            values,
            part,
            default_sof,
        });
    }
    let element = element_by_symbol(first).ok()?.symbol.clone();
    let values: Vec<String> = words
        .get(1..)?
        .iter()
        .map(|word| (*word).to_owned())
        .collect();
    if values.len() < 3 || values.iter().take(3).any(|value| !is_value_token(value)) {
        return None;
    }
    Some(AtomLine {
        label: first.to_owned(),
        element: Some(element),
        sfac: None,
        values,
        part,
        default_sof,
    })
}

fn parse_cell(words: &[&str]) -> Result<[f64; 6], Error> {
    let mut values = Vec::with_capacity(words.len());
    for word in words {
        values.push(parse_f64(word, "CELL value")?);
    }
    if values.len() < 6 {
        return Err(Error::InvalidInput(
            "SHELX CELL needs six numbers".to_owned(),
        ));
    }
    let skip = values.len().saturating_sub(6);
    let six: Vec<f64> = values.into_iter().skip(skip).collect();
    let mut cell = [0.0; 6];
    for (slot, value) in cell.iter_mut().zip(six) {
        *slot = value;
    }
    Ok(cell)
}

fn parse_sfac(words: &[&str]) -> Result<Vec<String>, Error> {
    let Some(first) = words.first().copied() else {
        return Err(Error::InvalidInput(
            "SHELX SFAC needs an element".to_owned(),
        ));
    };
    let long_form = words.iter().skip(1).any(|word| word.parse::<f64>().is_ok());
    if long_form {
        return Ok(vec![element_symbol(first)?]);
    }
    words.iter().map(|word| element_symbol(word)).collect()
}

fn element_symbol(text: &str) -> Result<String, Error> {
    Ok(element_by_symbol(text)?.symbol.clone())
}

fn is_value_token(token: &str) -> bool {
    !token.is_empty()
        && token.bytes().all(|byte| {
            byte.is_ascii_alphanumeric()
                || matches!(byte, b'.' | b'+' | b'-' | b'*' | b'/' | b'(' | b')' | b'_')
        })
}

fn field<'a>(words: &[&'a str], index: usize, keyword: &str) -> Result<&'a str, Error> {
    words
        .get(index)
        .copied()
        .ok_or_else(|| Error::InvalidInput(format!("SHELX {keyword} line is missing a value")))
}

fn value_at<'a>(values: &'a [String], index: usize, label: &str) -> Result<&'a str, Error> {
    values
        .get(index)
        .map(String::as_str)
        .ok_or_else(|| Error::InvalidInput(format!("SHELX atom {label} needs three coordinates")))
}

fn parse_f64(text: &str, what: &str) -> Result<f64, Error> {
    text.parse::<f64>()
        .map_err(|_| Error::InvalidInput(format!("bad SHELX {what}: {text:?}")))
}

fn parse_i32(text: &str, what: &str) -> Result<i32, Error> {
    text.parse::<i32>()
        .map_err(|_| Error::InvalidInput(format!("bad SHELX {what}: {text:?}")))
}

fn is_instruction(word: &str) -> bool {
    matches!(
        word,
        "ABIN"
            | "ACTA"
            | "AFIX"
            | "ANIS"
            | "ANSC"
            | "ANSR"
            | "BASF"
            | "BEDE"
            | "BIND"
            | "BLOC"
            | "BOND"
            | "BUMP"
            | "CELL"
            | "CGLS"
            | "CHIV"
            | "CONF"
            | "CONN"
            | "DAMP"
            | "DANG"
            | "DEFS"
            | "DELU"
            | "DFIX"
            | "DISP"
            | "EADP"
            | "END"
            | "EQIV"
            | "EXTI"
            | "EXYZ"
            | "FEND"
            | "FLAT"
            | "FMAP"
            | "FRAG"
            | "FREE"
            | "FVAR"
            | "GRID"
            | "HFIX"
            | "HKLF"
            | "HOPE"
            | "HTAB"
            | "ISOR"
            | "L.S."
            | "LATT"
            | "LAUE"
            | "LIST"
            | "LONE"
            | "MERG"
            | "MOLE"
            | "MORE"
            | "MOVE"
            | "MPLA"
            | "NCSY"
            | "NEUT"
            | "OMIT"
            | "PART"
            | "PLAN"
            | "PRIG"
            | "REM"
            | "RESI"
            | "RIGU"
            | "RTAB"
            | "SADI"
            | "SAME"
            | "SFAC"
            | "SHEL"
            | "SIMU"
            | "SIZE"
            | "SPEC"
            | "STIR"
            | "SUMP"
            | "SWAT"
            | "SYMM"
            | "TEMP"
            | "TIME"
            | "TITL"
            | "TWIN"
            | "TWST"
            | "UNIT"
            | "WGHT"
            | "WIGL"
            | "WPDB"
            | "XNPD"
            | "ZERR"
    )
}
