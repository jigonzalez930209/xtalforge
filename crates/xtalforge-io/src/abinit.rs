//! ABINIT input reader and writer (Task 13.2.2).
//!
//! Reads one variable assignment per line: `natom`, `ntypat`, `typat`,
//! `znucl`, `acell` (bohr, `angstrom` or `bohr` unit token), `rprim`, `xred`
//! and `xcart` (bohr by default). Repetition such as `3*10.0` is expanded.
//! Multi-line assignment continuation is not part of this documented subset.

#![allow(clippy::indexing_slicing)] // fixed-format parser, lengths checked

use serde_json::Value;
use xtalforge_core::elements::element_by_z;
use xtalforge_core::{Error, Lattice, Site, Structure, StructureBuilder};

use super::pwscf::BOHR;
use crate::molecule::parse_number;

/// Reads an ABINIT input file.
///
/// # Errors
///
/// Returns [`Error::InvalidInput`] for an unknown unit, a missing `acell` or
/// cell definition, inconsistent atom counts or a malformed number.
pub fn read_abinit(text: &str) -> Result<Structure, Error> {
    let assignments = assignments(text);
    let natom = integer(&assignments, "natom")?
        .ok_or_else(|| Error::InvalidInput("ABINIT input has no natom".to_owned()))?;
    let ntypat = integer(&assignments, "ntypat")?.unwrap_or(1);
    let typat = values(&assignments, "typat")?
        .ok_or_else(|| Error::InvalidInput("ABINIT input has no typat".to_owned()))?;
    if typat.len() < natom {
        return Err(Error::InvalidInput(format!(
            "ABINIT typat has {} entries for natom={natom}",
            typat.len()
        )));
    }
    let znucl = values(&assignments, "znucl")?
        .ok_or_else(|| Error::InvalidInput("ABINIT input has no znucl".to_owned()))?;
    if znucl.len() < ntypat {
        return Err(Error::InvalidInput(format!(
            "ABINIT znucl has {} entries for ntypat={ntypat}",
            znucl.len()
        )));
    }
    let (acell, acell_unit) = values_unit(&assignments, "acell")?
        .ok_or_else(|| Error::InvalidInput("ABINIT input has no acell".to_owned()))?;
    let acell_scale = match acell_unit.as_deref() {
        None | Some("bohr") | Some("au") => BOHR,
        Some("angstrom") | Some("ang") => 1.0,
        Some(other) => {
            return Err(Error::InvalidInput(format!(
                "unsupported ABINIT acell unit {other:?}"
            )));
        }
    };
    let acell: Vec<f64> = acell.iter().map(|value| value * acell_scale).collect();
    let rprim = values(&assignments, "rprim")?.unwrap_or_else(|| {
        vec![
            1.0, 0.0, 0.0, //
            0.0, 1.0, 0.0, //
            0.0, 0.0, 1.0,
        ]
    });
    if rprim.len() < 9 {
        return Err(Error::InvalidInput(
            "ABINIT rprim needs nine values".to_owned(),
        ));
    }
    let scale = |index: usize| acell.get(index).copied().unwrap_or(1.0);
    let lattice = Lattice::new([
        [
            rprim[0] * scale(0),
            rprim[1] * scale(0),
            rprim[2] * scale(0),
        ],
        [
            rprim[3] * scale(1),
            rprim[4] * scale(1),
            rprim[5] * scale(1),
        ],
        [
            rprim[6] * scale(2),
            rprim[7] * scale(2),
            rprim[8] * scale(2),
        ],
    ])?;
    let (positions, cartesian) = if let Some((values, unit)) = values_unit(&assignments, "xcart")? {
        let scale = match unit.as_deref() {
            None | Some("bohr") | Some("au") => BOHR,
            Some("angstrom") | Some("ang") => 1.0,
            Some(other) => {
                return Err(Error::InvalidInput(format!(
                    "unsupported ABINIT xcart unit {other:?}"
                )));
            }
        };
        let points: Vec<[f64; 3]> = values
            .chunks(3)
            .filter_map(|chunk| {
                Some([
                    chunk.first()? * scale,
                    chunk.get(1)? * scale,
                    chunk.get(2)? * scale,
                ])
            })
            .collect();
        (points, true)
    } else if let Some((values, _)) = values_unit(&assignments, "xred")? {
        let points: Vec<[f64; 3]> = values
            .chunks(3)
            .filter_map(|chunk| Some([*chunk.first()?, *chunk.get(1)?, *chunk.get(2)?]))
            .collect();
        (points, false)
    } else {
        return Err(Error::InvalidInput(
            "ABINIT input needs xred or xcart".to_owned(),
        ));
    };
    if positions.len() < natom {
        return Err(Error::InvalidInput(format!(
            "ABINIT has {} positions for natom={natom}",
            positions.len()
        )));
    }
    let mut builder = StructureBuilder::new("abinit", lattice);
    for (index, position) in positions.iter().take(natom).enumerate() {
        let type_index = typat.get(index).copied().unwrap_or(1.0) as usize;
        let z = znucl
            .get(type_index.saturating_sub(1))
            .copied()
            .unwrap_or(0.0) as u8;
        let element = element_by_z(z)?.symbol.clone();
        let frac = if cartesian {
            lattice.to_frac(*position)
        } else {
            *position
        };
        builder = builder.site(Site::new(
            format!("{element}{}", index + 1),
            &element,
            frac,
        )?);
    }
    builder = builder.metadata("abinit_ntypat", Value::from(ntypat));
    builder.build()
}

/// Writes a P1 structure as an ABINIT input file.
///
/// # Errors
///
/// Returns an error for an empty structure or a failed P1 expansion.
pub fn write_abinit(structure: &Structure) -> Result<String, Error> {
    let p1 = structure.p1()?;
    if p1.atoms.is_empty() {
        return Err(Error::InvalidInput(
            "an ABINIT input cannot hold a structure with no atoms".to_owned(),
        ));
    }
    let mut species: Vec<String> = Vec::new();
    let mut typat: Vec<usize> = Vec::new();
    for atom in &p1.atoms {
        let index = match species.iter().position(|known| known == &atom.element) {
            Some(index) => index,
            None => {
                species.push(atom.element.clone());
                species.len() - 1
            }
        };
        typat.push(index + 1);
    }
    let mut out = String::from("# ABINIT input written by xtalforge\n");
    out.push_str(&format!("natom {}\n", p1.atoms.len()));
    out.push_str(&format!("ntypat {}\n", species.len()));
    out.push_str(&format!(
        "typat {}\n",
        typat
            .iter()
            .map(usize::to_string)
            .collect::<Vec<_>>()
            .join(" ")
    ));
    let znucl: Vec<String> = species
        .iter()
        .map(|symbol| {
            xtalforge_core::elements::atomic_number(symbol)
                .map(|z| z.to_string())
                .unwrap_or_else(|_| "0".to_owned())
        })
        .collect();
    out.push_str(&format!("znucl {}\n", znucl.join(" ")));
    out.push_str("acell 1.0 1.0 1.0 angstrom\n");
    let matrix = structure.lattice().matrix();
    let rprim: Vec<String> = matrix
        .iter()
        .flatten()
        .map(|value| format!("{value:.10}"))
        .collect();
    out.push_str(&format!("rprim {}\n", rprim.join(" ")));
    let xred: Vec<String> = p1
        .atoms
        .iter()
        .flat_map(|atom| atom.frac.iter().map(|value| format!("{value:.10}")))
        .collect();
    out.push_str(&format!("xred {}\n", xred.join(" ")));
    out.push_str("ecut 20.0\n");
    Ok(out)
}

struct Assignment {
    key: String,
    tokens: Vec<String>,
}

fn assignments(text: &str) -> Vec<Assignment> {
    let mut out = Vec::new();
    for line in text.lines() {
        let line = line.split('#').next().unwrap_or("").trim();
        if line.is_empty() {
            continue;
        }
        let mut parts = line.split_whitespace();
        let Some(key) = parts.next() else {
            continue;
        };
        out.push(Assignment {
            key: key.to_ascii_lowercase(),
            tokens: parts.map(str::to_owned).collect(),
        });
    }
    out
}

fn tokens<'a>(assignments: &'a [Assignment], key: &str) -> Option<&'a [String]> {
    assignments
        .iter()
        .find(|assignment| assignment.key == key)
        .map(|assignment| assignment.tokens.as_slice())
}

fn values(assignments: &[Assignment], key: &str) -> Result<Option<Vec<f64>>, Error> {
    values_unit(assignments, key).map(|found| found.map(|(values, _)| values))
}

fn integer(assignments: &[Assignment], key: &str) -> Result<Option<usize>, Error> {
    let Some(values) = values(assignments, key)? else {
        return Ok(None);
    };
    let first = values
        .first()
        .copied()
        .ok_or_else(|| Error::InvalidInput(format!("ABINIT {key} has no value")))?;
    Ok(Some(first as usize))
}

type ValuesUnit = (Vec<f64>, Option<String>);

fn values_unit(assignments: &[Assignment], key: &str) -> Result<Option<ValuesUnit>, Error> {
    let Some(tokens) = tokens(assignments, key) else {
        return Ok(None);
    };
    let mut values = Vec::new();
    let mut unit = None;
    for token in tokens {
        if let Some((count, value)) = token.split_once('*') {
            let count: usize = count
                .parse()
                .map_err(|_| Error::InvalidInput(format!("bad ABINIT repetition {token:?}")))?;
            let value = parse_number(value, "ABINIT value")?;
            for _ in 0..count {
                values.push(value);
            }
        } else if let Ok(value) = token.parse::<f64>() {
            values.push(value);
        } else if token
            .chars()
            .next()
            .is_some_and(|c| c.is_ascii_alphabetic())
        {
            unit = Some(token.to_ascii_lowercase());
        } else {
            return Err(Error::InvalidInput(format!(
                "ABINIT {key} has a bad token {token:?}"
            )));
        }
    }
    for value in &values {
        if !value.is_finite() {
            return Err(Error::InvalidInput(format!(
                "ABINIT {key} has a non-finite value"
            )));
        }
    }
    Ok(Some((values, unit)))
}
