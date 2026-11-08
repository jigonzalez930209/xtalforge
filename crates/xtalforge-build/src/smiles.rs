//! SMILES graph parser (Task 6.1.1).

use serde::{Deserialize, Serialize};
use xtalforge_core::Error;

/// One atom of a parsed SMILES graph.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct SmilesAtom {
    /// Element symbol.
    pub element: String,
    /// Whether the atom is aromatic.
    pub aromatic: bool,
    /// Formal charge.
    pub charge: i32,
}

/// One bond of the graph.
#[derive(Debug, Clone, Copy, PartialEq, Serialize, Deserialize)]
pub struct SmilesBond {
    /// First atom.
    pub a: usize,
    /// Second atom.
    pub b: usize,
    /// Bond order (`1`, `2`, `3` or `1.5` for aromatic).
    pub order: f64,
}

/// A parsed molecular graph.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct SmilesMolecule {
    /// Atoms in parse order.
    pub atoms: Vec<SmilesAtom>,
    /// Bonds.
    pub bonds: Vec<SmilesBond>,
}

/// Parses the organic subset of SMILES into a graph.
///
/// Supports atoms `B C N O P S F Cl Br I`, aromatic lower-case atoms, branches
/// `()`, ring closures `1`-`9`, bond symbols `- = # : / \` and charge in
/// brackets (`[NH4+]`). The 3D embedding (Task 6.1.2) is not implemented here.
///
/// # Errors
///
/// Returns [`Error::InvalidInput`] for an empty string, an unbalanced branch or
/// ring closure, or an unknown atom.
pub fn parse_smiles(text: &str) -> Result<SmilesMolecule, Error> {
    if text.trim().is_empty() {
        return Err(Error::InvalidInput("SMILES is empty".to_owned()));
    }
    let mut atoms = Vec::new();
    let mut bonds = Vec::new();
    let mut stack = Vec::new();
    let mut previous: Option<usize> = None;
    let mut rings: std::collections::BTreeMap<char, usize> = std::collections::BTreeMap::new();
    let mut order = 1.0;
    let mut characters = text.chars().peekable();
    while let Some(character) = characters.next() {
        match character {
            ' ' | '\t' => {}
            '(' => {
                let Some(current) = previous else {
                    return Err(Error::InvalidInput("branch without an atom".to_owned()));
                };
                stack.push(current);
                order = 1.0;
            }
            ')' => {
                let Some(parent) = stack.pop() else {
                    return Err(Error::InvalidInput("unbalanced branch".to_owned()));
                };
                previous = Some(parent);
                order = 1.0;
            }
            '-' => order = 1.0,
            '=' => order = 2.0,
            '#' => order = 3.0,
            ':' => order = 1.5,
            '/' | '\\' => order = 1.0,
            '.' => {
                if previous.is_none() || !rings.is_empty() {
                    return Err(Error::InvalidInput(
                        "a component separator needs a finished component".to_owned(),
                    ));
                }
                previous = None;
                order = 1.0;
            }
            '1'..='9' => {
                let Some(current) = previous else {
                    return Err(Error::InvalidInput(
                        "ring closure without an atom".to_owned(),
                    ));
                };
                if let Some(other) = rings.remove(&character) {
                    bonds.push(SmilesBond {
                        a: other,
                        b: current,
                        order,
                    });
                    order = 1.0;
                } else {
                    rings.insert(character, current);
                }
            }
            '[' => {
                let mut inner = String::new();
                for next in characters.by_ref() {
                    if next == ']' {
                        break;
                    }
                    inner.push(next);
                }
                let atom = bracket_atom(&inner)?;
                let index = atoms.len();
                atoms.push(atom);
                attach(previous, index, order, &mut bonds);
                previous = Some(index);
                order = 1.0;
            }
            _ => {
                let (atom, width) = organic_atom(character, &mut characters)?;
                let index = atoms.len();
                atoms.push(atom);
                attach(previous, index, order, &mut bonds);
                previous = Some(index);
                order = 1.0;
                let _ = width;
            }
        }
    }
    if !stack.is_empty() || !rings.is_empty() {
        return Err(Error::InvalidInput(
            "unbalanced branch or ring closure".to_owned(),
        ));
    }
    if atoms.is_empty() {
        return Err(Error::InvalidInput("SMILES has no atoms".to_owned()));
    }
    Ok(SmilesMolecule { atoms, bonds })
}

fn attach(previous: Option<usize>, index: usize, order: f64, bonds: &mut Vec<SmilesBond>) {
    if let Some(previous) = previous {
        bonds.push(SmilesBond {
            a: previous,
            b: index,
            order,
        });
    }
}

fn organic_atom(
    first: char,
    characters: &mut std::iter::Peekable<std::str::Chars<'_>>,
) -> Result<(SmilesAtom, usize), Error> {
    if matches!(first, 'C' | 'c') && characters.peek() == Some(&'l') {
        characters.next();
        return Ok((
            SmilesAtom {
                element: "Cl".to_owned(),
                aromatic: false,
                charge: 0,
            },
            2,
        ));
    }
    if matches!(first, 'B' | 'b') && characters.peek() == Some(&'r') {
        characters.next();
        return Ok((
            SmilesAtom {
                element: "Br".to_owned(),
                aromatic: false,
                charge: 0,
            },
            2,
        ));
    }
    let (element, aromatic) = match first {
        'B' | 'C' | 'N' | 'O' | 'P' | 'S' | 'F' | 'I' => (first.to_string(), false),
        'b' | 'c' | 'n' | 'o' | 'p' | 's' => (first.to_ascii_uppercase().to_string(), true),
        other => {
            return Err(Error::InvalidInput(format!(
                "unsupported SMILES atom {other:?}"
            )));
        }
    };
    Ok((
        SmilesAtom {
            element,
            aromatic,
            charge: 0,
        },
        1,
    ))
}

fn bracket_atom(inner: &str) -> Result<SmilesAtom, Error> {
    let mut characters = inner.chars().peekable();
    let first = characters
        .next()
        .ok_or_else(|| Error::InvalidInput("empty bracket atom".to_owned()))?;
    if !first.is_ascii_alphabetic() {
        return Err(Error::InvalidInput("bad bracket element".to_owned()));
    }
    let aromatic = first.is_ascii_lowercase();
    let mut element = String::new();
    element.push(if aromatic {
        first.to_ascii_uppercase()
    } else {
        first
    });
    while element.len() < 2 {
        let Some(next) = characters.peek().copied() else {
            break;
        };
        if !next.is_ascii_lowercase() {
            break;
        }
        element.push(next);
        characters.next();
    }
    if !element.chars().all(|value| value.is_ascii_alphabetic()) {
        return Err(Error::InvalidInput("bad bracket element".to_owned()));
    }
    let mut charge = 0;
    let mut sign = 0;
    for character in characters {
        match character {
            '+' => {
                if sign == 0 {
                    charge += 1;
                }
                sign = 1;
            }
            '-' => {
                if sign == 0 {
                    charge -= 1;
                }
                sign = -1;
            }
            value if value.is_ascii_digit() && sign != 0 => {
                charge += sign * value.to_digit(10).unwrap_or(0) as i32;
                sign = 0;
            }
            _ => {}
        }
    }
    Ok(SmilesAtom {
        element,
        aromatic,
        charge,
    })
}
