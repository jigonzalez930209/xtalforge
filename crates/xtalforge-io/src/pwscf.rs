//! Quantum ESPRESSO `pw.x` input and output (Task 13.2.1).
//!
//! The reader accepts an input deck (namelist variables plus
//! `ATOMIC_SPECIES`, `CELL_PARAMETERS` and `ATOMIC_POSITIONS` cards) and the
//! text output of `pw.x` (crystal axes, Cartesian or crystal positions and the
//! final total energy). The writer emits a self-consistent `ibrav=0` input
//! deck; an output file is read-only because it describes a finished run.
//!
//! Supported input cells are `ibrav=0` with `CELL_PARAMETERS`, and
//! `ibrav=1,2,3,6,8`. Any other Bravais lattice returns a stable
//! `invalid_input` error that asks for `ibrav=0`.

use xtalforge_core::{Error, Structure};

mod input;
mod output;
mod write;

/// Bohr radius in Angstrom.
pub(crate) const BOHR: f64 = 0.529_177_210_903;

/// Rydberg in electronvolt.
pub(crate) const RYDBERG_EV: f64 = 13.605_693_122_994;

/// Reads a `pw.x` input deck or output file.
///
/// # Errors
///
/// Returns [`Error::InvalidInput`] for a missing cell, a missing atom card, an
/// unsupported Bravais lattice or a malformed number.
pub fn read_pwscf(text: &str) -> Result<Structure, Error> {
    if looks_like_output(text) {
        output::read(text)
    } else {
        input::read(text)
    }
}

/// Writes a P1 structure as a `pw.x` input deck.
///
/// # Errors
///
/// Returns [`Error::InvalidInput`] for an empty structure or a failed P1
/// expansion.
pub fn write_pwscf(structure: &Structure) -> Result<String, Error> {
    write::write(structure)
}

fn looks_like_output(text: &str) -> bool {
    [
        "Program PWSCF",
        "Self-consistent Calculation",
        "End of self-consistent calculation",
        "number of atoms/cell",
        "crystal axes:",
        "site n.     atom",
    ]
    .iter()
    .any(|marker| text.contains(marker))
}

/// Strips `!` and `#` comments outside quotes and keeps line breaks.
pub(crate) fn strip_comments(text: &str) -> String {
    let mut out = String::with_capacity(text.len());
    let mut quote: Option<char> = None;
    let mut comment = false;
    for character in text.chars() {
        if character == '\n' {
            out.push('\n');
            comment = false;
            continue;
        }
        if comment {
            continue;
        }
        if let Some(active) = quote {
            out.push(character);
            if character == active {
                quote = None;
            }
            continue;
        }
        match character {
            '\'' | '"' => {
                quote = Some(character);
                out.push(character);
            }
            '!' | '#' => comment = true,
            other => out.push(other),
        }
    }
    out
}

/// Finds a namelist variable value (case-insensitive, first `=` on a line).
pub(crate) fn variable(text: &str, name: &str) -> Option<String> {
    for line in strip_comments(text).lines() {
        let Some((left, right)) = line.split_once('=') else {
            continue;
        };
        if left.trim().eq_ignore_ascii_case(name) {
            let value = right
                .split([',', ' ', '\t'])
                .find(|token| !token.is_empty())
                .unwrap_or("")
                .trim()
                .trim_matches(['\'', '"']);
            if !value.is_empty() {
                return Some(value.to_owned());
            }
        }
    }
    None
}

/// Reads the three numbers of a lattice or position line.
pub(crate) fn three_numbers(line: &str) -> Option<[f64; 3]> {
    let mut values = line
        .split_whitespace()
        .filter_map(|value| value.parse::<f64>().ok());
    Some([values.next()?, values.next()?, values.next()?])
}

/// Reads the parenthesised vector of a `pw.x` output line.
pub(crate) fn parenthesised(line: &str) -> Option<[f64; 3]> {
    let equals = line.find('=')?;
    let rest = line.get(equals..)?;
    let start = rest.find('(')? + equals;
    let end = line.get(start..)?.find(')')? + start;
    three_numbers(line.get(start + 1..end)?)
}

/// Scale factor from card units to Angstrom.
pub(crate) fn unit_scale(unit: &str, alat: Option<f64>) -> Result<f64, Error> {
    match unit {
        "" | "alat" => alat.ok_or_else(|| {
            Error::InvalidInput(
                "a pw.x card in alat units needs celldm(1) or A in &system".to_owned(),
            )
        }),
        "angstrom" | "ang" => Ok(1.0),
        "bohr" | "a.u." | "au" => Ok(BOHR),
        other => Err(Error::InvalidInput(format!(
            "unsupported pw.x card unit {other:?}; use angstrom, bohr or alat"
        ))),
    }
}

/// Reads the unit token of a card header such as `ATOMIC_POSITIONS {crystal}`.
pub(crate) fn card_unit(header: &str) -> String {
    let rest = header
        .split_once(char::is_whitespace)
        .map_or("", |(_, rest)| rest)
        .trim();
    rest.trim_matches(|character| matches!(character, '{' | '}' | '(' | ')'))
        .split_whitespace()
        .next()
        .unwrap_or("")
        .to_ascii_lowercase()
}
