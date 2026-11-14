//! Helpers shared by the engine parsers (Task 18.2).
//!
//! Units at this boundary are Angstrom, eV, eV/Angstrom and GPa. Each parser
//! converts the engine's native units once, here, with the constants below.

use xtalforge_core::elements::parse_symbol;
use xtalforge_core::{Lattice, Site, Structure, StructureBuilder};

use crate::Error;
use crate::result::CalculationResult;

/// Hartree in electronvolt.
pub(crate) const HARTREE_EV: f64 = 27.211_386_245_988;
/// Rydberg in electronvolt.
pub(crate) const RYDBERG_EV: f64 = 13.605_693_122_994;
/// Bohr radius in Angstrom.
pub(crate) const BOHR_ANGSTROM: f64 = 0.529_177_210_903;
/// Hartree/bohr in eV/Angstrom.
pub(crate) const HARTREE_BOHR_EV_ANGSTROM: f64 = 51.422_067_476_3;
/// Rydberg/bohr in eV/Angstrom.
pub(crate) const RYDBERG_BOHR_EV_ANGSTROM: f64 = 25.711_043_095_4;
/// Rydberg/bohr^3 in GPa.
pub(crate) const RYDBERG_BOHR3_GPA: f64 = 14_710.0;
/// One kBar in GPa.
pub(crate) const KBAR_GPA: f64 = 0.1;
/// One bar in GPa (LAMMPS metal units).
pub(crate) const BAR_GPA: f64 = 0.000_1;
/// One atmosphere in GPa (LAMMPS real units).
pub(crate) const ATM_GPA: f64 = 0.000_101_325;
/// One kcal/mol in eV (LAMMPS real units).
pub(crate) const KCAL_MOL_EV: f64 = 0.043_364_115;

/// Parses one number, accepting Fortran `D` exponents and a trailing comma.
pub(crate) fn number(text: &str) -> Option<f64> {
    let cleaned: String = text
        .trim()
        .trim_matches(',')
        .chars()
        .map(|character| match character {
            'D' | 'd' => 'E',
            other => other,
        })
        .collect();
    if cleaned.is_empty() {
        return None;
    }
    cleaned.parse().ok()
}

/// Parses every number on a line, in order.
pub(crate) fn numbers(line: &str) -> Vec<f64> {
    line.split_whitespace().filter_map(number).collect()
}

/// Reads the first number after `marker`.
pub(crate) fn value_after(line: &str, marker: &str) -> Option<f64> {
    let (_, rest) = line.split_once(marker)?;
    numbers(rest).first().copied()
}

/// Parses `true`/`false`/`T`/`F`/`1`/`0` flags used by XML outputs.
pub(crate) fn flag(text: &str) -> Option<bool> {
    match text.trim().to_ascii_lowercase().as_str() {
        "true" | "t" | ".true." | "1" | "yes" => Some(true),
        "false" | "f" | ".false." | "0" | "no" => Some(false),
        _ => None,
    }
}

/// Builds a structure from a cell matrix, element symbols and fractional
/// coordinates, with one-based site labels.
pub(crate) fn build_structure(
    name: &str,
    matrix: [[f64; 3]; 3],
    elements: &[String],
    frac: &[[f64; 3]],
) -> Result<Structure, Error> {
    if elements.len() != frac.len() {
        return Err(Error::Parse(format!(
            "{name} has {} element symbols for {} positions",
            elements.len(),
            frac.len()
        )));
    }
    if elements.is_empty() {
        return Err(Error::Parse(format!("{name} has no atoms")));
    }
    let lattice = Lattice::new(matrix).map_err(|error| Error::Parse(format!("{name}: {error}")))?;
    let mut builder = StructureBuilder::new(name, lattice);
    for (index, (element, point)) in elements.iter().zip(frac).enumerate() {
        if !point.iter().all(|value| value.is_finite()) {
            return Err(Error::Parse(format!(
                "{name} site {} has a non-finite coordinate",
                index + 1
            )));
        }
        let symbol =
            parse_symbol(element).map_err(|error| Error::Parse(format!("{name}: {error}")))?;
        builder = builder.site(Site::new(
            format!("{symbol}{}", index + 1),
            &symbol,
            *point,
        )?);
    }
    builder
        .build()
        .map_err(|error| Error::Parse(format!("{name}: {error}")))
}

/// Largest force norm in a per-site force list.
pub(crate) fn max_force(forces: &[[f64; 3]]) -> f64 {
    forces
        .iter()
        .map(|force| force.iter().map(|value| value * value).sum::<f64>().sqrt())
        .fold(0.0, f64::max)
}

/// Converts a three-by-three stress block, rejecting incomplete rows.
pub(crate) fn stress_from_rows(rows: &[Vec<f64>], factor: f64) -> Option<[[f64; 3]; 3]> {
    let mut matrix = [[0.0; 3]; 3];
    for (row_index, row) in rows.iter().take(3).enumerate() {
        let target = matrix.get_mut(row_index)?;
        for (column, slot) in target.iter_mut().enumerate() {
            *slot = *row.get(column)? * factor;
        }
    }
    (rows.len() >= 3).then_some(matrix)
}

/// Appends a stable warning entry.
pub(crate) fn warn(result: &mut CalculationResult, code: &str, message: &str) {
    result.push_warning(code, message);
}

/// Splits a line into whitespace-separated fields.
pub(crate) fn fields(line: &str) -> Vec<&str> {
    line.split_whitespace().collect()
}

/// Finds the last line matching a predicate.
pub(crate) fn rposition<F>(lines: &[&str], predicate: F) -> Option<usize>
where
    F: Fn(&str) -> bool,
{
    lines.iter().rposition(|line| predicate(line))
}
