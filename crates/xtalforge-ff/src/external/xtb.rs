//! xTB coordinate input and output parsing (Task 4.9.3).

use serde::{Deserialize, Serialize};
use xtalforge_core::{Structure, element_info};

use super::text::{find_line, find_number, finite, push_line, validate_token};
use super::{EngineResult, HARTREE_PER_BOHR_TO_EV_PER_ANGSTROM, HARTREE_TO_EV};
use crate::Error;

/// xTB input options.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct XtbOptions {
    /// Method name recorded in the XYZ comment, for example `GFN2-xTB`.
    pub method: String,
    /// Total charge of the system in units of the elementary charge.
    pub charge: f64,
    /// Number of unpaired electrons.
    pub uhf: u32,
}

impl Default for XtbOptions {
    fn default() -> Self {
        Self {
            method: "GFN2-xTB".to_owned(),
            charge: 0.0,
            uhf: 0,
        }
    }
}

impl XtbOptions {
    /// Checks the xTB options.
    ///
    /// # Errors
    ///
    /// Returns [`Error::InvalidInput`] for an unsafe method name or a
    /// non-finite charge.
    pub fn validate(&self) -> Result<(), Error> {
        validate_token("xTB method", &self.method)?;
        if !self.charge.is_finite() {
            return Err(Error::InvalidInput("xTB charge must be finite".to_owned()));
        }
        Ok(())
    }
}

/// Writes an xTB coordinate file for `structure`.
///
/// The file starts with an XYZ block (count, comment, atoms) and continues
/// with `$coord`, `$chrg` and `$uhf` blocks.  xTB selects the method on its
/// command line, so `options.method` is recorded in the XYZ comment.
///
/// # Errors
///
/// Returns [`Error::InvalidInput`] for invalid options, an empty structure or
/// a dummy atom, and a core error when the P1 expansion or an element lookup
/// fails.
pub fn write_xtb_input(structure: &Structure, options: &XtbOptions) -> Result<String, Error> {
    options.validate()?;
    let cell = structure.p1()?;
    if cell.n_atoms() == 0 {
        return Err(Error::InvalidInput(
            "cannot write an xTB input for an empty structure".to_owned(),
        ));
    }
    let mut out = String::new();
    push_line(&mut out, cell.n_atoms().to_string());
    push_line(
        &mut out,
        format!("xtalforge {} method={}", structure.name(), options.method),
    );
    let positions = cell.lattice.to_cart_many(&cell.frac());
    for (atom, position) in cell.atoms.iter().zip(&positions) {
        if element_info(&atom.element)?.is_dummy() {
            return Err(Error::InvalidInput(
                "external engines do not accept the dummy X species".to_owned(),
            ));
        }
        push_line(
            &mut out,
            format!(
                "{:<2} {:.10} {:.10} {:.10}",
                atom.element, position[0], position[1], position[2]
            ),
        );
    }
    push_line(&mut out, "$coord");
    for (atom, position) in cell.atoms.iter().zip(&positions) {
        push_line(
            &mut out,
            format!(
                "{:.10} {:.10} {:.10} {:<2}",
                position[0], position[1], position[2], atom.element
            ),
        );
    }
    push_line(&mut out, "$end");
    push_line(&mut out, format!("$chrg {:.6}", options.charge));
    push_line(&mut out, format!("$uhf {}", options.uhf));
    Ok(out)
}

/// Parses an xTB output file.
///
/// `TOTAL ENERGY` is in Hartree and converted to eV.  A `TOTAL CHARGE` line is
/// validated when present.  Force rows printed under a line containing
/// `forces` are in Hartree/Bohr and converted to eV/Angstrom; no such block
/// yields an empty vector.
///
/// # Errors
///
/// Returns [`Error::InvalidInput`] when `TOTAL ENERGY` is missing or
/// malformed, or when `TOTAL CHARGE` is present but malformed.
pub fn parse_xtb_output(text: &str) -> Result<EngineResult, Error> {
    let lines: Vec<&str> = text.lines().collect();
    let energy = find_number(&lines, "TOTAL ENERGY")
        .ok_or_else(|| Error::InvalidInput("xTB output has no 'TOTAL ENERGY' line".to_owned()))?
        * HARTREE_TO_EV;
    finite("xTB energy", energy)?;
    if find_line(&lines, "TOTAL CHARGE").is_some() {
        let charge = find_number(&lines, "TOTAL CHARGE").ok_or_else(|| {
            Error::InvalidInput("xTB 'TOTAL CHARGE' line has no number".to_owned())
        })?;
        finite("xTB total charge", charge)?;
    }
    let forces = lines
        .iter()
        .position(|line| line.to_ascii_lowercase().contains("forces"))
        .map(|index| {
            super::text::parse_vectors(&lines, index + 1, HARTREE_PER_BOHR_TO_EV_PER_ANGSTROM)
        })
        .unwrap_or_default();
    Ok(EngineResult {
        energy,
        forces,
        stress: None,
    })
}
