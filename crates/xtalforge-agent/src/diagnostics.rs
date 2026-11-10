//! Inspection and diagnostics (Task 8.1.3).
//!
//! `diagnose` returns what an assistant needs in one call: overlapping atoms,
//! missing hydrogens, unusual coordination, partial occupancy, charge
//! imbalance and symmetry problems, each with a severity, the atoms involved
//! and the verb that would fix it.
//!
//! Ported from Crystal-Builder `xtal/agent/diagnostics.py` (MIT), with the
//! checks expressed against the typed structures of this library.

use std::collections::BTreeMap;

use serde::{Deserialize, Serialize};
use xtalforge_core::elements::element_by_symbol;
use xtalforge_core::{Element, Structure};

use crate::Error;

/// How serious a finding is.
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Serialize, Deserialize)]
#[serde(rename_all = "snake_case")]
pub enum Severity {
    /// Worth knowing; the structure is usable.
    Info,
    /// Likely a mistake in the model or the input file.
    Warning,
    /// The structure cannot be interpreted as it stands.
    Error,
}

impl Severity {
    /// Lower-case name used in JSON and text output.
    #[must_use]
    pub fn as_str(self) -> &'static str {
        match self {
            Self::Info => "info",
            Self::Warning => "warning",
            Self::Error => "error",
        }
    }
}

/// One finding.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct Diagnostic {
    /// Stable `snake_case` code.
    pub code: String,
    /// How serious it is.
    pub severity: Severity,
    /// One-line description.
    pub message: String,
    /// Atom indices in the P1 cell, when the finding names atoms.
    pub atoms: Vec<usize>,
    /// Verb that would fix it, when there is one.
    pub suggestion: Option<String>,
}

/// Everything `diagnose` found.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct DiagnosticReport {
    /// Atoms in the P1 cell.
    pub n_atoms: usize,
    /// Findings, most severe first.
    pub diagnostics: Vec<Diagnostic>,
    /// Number of findings per severity.
    pub counts: BTreeMap<String, usize>,
}

impl DiagnosticReport {
    /// Whether any finding is an error.
    #[must_use]
    pub fn has_errors(&self) -> bool {
        self.diagnostics
            .iter()
            .any(|diagnostic| diagnostic.severity == Severity::Error)
    }

    /// Finds one code.
    #[must_use]
    pub fn find(&self, code: &str) -> Option<&Diagnostic> {
        self.diagnostics
            .iter()
            .find(|diagnostic| diagnostic.code == code)
    }
}

fn covalent_radius(element: &str) -> f64 {
    element_by_symbol(element).map_or(0.8, |info| info.covalent_radius)
}

/// Typical maximum coordination for the light elements the check covers.
fn typical_valence(element: &Element) -> Option<usize> {
    match element.symbol.as_str() {
        "H" | "D" => Some(1),
        "C" => Some(4),
        "N" => Some(4),
        "O" => Some(2),
        "F" | "Cl" | "Br" | "I" => Some(1),
        "S" => Some(6),
        "P" => Some(5),
        "B" => Some(4),
        _ => None,
    }
}

fn push(
    diagnostics: &mut Vec<Diagnostic>,
    code: &'static str,
    severity: Severity,
    message: String,
    atoms: Vec<usize>,
    suggestion: Option<&'static str>,
) {
    diagnostics.push(Diagnostic {
        code: code.to_owned(),
        severity,
        message,
        atoms,
        suggestion: suggestion.map(str::to_owned),
    });
}

/// Inspects a structure and reports what is wrong with it.
///
/// # Errors
///
/// Returns the P1 expansion or bond perception error.
pub fn diagnose(structure: &Structure) -> Result<DiagnosticReport, Error> {
    let p1 = structure.p1()?;
    let graph = structure.bond_graph()?;
    let mut diagnostics = Vec::new();
    if p1.atoms.is_empty() {
        push(
            &mut diagnostics,
            "empty_structure",
            Severity::Error,
            "the structure has no atoms".to_owned(),
            Vec::new(),
            None,
        );
        return Ok(report(p1.atoms.len(), diagnostics));
    }

    // Overlapping atoms: the shortest periodic contact, excluding a site
    // with itself. Anything below 0.4 A is physically impossible.
    let mut overlaps = Vec::new();
    let mut shortest = f64::INFINITY;
    for (i, first) in p1.atoms.iter().enumerate() {
        for (j, second) in p1.atoms.iter().enumerate().skip(i + 1) {
            let distance = p1.lattice.distance(first.frac, second.frac);
            if distance < shortest {
                shortest = distance;
            }
            if distance < 0.4 {
                overlaps.push((i, j));
            }
        }
    }
    if !overlaps.is_empty() {
        let mut atoms: Vec<usize> = overlaps.iter().flat_map(|(i, j)| [*i, *j]).collect();
        atoms.sort_unstable();
        atoms.dedup();
        push(
            &mut diagnostics,
            "overlapping_atoms",
            Severity::Error,
            format!(
                "{} pair(s) of atoms closer than 0.4 A; shortest contact {:.3} A",
                overlaps.len(),
                shortest
            ),
            atoms,
            None,
        );
    }

    // Partial occupancy: a model choice that changes the density and the
    // charge, so it is reported even when it is intentional.
    let partial: Vec<usize> = p1
        .atoms
        .iter()
        .enumerate()
        .filter(|(_, atom)| atom.occupancy < 0.999)
        .map(|(index, _)| index)
        .collect();
    if !partial.is_empty() {
        push(
            &mut diagnostics,
            "partial_occupancy",
            Severity::Info,
            format!(
                "{} site(s) have occupancy below 1; the formula counts them differently",
                partial.len()
            ),
            partial,
            None,
        );
    }

    // Charge: the structure's stored total charge.
    let total_charge = structure.summary()?.total_charge;
    if total_charge.abs() > 1e-6 {
        push(
            &mut diagnostics,
            "charge_imbalance",
            Severity::Warning,
            format!("the stored total charge is {total_charge:+.3} e, not zero"),
            Vec::new(),
            Some("info"),
        );
    }

    // Coordination: a dangling light atom, or a light atom with fewer bonds
    // than its valence, usually means the deposited model has no hydrogens.
    let mut undercoordinated = Vec::new();
    let mut dangling = Vec::new();
    let mut overcoordinated = Vec::new();
    for (index, atom) in p1.atoms.iter().enumerate() {
        let Ok(element) = element_by_symbol(&atom.element) else {
            continue;
        };
        let bonds = graph.degree(index);
        if bonds == 0 {
            if !element.is_metal && !element.is_dummy() {
                dangling.push(index);
            }
            continue;
        }
        if let Some(valence) = typical_valence(element)
            && bonds < valence
            && matches!(element.symbol.as_str(), "C" | "N" | "O")
        {
            undercoordinated.push(index);
        }
        if bonds > 12 {
            overcoordinated.push(index);
        }
    }
    if !dangling.is_empty() {
        push(
            &mut diagnostics,
            "dangling_atoms",
            Severity::Warning,
            format!(
                "{} non-metal atom(s) have no bonds; check the bond rules or the site labels",
                dangling.len()
            ),
            dangling,
            Some("bonds"),
        );
    }
    if !undercoordinated.is_empty() {
        push(
            &mut diagnostics,
            "missing_hydrogens",
            Severity::Info,
            format!(
                "{} C/N/O atom(s) have fewer bonds than their valence; the model probably lacks hydrogens",
                undercoordinated.len()
            ),
            undercoordinated,
            Some("hydrogens"),
        );
    }
    if !overcoordinated.is_empty() {
        push(
            &mut diagnostics,
            "unusual_coordination",
            Severity::Warning,
            format!(
                "{} atom(s) have more than 12 bonds; the bond rules are probably too generous",
                overcoordinated.len()
            ),
            overcoordinated,
            Some("bonds"),
        );
    }

    // Symmetry: a failure is reported as info because a deliberately
    // asymmetric model is valid.
    match xtalforge_symmetry::detect_symmetry(structure, Default::default()) {
        Ok(dataset) => {
            if dataset.number == 1 && p1.atoms.len() > 3 {
                push(
                    &mut diagnostics,
                    "no_symmetry",
                    Severity::Info,
                    "the structure is P1; a CIF that declared a space group may have lost it"
                        .to_owned(),
                    Vec::new(),
                    Some("symmetry"),
                );
            }
        }
        Err(error) => push(
            &mut diagnostics,
            "symmetry_failed",
            Severity::Info,
            format!("symmetry detection did not finish: {error}"),
            Vec::new(),
            None,
        ),
    }

    diagnostics.sort_by(|left, right| {
        right
            .severity
            .cmp(&left.severity)
            .then_with(|| left.code.cmp(&right.code))
    });
    Ok(report(p1.atoms.len(), diagnostics))
}

fn report(n_atoms: usize, diagnostics: Vec<Diagnostic>) -> DiagnosticReport {
    let mut counts = BTreeMap::new();
    for diagnostic in &diagnostics {
        *counts
            .entry(diagnostic.severity.as_str().to_owned())
            .or_insert(0) += 1;
    }
    for severity in ["info", "warning", "error"] {
        counts.entry(severity.to_owned()).or_insert(0);
    }
    DiagnosticReport {
        n_atoms,
        diagnostics,
        counts,
    }
}

/// Shortest periodic contact in a structure, in Angstrom.
///
/// # Errors
///
/// Returns the P1 expansion error.
pub fn shortest_contact(structure: &Structure) -> Result<f64, Error> {
    let p1 = structure.p1()?;
    let mut shortest = f64::INFINITY;
    for (i, first) in p1.atoms.iter().enumerate() {
        for second in p1.atoms.iter().skip(i + 1) {
            shortest = shortest.min(p1.lattice.distance(first.frac, second.frac));
        }
    }
    if shortest.is_finite() {
        Ok(shortest)
    } else {
        Err(Error::Core(xtalforge_core::Error::InvalidInput(
            "a contact needs at least two atoms".to_owned(),
        )))
    }
}

/// Mean covalent-radius sum used by the bond rules, kept for callers that
/// want to compare a contact with the perception threshold.
///
/// # Errors
///
/// Returns [`xtalforge_core::Error::InvalidInput`] for an unknown element.
pub fn bond_threshold(first: &str, second: &str, scale: f64) -> Result<f64, Error> {
    Ok((covalent_radius(first) + covalent_radius(second)) * scale)
}
