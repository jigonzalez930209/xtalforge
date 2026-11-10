//! Structure inspection (Task 8.1.3).
//!
//! One call that answers "what is this?": the summary, the cell, the
//! composition, the site table and the bonding census. The problems are in
//! [`crate::diagnostics`]; the session's `info` verb returns both.

use std::collections::BTreeMap;

use serde::{Deserialize, Serialize};
use xtalforge_core::elements::element_by_symbol;
use xtalforge_core::{Structure, StructureSummary};

use crate::Error;

/// Cell parameters in the usual units.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct CellSummary {
    /// `a`, `b`, `c` in Angstrom.
    pub lengths: [f64; 3],
    /// `alpha`, `beta`, `gamma` in degrees.
    pub angles: [f64; 3],
    /// Volume in cubic Angstrom.
    pub volume: f64,
}

/// One row of the site table.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct SiteRow {
    /// Index in the asymmetric unit.
    pub index: usize,
    /// Site label.
    pub label: String,
    /// Element symbol.
    pub element: String,
    /// Fractional coordinates.
    pub frac: [f64; 3],
    /// Occupancy.
    pub occupancy: f64,
    /// Stored charge, when present.
    pub charge: Option<f64>,
    /// Bonds perceived for the first P1 image of this site.
    pub bonds: usize,
}

/// One element of the composition.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct ElementCount {
    /// Element symbol.
    pub element: String,
    /// Occupancy-weighted atoms in the P1 cell.
    pub count: f64,
    /// Mass fraction in the cell.
    pub mass_fraction: f64,
}

/// Everything `inspect` reports.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct Inspection {
    /// Name, formula, density and charge.
    pub summary: StructureSummaryView,
    /// Cell parameters.
    pub cell: CellSummary,
    /// Atoms in the P1 cell.
    pub n_p1_atoms: usize,
    /// Perceived bonds in the P1 cell.
    pub n_bonds: usize,
    /// Composition, heaviest first.
    pub composition: Vec<ElementCount>,
    /// Site table.
    pub sites: Vec<SiteRow>,
}

/// Serialisable mirror of the core summary (the core type is `Serialize`
/// only, and an answer must round-trip through JSON).
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct StructureSummaryView {
    /// Structure name.
    pub name: String,
    /// Formula.
    pub formula: String,
    /// Hall symbol.
    pub hall: String,
    /// Number of asymmetric sites.
    pub n_sites: usize,
    /// Volume in cubic Angstrom.
    pub volume: f64,
    /// Density in g/cm3.
    pub density: f64,
    /// Total charge.
    pub total_charge: f64,
    /// Whether the structure is P1.
    pub is_p1: bool,
}

impl From<StructureSummary> for StructureSummaryView {
    fn from(summary: StructureSummary) -> Self {
        Self {
            name: summary.name,
            formula: summary.formula,
            hall: summary.hall,
            n_sites: summary.n_sites,
            volume: summary.volume,
            density: summary.density,
            total_charge: summary.total_charge,
            is_p1: summary.is_p1,
        }
    }
}

/// Inspects a structure.
///
/// # Errors
///
/// Returns the P1 expansion or bond perception error.
pub fn inspect_structure(structure: &Structure) -> Result<Inspection, Error> {
    let summary = StructureSummaryView::from(structure.summary()?);
    let p1 = structure.p1()?;
    let graph = structure.bond_graph()?;
    let parameters = structure.lattice().parameters();

    let mut counts: BTreeMap<String, f64> = BTreeMap::new();
    let mut masses: BTreeMap<String, f64> = BTreeMap::new();
    let mut total_mass = 0.0;
    for atom in &p1.atoms {
        let mass = element_by_symbol(&atom.element).map_or(0.0, |info| info.mass);
        *counts.entry(atom.element.clone()).or_insert(0.0) += atom.occupancy;
        *masses.entry(atom.element.clone()).or_insert(0.0) += mass * atom.occupancy;
        total_mass += mass * atom.occupancy;
    }
    let mut composition: Vec<ElementCount> = counts
        .into_iter()
        .map(|(element, count)| ElementCount {
            mass_fraction: if total_mass > 0.0 {
                masses.get(&element).copied().unwrap_or(0.0) / total_mass
            } else {
                0.0
            },
            element,
            count,
        })
        .collect();
    composition.sort_by(|left, right| {
        right
            .mass_fraction
            .total_cmp(&left.mass_fraction)
            .then_with(|| left.element.cmp(&right.element))
    });

    let mut sites = Vec::with_capacity(structure.sites().len());
    for (index, site) in structure.sites().iter().enumerate() {
        let first_image = p1
            .atoms
            .iter()
            .position(|atom| atom.site_idx == index)
            .map_or(0, |atom| graph.degree(atom));
        sites.push(SiteRow {
            index,
            label: site.label().to_owned(),
            element: site.element().to_owned(),
            frac: site.frac(),
            occupancy: site.occupancy(),
            charge: site.charge(),
            bonds: first_image,
        });
    }

    Ok(Inspection {
        summary,
        cell: CellSummary {
            lengths: [parameters[0], parameters[1], parameters[2]],
            angles: [parameters[3], parameters[4], parameters[5]],
            volume: structure.lattice().volume(),
        },
        n_p1_atoms: p1.atoms.len(),
        n_bonds: graph.bonds.len(),
        composition,
        sites,
    })
}
