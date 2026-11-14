//! Substitution-based generation (Task 20.3.3).
//!
//! Given a structure (typically a prototype), the oxidation state of each
//! element is estimated with the bond-valence machinery of Task 9.4 and every
//! ion in the embedded Shannon table with the same charge and a radius close
//! to the original is proposed. The candidate list is deterministic: elements
//! keep their order and candidates are sorted by radius mismatch, then by
//! symbol.
//!
//! The result is a composition proposal, not a relaxed structure. Apply a
//! candidate with [`apply_substitution`] and relax it with the force field if
//! needed.

use std::collections::BTreeMap;

use serde::{Deserialize, Serialize};
use xtalforge_analysis::{BvsOptions, oxidation_states};
use xtalforge_core::{Error, Structure, StructureBuilder};

/// One ionic radius.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct IonicRadius {
    /// Element symbol.
    pub element: String,
    /// Formal charge.
    pub charge: i32,
    /// Effective ionic radius in Angstrom (coordination VI).
    pub radius: f64,
}

const RADII_JSON: &str = include_str!("../data/ionic_radii.json");

#[derive(Debug, Deserialize)]
struct RadiusFile {
    ions: Vec<IonicRadius>,
}

/// The embedded Shannon ionic-radius table.
///
/// # Errors
///
/// Returns [`Error::DataJson`] when the file cannot be parsed.
pub fn ionic_radii() -> Result<Vec<IonicRadius>, Error> {
    let file: RadiusFile = serde_json::from_str(RADII_JSON).map_err(|error| Error::DataJson {
        name: "xtalforge-generation/ionic_radii.json".to_owned(),
        reason: error.to_string(),
    })?;
    Ok(file.ions)
}

/// Options for [`substitution_candidates`].
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct SubstitutionOptions {
    /// Maximum relative radius difference, for example `0.15` for 15 %.
    pub radius_tolerance: f64,
    /// Whether candidates must keep the estimated oxidation state.
    pub same_charge: bool,
    /// Maximum number of candidates.
    pub limit: usize,
}

impl Default for SubstitutionOptions {
    fn default() -> Self {
        Self {
            radius_tolerance: 0.15,
            same_charge: true,
            limit: 32,
        }
    }
}

/// One proposed substitution.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct SubstitutionCandidate {
    /// Element replaced.
    pub from: String,
    /// Oxidation state used for the rule.
    pub from_charge: i32,
    /// Proposed element.
    pub to: String,
    /// Charge of the proposed ion.
    pub to_charge: i32,
    /// Radius of the original ion.
    pub from_radius: f64,
    /// Radius of the proposed ion.
    pub to_radius: f64,
    /// Relative radius difference.
    pub mismatch: f64,
}

/// Proposes chemically similar substitutions for every element.
///
/// # Errors
///
/// Returns the errors of the bond-valence estimator and of the embedded table.
pub fn substitution_candidates(
    structure: &Structure,
    options: &SubstitutionOptions,
) -> Result<Vec<SubstitutionCandidate>, Error> {
    if !options.radius_tolerance.is_finite() || options.radius_tolerance < 0.0 {
        return Err(Error::InvalidInput(
            "radius tolerance must be finite and non-negative".to_owned(),
        ));
    }
    if options.limit == 0 {
        return Err(Error::InvalidInput(
            "the candidate limit must be at least one".to_owned(),
        ));
    }
    let report = oxidation_states(structure, &BvsOptions::default())?;
    let radii = ionic_radii()?;
    let p1 = structure.p1()?;
    let mut elements: Vec<String> = Vec::new();
    for atom in &p1.atoms {
        if !elements.contains(&atom.element) {
            elements.push(atom.element.clone());
        }
    }
    let mut candidates = Vec::new();
    for element in &elements {
        let charge = report.by_element.get(element).copied().unwrap_or(0);
        let Some(original) = radii
            .iter()
            .find(|ion| ion.element == *element && ion.charge == charge)
            .or_else(|| radii.iter().find(|ion| ion.element == *element))
        else {
            continue;
        };
        for ion in &radii {
            if ion.element == *element {
                continue;
            }
            if options.same_charge && charge != 0 && ion.charge != charge {
                continue;
            }
            let mismatch = (ion.radius - original.radius).abs() / original.radius.max(1e-6);
            if mismatch > options.radius_tolerance {
                continue;
            }
            candidates.push(SubstitutionCandidate {
                from: element.clone(),
                from_charge: charge,
                to: ion.element.clone(),
                to_charge: ion.charge,
                from_radius: original.radius,
                to_radius: ion.radius,
                mismatch,
            });
        }
    }
    candidates.sort_by(|left, right| {
        left.mismatch
            .total_cmp(&right.mismatch)
            .then_with(|| left.from.cmp(&right.from))
            .then_with(|| left.to.cmp(&right.to))
            .then_with(|| left.to_charge.cmp(&right.to_charge))
    });
    candidates.truncate(options.limit);
    Ok(candidates)
}

/// Applies one substitution to every site of an element.
///
/// Bonds and metadata are preserved; labels keep their suffix so the mapping
/// back to the original structure stays visible.
///
/// # Errors
///
/// Returns [`Error::InvalidInput`] for an unknown target element or a builder
/// error.
pub fn apply_substitution(structure: &Structure, from: &str, to: &str) -> Result<Structure, Error> {
    xtalforge_core::elements::element_by_symbol(to)?;
    let mut builder =
        StructureBuilder::new(format!("{}-{}", structure.name(), to), *structure.lattice())
            .hall(structure.hall())
            .metadata("generator", serde_json::json!("substitution"))
            .metadata("substituted_from", serde_json::json!(from))
            .metadata("substituted_to", serde_json::json!(to));
    for operation in structure.operations() {
        builder = builder.operation(*operation);
    }
    for (key, value) in structure.metadata() {
        builder = builder.metadata(key.clone(), value.clone());
    }
    for bond in structure.explicit_bonds() {
        builder = builder.explicit_bond(*bond);
    }
    for bond in structure.topology_bonds() {
        builder = builder.topology_bond(*bond);
    }
    for site in structure.sites() {
        if site.element() == from {
            builder = builder.site(
                xtalforge_core::Site::new(site.label(), to, site.frac())?
                    .with_occupancy(site.occupancy())?,
            );
        } else {
            builder = builder.site(site.clone());
        }
    }
    builder.build()
}

/// A JSON-ready substitution report.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct SubstitutionReport {
    /// Candidates in ranked order.
    pub candidates: Vec<SubstitutionCandidate>,
    /// Distinct proposed compositions, for example `Na->K`.
    pub compositions: Vec<String>,
    /// Whether the oxidation-state estimate was charge neutral.
    pub neutral: bool,
}

impl SubstitutionReport {
    /// Serializes the report to JSON.
    ///
    /// # Errors
    ///
    /// Returns [`Error::InvalidInput`] when serialization fails.
    pub fn to_json(&self) -> Result<String, Error> {
        serde_json::to_string(self).map_err(|error| {
            Error::InvalidInput(format!("could not serialize substitution report: {error}"))
        })
    }
}

/// Runs [`substitution_candidates`] and groups the result by composition.
///
/// # Errors
///
/// Returns the errors of [`substitution_candidates`] and the bond-valence
/// estimator.
pub fn substitution_report(
    structure: &Structure,
    options: &SubstitutionOptions,
) -> Result<SubstitutionReport, Error> {
    let candidates = substitution_candidates(structure, options)?;
    let report = oxidation_states(structure, &BvsOptions::default())?;
    let mut compositions: Vec<String> = candidates
        .iter()
        .map(|candidate| format!("{}->{}", candidate.from, candidate.to))
        .collect();
    compositions.dedup();
    Ok(SubstitutionReport {
        candidates,
        compositions,
        neutral: report.neutral,
    })
}

/// Applies every candidate that changes one distinct element.
///
/// # Errors
///
/// Returns the errors of [`substitution_candidates`] and
/// [`apply_substitution`].
pub fn substituted_structures(
    structure: &Structure,
    options: &SubstitutionOptions,
) -> Result<Vec<(SubstitutionCandidate, Structure)>, Error> {
    let candidates = substitution_candidates(structure, options)?;
    let mut out = Vec::new();
    for candidate in candidates {
        let substituted = apply_substitution(structure, &candidate.from, &candidate.to)?;
        out.push((candidate, substituted));
    }
    Ok(out)
}

/// Groups candidates by `from` element.
#[must_use]
pub fn group_by_element(
    candidates: &[SubstitutionCandidate],
) -> BTreeMap<String, Vec<&SubstitutionCandidate>> {
    let mut grouped: BTreeMap<String, Vec<&SubstitutionCandidate>> = BTreeMap::new();
    for candidate in candidates {
        grouped
            .entry(candidate.from.clone())
            .or_default()
            .push(candidate);
    }
    grouped
}
