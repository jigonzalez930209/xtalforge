//! Prototype library and identification (Task 20.3).
//!
//! The data file `data/prototypes.json` is a curated subset of the
//! Strukturbericht/AFLOW prototype labels. Each entry records the space group,
//! the occupied Wyckoff letters, the ideal free parameters and the lattice
//! ratios; it contains no copied AFLOW records, only the published
//! crystallographic facts. [`THIRD_PARTY.md`](../../THIRD_PARTY.md) records
//! the source.
//!
//! [`build_prototype`] builds any entry (optionally scaled), and
//! [`identify_prototype`] assigns the label of an input structure from its
//! space group, Wyckoff letters and stoichiometry and confirms the assignment
//! with the structure matcher.

use std::collections::BTreeMap;

use serde::{Deserialize, Serialize};
use xtalforge_core::{Error, Lattice, Site, Structure, StructureBuilder};
use xtalforge_symmetry::{
    MatcherOptions, SpaceGroup, SymmetryOptions, WyckoffPosition, detect_symmetry,
    match_structures, wyckoff_assignments,
};

use crate::geometry::evaluate_template;

const PROTOTYPES_JSON: &str = include_str!("../data/prototypes.json");

#[derive(Debug, Clone, Deserialize)]
struct PrototypeFile {
    prototypes: Vec<Prototype>,
}

/// One prototype entry.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct Prototype {
    /// Short Strukturbericht/AFLOW label, for example `B1`.
    pub label: String,
    /// Common name, for example `NaCl`.
    pub name: String,
    /// ITA space-group number.
    pub space_group: u32,
    /// Lattice parameters in Angstrom and degrees.
    pub lattice: PrototypeLattice,
    /// Occupied sites in asymmetric-unit order.
    pub sites: Vec<PrototypeSite>,
}

/// Lattice parameters of a prototype.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct PrototypeLattice {
    /// `a` in Angstrom.
    pub a: f64,
    /// `c / a` ratio; `1.0` when the cell is cubic.
    #[serde(default = "one")]
    pub c_over_a: f64,
    /// `alpha` in degrees.
    #[serde(default = "ninety")]
    pub alpha: f64,
    /// `beta` in degrees.
    #[serde(default = "ninety")]
    pub beta: f64,
    /// `gamma` in degrees.
    #[serde(default = "ninety")]
    pub gamma: f64,
}

fn one() -> f64 {
    1.0
}

fn ninety() -> f64 {
    90.0
}

/// One asymmetric site of a prototype.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct PrototypeSite {
    /// Element symbol.
    pub element: String,
    /// Wyckoff letter.
    pub wyckoff: String,
    /// Free parameters `x`, `y`, `z` for templates that use them.
    #[serde(default)]
    pub free: [f64; 3],
}

/// Options for prototype building and identification.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct PrototypeOptions {
    /// Uniform scale applied to the prototype lattice.
    pub scale: f64,
    /// Matcher tolerances used by [`identify_prototype`].
    pub matcher: MatcherOptions,
}

impl Default for PrototypeOptions {
    fn default() -> Self {
        Self {
            scale: 1.0,
            matcher: MatcherOptions::default(),
        }
    }
}

/// The curated prototype library.
///
/// # Errors
///
/// Returns [`Error::DataJson`] when the embedded file is invalid.
pub fn prototype_library() -> Result<Vec<Prototype>, Error> {
    let file: PrototypeFile =
        serde_json::from_str(PROTOTYPES_JSON).map_err(|error| Error::DataJson {
            name: "xtalforge-generation/prototypes.json".to_owned(),
            reason: error.to_string(),
        })?;
    Ok(file.prototypes)
}

/// Prototype labels in file order.
///
/// # Errors
///
/// Returns the errors of [`prototype_library`].
pub fn prototype_names() -> Result<Vec<String>, Error> {
    Ok(prototype_library()?
        .into_iter()
        .map(|prototype| prototype.label)
        .collect())
}

/// Builds one prototype by label.
///
/// # Errors
///
/// Returns [`Error::InvalidInput`] for an unknown label or invalid scale, and
/// the errors of the structure builder.
pub fn build_prototype(label: &str, options: &PrototypeOptions) -> Result<Structure, Error> {
    if !options.scale.is_finite() || options.scale <= 0.0 {
        return Err(Error::InvalidInput(
            "the prototype scale must be positive and finite".to_owned(),
        ));
    }
    let library = prototype_library()?;
    let prototype = library
        .iter()
        .find(|prototype| prototype.label.eq_ignore_ascii_case(label))
        .ok_or_else(|| {
            Error::InvalidInput(format!(
                "unknown prototype `{label}`; known labels: {}",
                library
                    .iter()
                    .map(|entry| entry.label.as_str())
                    .collect::<Vec<_>>()
                    .join(" ")
            ))
        })?;
    build_entry(prototype, options.scale)
}

fn build_entry(prototype: &Prototype, scale: f64) -> Result<Structure, Error> {
    let group = SpaceGroup::from_number(prototype.space_group)?;
    let positions = WyckoffPosition::list(group.hall_number)?;
    let lattice = Lattice::from_parameters(
        prototype.lattice.a * scale,
        prototype.lattice.a * scale,
        prototype.lattice.a * prototype.lattice.c_over_a * scale,
        prototype.lattice.alpha,
        prototype.lattice.beta,
        prototype.lattice.gamma,
    )?;
    let mut builder = StructureBuilder::new(prototype.label.clone(), lattice)
        .hall(group.hall_symbol.clone())
        .metadata("prototype", serde_json::json!(prototype.label))
        .metadata("prototype_name", serde_json::json!(prototype.name));
    for operation in &group.operations {
        builder = builder.operation(*operation);
    }
    for (index, site) in prototype.sites.iter().enumerate() {
        let Some(position) = positions
            .iter()
            .find(|position| position.letter.to_string() == site.wyckoff)
        else {
            return Err(Error::InvalidInput(format!(
                "prototype `{}` uses Wyckoff letter `{}` absent from space group {}",
                prototype.label, site.wyckoff, prototype.space_group
            )));
        };
        let point = evaluate_template(&position.coordinates, site.free)?;
        builder = builder.site(Site::new(
            format!("{}{}", site.element, index + 1),
            &site.element,
            point,
        )?);
    }
    builder.build()
}

/// Result of a prototype identification.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct PrototypeIdentification {
    /// Matched label, or `None` when nothing matched.
    pub label: Option<String>,
    /// Common name of the match.
    pub name: Option<String>,
    /// Detected ITA number.
    pub space_group: u32,
    /// Wyckoff signature of the input, for example `4a 4b`.
    pub signature: String,
    /// Matcher RMS in Angstrom for the best candidate.
    pub rms: Option<f64>,
    /// Candidate labels with the same space group and signature.
    pub candidates: Vec<String>,
}

/// Assigns the prototype label of an input structure.
///
/// The signature is the space group number, the sorted `(multiplicity,
/// Wyckoff letter)` pairs and the sorted stoichiometry. Candidates with that
/// signature are built and matched with the structure matcher; the smallest
/// RMS wins.
///
/// # Errors
///
/// Returns the errors of the symmetry search and the matcher.
pub fn identify_prototype(
    structure: &Structure,
    options: &PrototypeOptions,
) -> Result<PrototypeIdentification, Error> {
    let dataset = detect_symmetry(structure, SymmetryOptions::default())?;
    let assignments = wyckoff_assignments(structure, SymmetryOptions::default())?;
    let mut pairs: Vec<String> = assignments
        .iter()
        .map(|assignment| format!("{}{}", assignment.multiplicity, assignment.letter))
        .collect();
    pairs.sort();
    pairs.dedup();
    let signature = pairs.join(" ");
    let composition = stoichiometry_signature(structure)?;
    let library = prototype_library()?;
    let mut candidates = Vec::new();
    for prototype in &library {
        if prototype.space_group != dataset.number {
            continue;
        }
        if prototype_signature(prototype)? != signature {
            continue;
        }
        let built = build_entry(prototype, 1.0)?;
        if stoichiometry_signature(&built)? != composition {
            continue;
        }
        candidates.push(prototype.clone());
    }
    let mut best: Option<(f64, Prototype)> = None;
    for prototype in &candidates {
        let build = match build_entry(prototype, options.scale) {
            Ok(build) => build,
            Err(_) => continue,
        };
        let mapped = remap_species(structure, &build)?;
        let matched = match_structures(structure, &mapped, &options.matcher)?;
        if matched.matched && best.as_ref().is_none_or(|(rms, _)| matched.rms < *rms) {
            best = Some((matched.rms, prototype.clone()));
        }
    }
    Ok(PrototypeIdentification {
        label: best.as_ref().map(|(_, prototype)| prototype.label.clone()),
        name: best.as_ref().map(|(_, prototype)| prototype.name.clone()),
        space_group: dataset.number,
        signature,
        rms: best.map(|(rms, _)| rms),
        candidates: candidates
            .iter()
            .map(|prototype| prototype.label.clone())
            .collect(),
    })
}

fn prototype_signature(prototype: &Prototype) -> Result<String, Error> {
    let positions =
        WyckoffPosition::list(SpaceGroup::from_number(prototype.space_group)?.hall_number)?;
    let mut pairs = Vec::new();
    for site in &prototype.sites {
        let Some(position) = positions
            .iter()
            .find(|position| position.letter.to_string() == site.wyckoff)
        else {
            return Ok(String::new());
        };
        pairs.push(format!("{}{}", position.multiplicity, position.letter));
    }
    pairs.sort();
    pairs.dedup();
    Ok(pairs.join(" "))
}

fn stoichiometry_signature(structure: &Structure) -> Result<String, Error> {
    let p1 = structure.p1()?;
    let mut counts: BTreeMap<String, usize> = BTreeMap::new();
    for atom in &p1.atoms {
        *counts.entry(atom.element.clone()).or_default() += 1;
    }
    let mut parts: Vec<usize> = counts.into_values().collect();
    parts.sort_unstable();
    Ok(parts
        .into_iter()
        .map(|count| count.to_string())
        .collect::<Vec<_>>()
        .join(" "))
}

/// Rewrites the built prototype's species to the input's species in
/// first-appearance order, keeping the stoichiometry.
fn remap_species(input: &Structure, prototype: &Structure) -> Result<Structure, Error> {
    let mut input_order: Vec<String> = Vec::new();
    let input_p1 = input.p1()?;
    for atom in &input_p1.atoms {
        if !input_order.contains(&atom.element) {
            input_order.push(atom.element.clone());
        }
    }
    let mut prototype_order: Vec<String> = Vec::new();
    let prototype_p1 = prototype.p1()?;
    for atom in &prototype_p1.atoms {
        if !prototype_order.contains(&atom.element) {
            prototype_order.push(atom.element.clone());
        }
    }
    let mut mapping = BTreeMap::new();
    for (index, element) in prototype_order.iter().enumerate() {
        if let Some(target) = input_order.get(index) {
            mapping.insert(element.clone(), target.clone());
        }
    }
    let mut builder =
        StructureBuilder::new(prototype.name(), *prototype.lattice()).hall(prototype.hall());
    for operation in prototype.operations() {
        builder = builder.operation(*operation);
    }
    for site in prototype.sites() {
        let element = mapping
            .get(site.element())
            .cloned()
            .unwrap_or_else(|| site.element().to_owned());
        builder = builder.site(
            Site::new(site.label(), &element, site.frac())?.with_occupancy(site.occupancy())?,
        );
    }
    builder.build()
}
