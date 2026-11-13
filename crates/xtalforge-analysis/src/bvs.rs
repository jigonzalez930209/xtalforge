#![allow(clippy::indexing_slicing, clippy::needless_range_loop)] // fixed-size numeric kernels

//! Bond-valence sums and oxidation-state estimation (Task 9.3.3).
//!
//! Parameters follow O'Keefe and Brese, *JACS* 113 (1991) 3226; the
//! data-mined ICSD distributions are the ones pymatgen ships (MIT) and are
//! recorded in `THIRD_PARTY.md`. The estimator picks one oxidation state per
//! element by a maximum-a-posteriori search with a charge-balance penalty.

use std::collections::BTreeMap;
use std::sync::OnceLock;

use serde::{Deserialize, Serialize};

use xtalforge_core::data::EmbeddedBundle;
use xtalforge_core::{Error, Structure, neighbors};

const PARAMS_JSON: &str = include_str!("../data/bvparam_1991.json");
/// The gzip table is decompressed and checksum-verified on first use
/// (Task 15.3.2).
static ICSD_JSON: EmbeddedBundle = EmbeddedBundle::new(
    "xtalforge-analysis/icsd_bv.json.gz",
    include_bytes!("../data/icsd_bv.json.gz"),
    "aab60250ec0e66e05d5a19843d10dc6704eca9802d3f7007667d37055574b5d7",
);
const ELECTRONEGATIVE: [&str; 16] = [
    "H", "B", "C", "Si", "N", "P", "As", "Sb", "O", "S", "Se", "Te", "F", "Cl", "Br", "I",
];

/// One bond-valence parameter pair.
#[derive(Debug, Clone, Copy, PartialEq, Deserialize)]
struct BvParameter {
    c: f64,
    r: f64,
}

/// One ICSD bond-valence distribution.
#[derive(Debug, Clone, Copy, PartialEq, Deserialize)]
struct IcsdData {
    mean: f64,
    std: f64,
}

#[derive(Debug, Clone, Deserialize)]
struct IcsdFile {
    bvsum: BTreeMap<String, IcsdData>,
    occurrence: BTreeMap<String, f64>,
}

/// Options for the bond-valence analysis.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct BvsOptions {
    /// Neighbor cutoff in Angstrom.
    pub cutoff: f64,
    /// Distance scale factor applied before the exponential.
    pub scale: f64,
}

impl Default for BvsOptions {
    fn default() -> Self {
        Self {
            cutoff: 5.0,
            scale: 1.0,
        }
    }
}

/// One atom of a bond-valence report.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct BondValenceSite {
    /// P1 atom index.
    pub index: usize,
    /// Element symbol.
    pub element: String,
    /// Signed bond-valence sum.
    pub bvs: f64,
    /// Number of bonds that contributed.
    pub bonds: usize,
}

/// Bond-valence sums for every atom.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct BondValenceReport {
    /// Per-atom sums.
    pub sites: Vec<BondValenceSite>,
}

/// One estimated oxidation state.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct OxidationSite {
    /// P1 atom index.
    pub index: usize,
    /// Element symbol.
    pub element: String,
    /// Estimated oxidation state.
    pub state: i32,
    /// Bond-valence sum used by the estimate.
    pub bvs: f64,
}

/// Oxidation-state estimation result.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct OxidationReport {
    /// Per-atom states.
    pub sites: Vec<OxidationSite>,
    /// One state per element.
    pub by_element: BTreeMap<String, i32>,
    /// Total charge of the cell.
    pub total_charge: f64,
    /// Sum of the negative log posterior and charge penalty.
    pub score: f64,
    /// Whether the total charge is neutral within `1e-4`.
    pub neutral: bool,
}

/// Computes bond-valence sums for every P1 atom.
///
/// # Errors
///
/// Returns [`Error::InvalidInput`] for an invalid cutoff or scale, or the
/// errors of the P1 expansion and neighbor search.
pub fn bond_valence_sums(
    structure: &Structure,
    options: &BvsOptions,
) -> Result<BondValenceReport, Error> {
    if !options.cutoff.is_finite() || options.cutoff <= 0.0 {
        return Err(Error::InvalidInput(
            "bond-valence cutoff must be positive and finite".to_owned(),
        ));
    }
    if !options.scale.is_finite() || options.scale <= 0.0 {
        return Err(Error::InvalidInput(
            "bond-valence distance scale must be positive and finite".to_owned(),
        ));
    }
    let parameters = parameters()?;
    let p1 = structure.p1()?;
    let snapshot = p1_structure(&p1)?;
    let pairs = neighbors(&snapshot, options.cutoff)?;
    let mut sites: Vec<BondValenceSite> = p1
        .atoms
        .iter()
        .enumerate()
        .map(|(index, atom)| BondValenceSite {
            index,
            element: atom.element.clone(),
            bvs: 0.0,
            bonds: 0,
        })
        .collect();
    for pair in pairs {
        let (Some(first), Some(second)) = (p1.atoms.get(pair.i), p1.atoms.get(pair.j)) else {
            continue;
        };
        let Some(value) = pair_valence(
            &first.element,
            &second.element,
            pair.distance * options.scale,
            parameters,
        ) else {
            continue;
        };
        if let Some(entry) = sites.get_mut(pair.i) {
            entry.bvs += value;
            entry.bonds += 1;
        }
        if let Some(entry) = sites.get_mut(pair.j) {
            entry.bvs -= value;
            entry.bonds += 1;
        }
    }
    Ok(BondValenceReport { sites })
}

/// Estimates one oxidation state per element from bond-valence sums.
///
/// # Errors
///
/// Returns the errors of [`bond_valence_sums`] or [`Error::InvalidInput`]
/// when an element has no ICSD distribution.
pub fn oxidation_states(
    structure: &Structure,
    options: &BvsOptions,
) -> Result<OxidationReport, Error> {
    let report = bond_valence_sums(structure, options)?;
    let data = icsd()?;
    let mut grouped: BTreeMap<String, Vec<f64>> = BTreeMap::new();
    for site in &report.sites {
        grouped
            .entry(site.element.clone())
            .or_default()
            .push(site.bvs);
    }
    let mut candidates: BTreeMap<String, Vec<(i32, f64)>> = BTreeMap::new();
    for (element, values) in &grouped {
        let mean = values.iter().sum::<f64>() / values.len().max(1) as f64;
        let mut options: Vec<(i32, f64)> = Vec::new();
        for (key, distribution) in &data.bvsum {
            let Some((symbol, state)) = parse_species(key) else {
                continue;
            };
            if symbol != *element || distribution.std <= 0.0 {
                continue;
            }
            let prior = data.occurrence.get(key).copied().unwrap_or(0.0);
            if prior <= 0.0 {
                continue;
            }
            let probability = (-((mean - distribution.mean).powi(2))
                / (2.0 * distribution.std * distribution.std))
                .exp()
                / distribution.std
                * prior;
            options.push((state, -probability.max(1e-300).ln()));
        }
        if options.is_empty() {
            return Err(Error::InvalidInput(format!(
                "no oxidation-state distribution for {element}"
            )));
        }
        options.sort_by(|first, second| {
            second
                .1
                .partial_cmp(&first.1)
                .unwrap_or(std::cmp::Ordering::Equal)
                .then_with(|| first.0.cmp(&second.0))
        });
        options.truncate(5);
        candidates.insert(element.clone(), options);
    }
    let elements: Vec<String> = grouped.keys().cloned().collect();
    let mut best: Option<(f64, BTreeMap<String, i32>)> = None;
    let mut current = BTreeMap::new();
    let mut nodes = 0usize;
    search(
        &elements,
        &candidates,
        &mut current,
        0,
        0.0,
        0.0,
        &mut best,
        &mut nodes,
    );
    let Some((score, states)) = best else {
        return Err(Error::InvalidInput(
            "no oxidation-state combination was found".to_owned(),
        ));
    };
    let mut sites = Vec::with_capacity(report.sites.len());
    let mut total_charge = 0.0;
    for site in &report.sites {
        let state = states.get(&site.element).copied().unwrap_or(0);
        total_charge += f64::from(state);
        sites.push(OxidationSite {
            index: site.index,
            element: site.element.clone(),
            state,
            bvs: site.bvs,
        });
    }
    Ok(OxidationReport {
        sites,
        by_element: states,
        total_charge,
        score,
        neutral: total_charge.abs() < 1e-4,
    })
}

#[allow(clippy::too_many_arguments)] // recursive MAP search context
fn search(
    elements: &[String],
    candidates: &BTreeMap<String, Vec<(i32, f64)>>,
    current: &mut BTreeMap<String, i32>,
    position: usize,
    charge: f64,
    score: f64,
    best: &mut Option<(f64, BTreeMap<String, i32>)>,
    nodes: &mut usize,
) {
    *nodes += 1;
    if *nodes > 500_000 {
        return;
    }
    if position >= elements.len() {
        let penalty = 10.0 * charge.abs();
        let total = score + penalty;
        let better = best
            .as_ref()
            .is_none_or(|(current_score, _)| total < *current_score);
        if better {
            *best = Some((total, current.clone()));
        }
        return;
    }
    let Some(element) = elements.get(position) else {
        return;
    };
    let Some(options) = candidates.get(element) else {
        return;
    };
    let maximum_remaining: i32 = elements
        .iter()
        .skip(position + 1)
        .filter_map(|name| candidates.get(name))
        .map(|entries| {
            entries
                .iter()
                .map(|(state, _)| state.abs())
                .max()
                .unwrap_or(0)
        })
        .sum();
    for (state, log_probability) in options {
        let next_charge = charge + f64::from(*state);
        if next_charge.abs() > f64::from(maximum_remaining) + 1.0 {
            continue;
        }
        current.insert(element.clone(), *state);
        search(
            elements,
            candidates,
            current,
            position + 1,
            next_charge,
            score + log_probability,
            best,
            nodes,
        );
        current.remove(element);
    }
}

fn pair_valence(
    first: &str,
    second: &str,
    distance: f64,
    parameters: &BTreeMap<String, BvParameter>,
) -> Option<f64> {
    if first == second {
        return None;
    }
    let first_electronegative = ELECTRONEGATIVE.contains(&first);
    let second_electronegative = ELECTRONEGATIVE.contains(&second);
    if !first_electronegative && !second_electronegative {
        return None;
    }
    let first_parameter = parameters.get(first)?;
    let second_parameter = parameters.get(second)?;
    let c1 = first_parameter.c;
    let c2 = second_parameter.c;
    let r1 = first_parameter.r;
    let r2 = second_parameter.r;
    let denominator = c1 * r1 + c2 * r2;
    if denominator <= 0.0 {
        return None;
    }
    let r0 = r1 + r2 - r1 * r2 * (c1.sqrt() - c2.sqrt()).powi(2) / denominator;
    let magnitude = ((r0 - distance) / 0.31).exp();
    let first_x = xtalforge_core::element_info(first)
        .ok()
        .and_then(|element| element.electronegativity)
        .unwrap_or(0.0);
    let second_x = xtalforge_core::element_info(second)
        .ok()
        .and_then(|element| element.electronegativity)
        .unwrap_or(0.0);
    let sign = if first_x < second_x { 1.0 } else { -1.0 };
    Some(magnitude * sign)
}

fn parse_species(key: &str) -> Option<(String, i32)> {
    let mut letters = String::new();
    let mut rest = key;
    for (index, character) in key.char_indices() {
        if character.is_ascii_alphabetic() && rest == key {
            letters.push(character);
        } else {
            rest = key.get(index..).unwrap_or("");
            break;
        }
    }
    if letters.is_empty() {
        return None;
    }
    let sign = if rest.contains('-') { -1 } else { 1 };
    let digits: String = rest
        .chars()
        .filter(|character| character.is_ascii_digit())
        .collect();
    let magnitude = if digits.is_empty() {
        1
    } else {
        digits.parse::<i32>().ok()?
    };
    Some((letters, sign * magnitude))
}

fn parameters() -> Result<&'static BTreeMap<String, BvParameter>, Error> {
    static PARAMS: OnceLock<Result<BTreeMap<String, BvParameter>, Error>> = OnceLock::new();
    let cached = PARAMS.get_or_init(|| {
        serde_json::from_str(PARAMS_JSON)
            .map_err(|error| Error::InvalidInput(format!("invalid bond-valence data: {error}")))
    });
    match cached {
        Ok(parameters) => Ok(parameters),
        Err(error) => Err(error.clone()),
    }
}

fn icsd() -> Result<&'static IcsdFile, Error> {
    static DATA: OnceLock<Result<IcsdFile, Error>> = OnceLock::new();
    let cached = DATA.get_or_init(|| {
        let payload = ICSD_JSON.bytes()?;
        serde_json::from_slice(payload)
            .map_err(|error| Error::InvalidInput(format!("invalid ICSD data: {error}")))
    });
    match cached {
        Ok(data) => Ok(data),
        Err(error) => Err(error.clone()),
    }
}

fn p1_structure(p1: &xtalforge_core::P1Cell) -> Result<Structure, Error> {
    let mut builder = xtalforge_core::StructureBuilder::new("p1", p1.lattice);
    for atom in &p1.atoms {
        builder = builder.site(
            xtalforge_core::Site::new(atom.label.clone(), &atom.element, atom.frac)?
                .with_occupancy(atom.occupancy)?,
        );
    }
    builder.build()
}
