//! Hydrogen placement (Task 4.6.1).
//!
//! X-ray refinements usually omit hydrogen atoms because they scatter so few
//! electrons. The missing atoms are rebuilt from the perceived bond graph:
//! the target valence decides how many are needed and the existing bonded
//! directions decide where they go. Geometry is periodic: every direction is
//! taken through the lattice minimum image.

use std::collections::{BTreeMap, BTreeSet};

use serde::{Deserialize, Serialize};
use xtalforge_core::elements::element_by_symbol;
use xtalforge_core::{BondGraph, Error, Lattice, P1Cell, Site, Structure, StructureBuilder};

pub(crate) mod geometry;

use geometry::{hydrogen_directions, normalized};

/// Default carbon-hydrogen distance in Angstrom.
const DEFAULT_C_H: f64 = 1.09;
/// Default nitrogen-hydrogen distance in Angstrom.
const DEFAULT_N_H: f64 = 1.01;
/// Default oxygen-hydrogen distance in Angstrom.
const DEFAULT_O_H: f64 = 0.96;
/// Ideal tetrahedral angle in degrees.
const TETRAHEDRAL: f64 = 109.47;
/// Ideal trigonal angle in degrees.
const TRIGONAL: f64 = 120.0;
/// Ideal bent angle in degrees.
const BENT: f64 = 104.5;
/// Below this length a direction is considered undefined.
const DEGENERATE: f64 = 1e-9;

/// Options for [`add_hydrogens`].
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct HydrogenOptions {
    /// X-H bond length in Angstrom.
    ///
    /// The default is the C-H length, 1.09 A. While the field keeps that
    /// value the element-specific defaults apply: N-H 1.01 A and O-H
    /// 0.96 A (boron, phosphorus and sulfur use this field). Setting any
    /// other value uses it for every element.
    pub distance: f64,
    /// Ideal X-H bond angle in degrees, keyed by canonical element symbol.
    ///
    /// The defaults are tetrahedral for carbon (109.47), trigonal for
    /// nitrogen (120) and bent for oxygen (104.5). Boron and phosphorus
    /// fall back to 120 and sulfur to 104.5 when no entry is given.
    pub angles: BTreeMap<String, f64>,
}

impl Default for HydrogenOptions {
    fn default() -> Self {
        let mut angles = BTreeMap::new();
        angles.insert("C".to_owned(), TETRAHEDRAL);
        angles.insert("N".to_owned(), TRIGONAL);
        angles.insert("O".to_owned(), BENT);
        Self {
            distance: DEFAULT_C_H,
            angles,
        }
    }
}

/// Adds the hydrogens missing from the perceived coordination.
///
/// Every asymmetric site with a target valence (C 4, N 3, O 2, B 3, P 3,
/// S 2) receives one hydrogen per unused valence, placed in directions
/// derived from its existing bonded neighbours through the lattice minimum
/// image:
///
/// - no neighbour: `+x`, `+y`, `+z`, or the four tetrahedral directions;
/// - one neighbour: the ideal angle, distributed on a deterministic cone;
/// - two neighbours: the ideal angle completing the trigonal or bent
///   geometry;
/// - three neighbours: opposite the sum of the unit bond vectors.
///
/// The original name, lattice, Hall symbol, operations, explicit and
/// topology bonds, sites and metadata are preserved. New sites are labelled
/// `H1`, `H2`, ... with the first number not already used by the structure.
/// Atoms that already carry their full valence, hydrogen and deuterium, and
/// elements without a target valence are left alone.
///
/// # Errors
///
/// Returns [`Error::InvalidInput`] for a non-positive or non-finite distance,
/// an angle outside `(0, 180)`, an element in `angles` without a target
/// valence, an unknown element in the structure, or a zero-length direction
/// basis (coincident atoms or exactly cancelling bond vectors).
pub fn add_hydrogens(structure: &Structure, options: &HydrogenOptions) -> Result<Structure, Error> {
    validate_options(options)?;
    let angles = canonical_angles(options)?;
    let p1 = structure.p1()?;
    let graph = structure.bond_graph()?;

    let mut builder = StructureBuilder::new(structure.name(), *structure.lattice())
        .hall(structure.hall().to_owned());
    for operation in structure.operations() {
        builder = builder.operation(*operation);
    }
    for bond in structure.explicit_bonds() {
        builder = builder.explicit_bond(*bond);
    }
    for bond in structure.topology_bonds() {
        builder = builder.topology_bond(*bond);
    }
    for (key, value) in structure.metadata() {
        builder = builder.metadata(key.clone(), value.clone());
    }
    let mut used_labels = BTreeSet::new();
    for site in structure.sites() {
        element_by_symbol(site.element())?;
        used_labels.insert(site.label().to_owned());
        builder = builder.site(site.clone());
    }

    let mut counter = 0_usize;
    for (site_index, site) in structure.sites().iter().enumerate() {
        let anchor = Lattice::wrap(site.frac());
        let atom_index = p1
            .atoms
            .iter()
            .position(|atom| {
                atom.site_idx == site_index && p1.lattice.distance(atom.frac, anchor) <= 0.05
            })
            .or_else(|| p1.atoms.iter().position(|atom| atom.site_idx == site_index));
        let Some(atom_index) = atom_index else {
            continue;
        };
        let Some(valence) = target_valence(site.element()) else {
            continue;
        };
        let (directions, used_order) = bonded_directions(&p1, &graph, atom_index)?;
        let missing = usize::from(valence).saturating_sub(round_order(used_order));
        if missing == 0 {
            continue;
        }
        let angle = ideal_angle(site.element(), &angles);
        let targets = hydrogen_directions(&directions, missing, angle)?;
        let length = bond_length(site.element(), options);
        let Some(atom) = p1.atoms.get(atom_index) else {
            continue;
        };
        let origin = p1.lattice.to_cart(atom.frac);
        for direction in targets {
            let cart = [
                origin[0] + direction[0] * length,
                origin[1] + direction[1] * length,
                origin[2] + direction[2] * length,
            ];
            let label = next_label(&mut used_labels, &mut counter);
            builder = builder.site(Site::new(label, "H", p1.lattice.to_frac(cart))?);
        }
    }
    builder.build()
}

/// Target valence of the elements this task covers.
fn target_valence(element: &str) -> Option<u8> {
    match element {
        "C" => Some(4),
        "N" | "B" | "P" => Some(3),
        "O" | "S" => Some(2),
        _ => None,
    }
}

/// Checks the option values and every requested angle element.
fn validate_options(options: &HydrogenOptions) -> Result<(), Error> {
    if !options.distance.is_finite() || options.distance <= 0.0 {
        return Err(Error::InvalidInput(
            "hydrogen distance must be positive and finite".to_owned(),
        ));
    }
    for (element, angle) in &options.angles {
        if !angle.is_finite() || *angle <= 0.0 || *angle >= 180.0 {
            return Err(Error::InvalidInput(format!(
                "hydrogen angle for {element} must be inside (0, 180) degrees, got {angle}"
            )));
        }
        let symbol = element_by_symbol(element)?.symbol.clone();
        if target_valence(&symbol).is_none() {
            return Err(Error::InvalidInput(format!(
                "element {symbol} has no target hydrogen valence"
            )));
        }
    }
    Ok(())
}

/// Canonicalizes the angle keys and rejects two spellings of one element.
fn canonical_angles(options: &HydrogenOptions) -> Result<BTreeMap<String, f64>, Error> {
    let mut canonical = BTreeMap::new();
    for (element, angle) in &options.angles {
        let symbol = element_by_symbol(element)?.symbol.clone();
        if canonical.insert(symbol.clone(), *angle).is_some() {
            return Err(Error::InvalidInput(format!(
                "hydrogen angle for {symbol} is given more than once"
            )));
        }
    }
    Ok(canonical)
}

/// Bond length for one element.
fn bond_length(element: &str, options: &HydrogenOptions) -> f64 {
    if options.distance != DEFAULT_C_H {
        return options.distance;
    }
    match element {
        "N" => DEFAULT_N_H,
        "O" => DEFAULT_O_H,
        _ => options.distance,
    }
}

/// Ideal X-H angle for one element.
fn ideal_angle(element: &str, angles: &BTreeMap<String, f64>) -> f64 {
    if let Some(angle) = angles.get(element) {
        return *angle;
    }
    match element {
        "B" | "P" => TRIGONAL,
        "S" => BENT,
        _ => TETRAHEDRAL,
    }
}

/// Hydrogens implied by a SMILES atom with `used_order` already bonded.
///
/// The sum of the bond orders (aromatic bonds count as `1.5`) is rounded to
/// whole bonds exactly as [`add_hydrogens`] does. The formal charge shifts the
/// target valence, so `[NH4+]` keeps four hydrogens and `[OH-]` one. Returns
/// `None` for elements outside the target-valence table.
pub(crate) fn inference_missing(element: &str, charge: i32, used_order: f64) -> Option<usize> {
    let target = i32::from(target_valence(element)?) + charge;
    if target <= 0 {
        return Some(0);
    }
    Some((target as usize).saturating_sub(round_order(used_order)))
}

/// Ideal X-H length in Angstrom and X-H angle in degrees for inference.
pub(crate) fn inference_geometry(element: &str) -> (f64, f64) {
    let options = HydrogenOptions::default();
    (
        bond_length(element, &options),
        ideal_angle(element, &options.angles),
    )
}

/// Number of valence units already used, rounded to whole bonds.
fn round_order(used: f64) -> usize {
    let rounded = used.round();
    if rounded.is_finite() && rounded > 0.0 {
        rounded as usize
    } else {
        0
    }
}

/// First free `H<n>` label, in deterministic order.
fn next_label(used: &mut BTreeSet<String>, counter: &mut usize) -> String {
    loop {
        *counter += 1;
        let label = format!("H{counter}");
        if used.insert(label.clone()) {
            return label;
        }
    }
}

/// Unit vectors from one P1 atom to its bonded neighbours and their orders.
fn bonded_directions(
    p1: &P1Cell,
    graph: &BondGraph,
    index: usize,
) -> Result<(Vec<[f64; 3]>, f64), Error> {
    let Some(atom) = p1.atoms.get(index) else {
        return Err(Error::InvalidInput(format!(
            "atom {index} is outside the P1 cell"
        )));
    };
    let mut directions = Vec::new();
    let mut used_order = 0.0;
    for bond_index in graph.bonds_of(index) {
        let Some(bond) = graph.bonds.get(*bond_index) else {
            continue;
        };
        used_order += bond.order;
        let other = if bond.i == index { bond.j } else { bond.i };
        let Some(neighbour) = p1.atoms.get(other) else {
            continue;
        };
        let delta = [
            neighbour.frac[0] - atom.frac[0],
            neighbour.frac[1] - atom.frac[1],
            neighbour.frac[2] - atom.frac[2],
        ];
        let (image, distance) = p1.lattice.minimum_image(delta);
        if distance <= DEGENERATE {
            return Err(Error::InvalidInput(format!(
                "zero-length direction basis at atom {index}"
            )));
        }
        directions.push(normalized(p1.lattice.to_cart(image))?);
    }
    Ok((directions, used_order))
}
