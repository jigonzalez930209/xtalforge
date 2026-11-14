#![allow(clippy::indexing_slicing)] // fixed three-component vectors and bounded orbit arrays

//! Low-dimensional generators (Task 20.1.3).
//!
//! # Layers
//!
//! [`generate_layer_crystal`] uses the 80 layer groups of moyo. The cell is
//! slab-shaped: `a` and `b` follow the two-dimensional Bravais class and `c`
//! is a 20 A vacuum gap, so the layer detector of Task 2.8.1 sees an isolated
//! layer. Candidates are accepted when the periodic minimum distance holds and
//! [`detect_layer_symmetry`] returns
//! the requested layer-group number.
//!
//! # Rods
//!
//! moyo has no rod-group tables. xtalforge ships a documented catalogue of
//! fourteen common rod groups as explicit operation sets (identity, inversion,
//! two-, three-, four- and six-fold rotations about the chain axis and their
//! mirrors). A rod is generated in a wide cell (`a = b = 20 A`) with the chain
//! repeated along `c`; the requested operation set is verified directly on the
//! expanded atoms, so this is a real symmetry check rather than a 3D
//! space-group approximation.

use std::collections::BTreeSet;

use serde::{Deserialize, Serialize};
use serde_json::json;
use xtalforge_core::{Error, Lattice, Site, Structure, StructureBuilder};
use xtalforge_symmetry::{
    LayerGroup, LayerWyckoffPosition, SymmetryOptions, detect_layer_symmetry,
};

use crate::geometry::{evaluate_template, minimum_distance, operations_preserve, orbit, wrap};
use crate::rng::SplitMix64;
use crate::wyckoff::{Selection, choose_orbit_letters};

/// An optical gap used for the aperiodic stacking direction.
pub const LAYER_VACUUM: f64 = 20.0;

/// Options for [`generate_layer_crystal`].
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct LayerCrystalOptions {
    /// Layer-group number (`1` to `80`).
    pub layer_group: u32,
    /// Element symbols and counts in the conventional layer cell.
    pub composition: Vec<(String, u32)>,
    /// Target area per atom in square Angstrom.
    pub area_per_atom: f64,
    /// Smallest allowed periodic distance in Angstrom.
    pub min_distance: f64,
    /// Number of attempts before the generator reports failure.
    pub max_attempts: usize,
    /// Explicit seed.
    pub seed: u64,
}

impl LayerCrystalOptions {
    /// Options with the documented defaults (10 A^2 per atom, 1.0 A minimum).
    #[must_use]
    pub fn new(layer_group: u32, composition: Vec<(String, u32)>, seed: u64) -> Self {
        Self {
            layer_group,
            composition,
            area_per_atom: 10.0,
            min_distance: 1.0,
            max_attempts: 64,
            seed,
        }
    }
}

/// Result of one layer generation.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct LowDimOptions {
    /// Group number.
    pub group: u32,
    /// Wyckoff letters used.
    pub wyckoffs: Vec<String>,
    /// Seed.
    pub seed: u64,
}

/// Result of one layer generation.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct LowDimCrystal {
    /// The generated layer.
    pub structure: Structure,
    /// Layer-group number.
    pub layer_group: u32,
    /// Hermann-Mauguin symbol.
    pub hm_symbol: String,
    /// Wyckoff letters used.
    pub wyckoffs: Vec<String>,
    /// Attempt that produced the structure.
    pub attempt: usize,
    /// Seed that produced the structure.
    pub seed: u64,
}

impl LowDimCrystal {
    /// Serializes the result to JSON.
    ///
    /// # Errors
    ///
    /// Returns [`Error::InvalidInput`] when serialization fails.
    pub fn to_json(&self) -> Result<String, Error> {
        serde_json::to_string(self).map_err(|error| {
            Error::InvalidInput(format!(
                "could not serialize low-dimensional crystal: {error}"
            ))
        })
    }
}

/// The 80 layer groups.
///
/// # Errors
///
/// Returns the errors of [`xtalforge_symmetry::layer_groups()`].
pub fn layer_groups_catalogue() -> Result<Vec<LayerGroup>, Error> {
    xtalforge_symmetry::layer_groups()
}

/// Layer groups, re-exported for the bindings.
///
/// # Errors
///
/// Returns the errors of [`xtalforge_symmetry::layer_groups()`].
pub fn layer_groups() -> Result<Vec<LayerGroup>, Error> {
    xtalforge_symmetry::layer_groups()
}

/// Generates one random layer crystal.
///
/// # Errors
///
/// Returns [`Error::InvalidInput`] for invalid options or exhausted attempts.
pub fn generate_layer_crystal(options: &LayerCrystalOptions) -> Result<LowDimCrystal, Error> {
    validate_layer(options)?;
    let group = LayerGroup::from_number(options.layer_group)?;
    let positions = LayerWyckoffPosition::list(group.hall_number)?;
    let multiplicities: Vec<(usize, char, bool)> = positions
        .iter()
        .map(|position| {
            (
                position.multiplicity,
                position.letter,
                position.coordinates.contains(['x', 'y', 'z']),
            )
        })
        .collect();
    let mut rng = SplitMix64::new(options.seed);
    let mut last_reason = String::from("no candidate was tried");
    for attempt in 1..=options.max_attempts {
        let Some(selection) = choose_orbit_letters(&options.composition, &multiplicities, &mut rng)
        else {
            last_reason = "the composition cannot be expressed with the layer Wyckoff \
                multiplicities"
                .to_owned();
            continue;
        };
        let total: usize = selection
            .iter()
            .filter_map(|chosen| positions.get(chosen.position))
            .map(|position| position.multiplicity)
            .sum();
        let lattice = layer_lattice(&group.bravais, total, options.area_per_atom, &mut rng)?;
        match build_layer(&group, &positions, &selection, lattice, options, &mut rng) {
            Ok(Some((structure, letters))) => {
                return Ok(LowDimCrystal {
                    structure,
                    layer_group: group.number,
                    hm_symbol: group.hm_symbol,
                    wyckoffs: letters,
                    attempt,
                    seed: options.seed,
                });
            }
            Ok(None) => {
                last_reason =
                    "a candidate collapsed its orbit or violated the minimum distance".to_owned();
            }
            Err(error) => last_reason = error.message(),
        }
    }
    Err(Error::InvalidInput(format!(
        "could not generate layer group {} after {} attempts: {last_reason}",
        options.layer_group, options.max_attempts
    )))
}

fn layer_lattice(
    bravais: &str,
    atoms: usize,
    area_per_atom: f64,
    rng: &mut SplitMix64,
) -> Result<Lattice, Error> {
    if !area_per_atom.is_finite() || area_per_atom <= 0.0 {
        return Err(Error::InvalidInput(
            "area per atom must be positive and finite".to_owned(),
        ));
    }
    let area = area_per_atom * atoms.max(1) as f64;
    let c = LAYER_VACUUM;
    let (a, b, gamma) = match bravais {
        "tp" => {
            let side = area.sqrt();
            (side, side, 90.0)
        }
        "hp" => {
            let side = (area / (3.0_f64.sqrt() / 2.0)).sqrt();
            (side, side, 120.0)
        }
        "op" | "oc" => {
            let ratio = rng.range(0.7, 1.4);
            let a = (area / ratio).sqrt();
            (a, a * ratio, 90.0)
        }
        _ => {
            let ratio = rng.range(0.7, 1.4);
            let a = (area / ratio).sqrt();
            (a, a * ratio, rng.range(72.0, 108.0))
        }
    };
    Lattice::from_parameters(a, b, c, 90.0, 90.0, gamma)
}

fn build_layer(
    group: &LayerGroup,
    positions: &[LayerWyckoffPosition],
    selection: &[Selection],
    lattice: Lattice,
    options: &LayerCrystalOptions,
    rng: &mut SplitMix64,
) -> Result<Option<(Structure, Vec<String>)>, Error> {
    let mut builder = StructureBuilder::new(format!("layer-{}", group.number), lattice)
        .hall("P 1")
        .metadata("generator", json!("layer"))
        .metadata("layer_group", json!(group.number))
        .metadata("seed", json!(options.seed));
    for operation in &group.operations {
        builder = builder.operation(*operation);
    }
    let mut letters = Vec::new();
    let mut label_index = 0usize;
    for chosen in selection {
        let Some(position) = positions.get(chosen.position) else {
            return Ok(None);
        };
        let free = [rng.uniform(), rng.uniform(), rng.uniform()];
        let point = wrap(evaluate_template(&position.coordinates, free)?);
        let images = orbit(&[point], &group.operations, &lattice);
        if images.len() != position.multiplicity {
            return Ok(None);
        }
        label_index += 1;
        let site = Site::new(
            format!("{}{label_index}", chosen.element),
            &chosen.element,
            point,
        )?;
        builder = builder.site(site);
        letters.push(position.letter.to_string());
    }
    let structure = builder.build()?;
    if minimum_distance(&structure)? < options.min_distance {
        return Ok(None);
    }
    if operations_preserve(&structure, &group.operations)? {
        let dataset = detect_layer_symmetry(&structure, SymmetryOptions::default().symprec)?;
        if dataset.number == group.number {
            return Ok(Some((structure, letters)));
        }
    }
    Ok(None)
}

fn validate_layer(options: &LayerCrystalOptions) -> Result<(), Error> {
    if !(1..=80).contains(&options.layer_group) {
        return Err(Error::InvalidInput(format!(
            "layer-group number must be in 1..=80, got {}",
            options.layer_group
        )));
    }
    validate_common(
        &options.composition,
        options.min_distance,
        options.max_attempts,
    )
}

pub(crate) fn validate_common(
    composition: &[(String, u32)],
    min_distance: f64,
    max_attempts: usize,
) -> Result<(), Error> {
    if composition.is_empty() {
        return Err(Error::InvalidInput(
            "the composition needs at least one element".to_owned(),
        ));
    }
    let mut seen = BTreeSet::new();
    let mut total = 0u32;
    for (element, count) in composition {
        let info = xtalforge_core::elements::element_by_symbol(element)?;
        if !seen.insert(info.symbol.clone()) {
            return Err(Error::InvalidInput(format!(
                "element `{element}` is listed twice"
            )));
        }
        if *count == 0 {
            return Err(Error::InvalidInput(format!(
                "element `{element}` has a zero count"
            )));
        }
        total += *count;
    }
    if total > 400 {
        return Err(Error::InvalidInput(format!(
            "the conventional cell is capped at 400 atoms, got {total}"
        )));
    }
    if !min_distance.is_finite() || min_distance < 0.0 {
        return Err(Error::InvalidInput(
            "minimum distance must be finite and non-negative".to_owned(),
        ));
    }
    if max_attempts == 0 {
        return Err(Error::InvalidInput(
            "max_attempts must be at least one".to_owned(),
        ));
    }
    Ok(())
}
