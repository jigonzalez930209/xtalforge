#![allow(clippy::indexing_slicing)] // fixed three-component vectors and bounded atom indices

//! Rigid molecular crystals (Task 20.1.2).
//!
//! A molecule read from SMILES, a file or any P1 [`Structure`] is treated as a
//! rigid body: its internal geometry never changes. Each independent molecule
//! (Z') is placed on a Wyckoff position with a random orientation drawn from
//! the seeded generator, and the space-group operations produce the rest of
//! the cell.
//!
//! Special positions are accepted only when the molecule is compatible with
//! the site symmetry. The test is geometric: every operation that maps a
//! molecule centre onto itself must map the molecule onto itself within
//! `tolerance`, and no two symmetry-related molecules may approach closer
//! than `min_distance`. Candidates that fail either test are redrawn.
//!
//! The output keeps the requested Hall symbol and operations, so the normal
//! P1 expansion reproduces the crystal.

use std::collections::BTreeSet;

use serde::{Deserialize, Serialize};
use serde_json::json;
use xtalforge_core::{Error, Site, Structure, StructureBuilder};
use xtalforge_symmetry::{
    SpaceGroup, SymmetryOptions, WyckoffPosition, crystal_system, detect_symmetry,
};

use crate::geometry::{
    SAME_POINT, apply_matrix, evaluate_template, lattice_for_system, minimum_intergroup_distance,
    quaternion_matrix, wrap,
};
use crate::rng::SplitMix64;

/// Options for [`generate_molecular_crystal`].
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct MolecularCrystalOptions {
    /// ITA space-group number (`1` to `230`).
    pub space_group: u32,
    /// Rigid molecule; its own lattice is ignored.
    pub molecule: Structure,
    /// Independent molecules in the asymmetric unit (Z').
    pub copies: u32,
    /// Target volume per atom in cubic Angstrom.
    pub volume_per_atom: f64,
    /// Smallest allowed intermolecular distance in Angstrom.
    pub min_distance: f64,
    /// Distance tolerance for site-symmetry compatibility in Angstrom.
    pub tolerance: f64,
    /// Number of attempts before the generator reports failure.
    pub max_attempts: usize,
    /// Explicit seed.
    pub seed: u64,
}

impl MolecularCrystalOptions {
    /// Options with Z' = 1 and the documented defaults.
    #[must_use]
    pub fn new(space_group: u32, molecule: Structure, seed: u64) -> Self {
        Self {
            space_group,
            molecule,
            copies: 1,
            volume_per_atom: 20.0,
            min_distance: 1.0,
            tolerance: 0.25,
            max_attempts: 64,
            seed,
        }
    }
}

/// A generated molecular crystal.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct MolecularCrystal {
    /// The crystal.
    pub structure: Structure,
    /// Requested ITA number.
    pub space_group: u32,
    /// Wyckoff letters used by the independent molecules.
    pub wyckoffs: Vec<String>,
    /// Number of independent molecules.
    pub copies: u32,
    /// Total molecules in the conventional cell.
    pub molecules_per_cell: usize,
    /// Seed that produced the crystal.
    pub seed: u64,
    /// Attempt that produced the crystal (1-based).
    pub attempt: usize,
}

impl MolecularCrystal {
    /// Serializes the result to JSON.
    ///
    /// # Errors
    ///
    /// Returns [`Error::InvalidInput`] when serialization fails.
    pub fn to_json(&self) -> Result<String, Error> {
        serde_json::to_string(self).map_err(|error| {
            Error::InvalidInput(format!("could not serialize molecular crystal: {error}"))
        })
    }
}

/// Generates a molecular crystal.
///
/// # Errors
///
/// Returns [`Error::InvalidInput`] for invalid options, a molecule without
/// atoms, or a seed range that cannot place the molecule without clashes.
pub fn generate_molecular_crystal(
    options: &MolecularCrystalOptions,
) -> Result<MolecularCrystal, Error> {
    validate(options)?;
    let group = SpaceGroup::from_number(options.space_group)?;
    let positions = WyckoffPosition::list(group.hall_number)?;
    let molecule = molecule_geometry(&options.molecule)?;
    let total_atoms = molecule.points.len() * options.copies as usize * 4;
    let volume = options.volume_per_atom * total_atoms.max(1) as f64;
    let mut rng = SplitMix64::new(options.seed);
    let mut last_reason = String::from("no candidate was tried");
    for attempt in 1..=options.max_attempts {
        let lattice = lattice_for_system(crystal_system(group.number), volume, &mut rng)?;
        match place_molecules(options, &group, &positions, &molecule, lattice, &mut rng) {
            Ok(Some(crystal)) => {
                let mut crystal = crystal;
                crystal.attempt = attempt;
                crystal.seed = options.seed;
                return Ok(crystal);
            }
            Ok(None) => {
                last_reason =
                    "the molecule clashed or was incompatible with the site symmetry".to_owned();
            }
            Err(error) => last_reason = error.message(),
        }
    }
    Err(Error::InvalidInput(format!(
        "could not place the molecule in space group {} after {} attempts: {last_reason}",
        options.space_group, options.max_attempts
    )))
}

struct Molecule {
    elements: Vec<String>,
    labels: Vec<String>,
    /// Cartesian coordinates relative to the centroid.
    points: Vec<[f64; 3]>,
}

fn molecule_geometry(structure: &Structure) -> Result<Molecule, Error> {
    let p1 = structure.p1()?;
    if p1.atoms.is_empty() {
        return Err(Error::InvalidInput("the molecule has no atoms".to_owned()));
    }
    let points = p1
        .lattice
        .to_cart_many(&p1.atoms.iter().map(|atom| atom.frac).collect::<Vec<_>>());
    let mut centre = [0.0; 3];
    for point in &points {
        for axis in 0..3 {
            centre[axis] += point[axis] / points.len() as f64;
        }
    }
    Ok(Molecule {
        elements: p1.atoms.iter().map(|atom| atom.element.clone()).collect(),
        labels: p1.atoms.iter().map(|atom| atom.label.clone()).collect(),
        points: points
            .iter()
            .map(|point| {
                [
                    point[0] - centre[0],
                    point[1] - centre[1],
                    point[2] - centre[2],
                ]
            })
            .collect(),
    })
}

fn place_molecules(
    options: &MolecularCrystalOptions,
    group: &SpaceGroup,
    positions: &[WyckoffPosition],
    molecule: &Molecule,
    lattice: xtalforge_core::Lattice,
    rng: &mut SplitMix64,
) -> Result<Option<MolecularCrystal>, Error> {
    let mut letters = BTreeSet::new();
    let mut chosen: Vec<(usize, [f64; 3], [[f64; 3]; 3])> = Vec::new();
    for _ in 0..options.copies {
        let Some(position) = positions.get(rng.below(positions.len())) else {
            return Ok(None);
        };
        let free = [rng.uniform(), rng.uniform(), rng.uniform()];
        let centre = wrap(evaluate_template(&position.coordinates, free)?);
        let rotation = quaternion_matrix(rng.quaternion());
        letters.insert(position.letter);
        chosen.push((position.multiplicity, centre, rotation));
    }
    // Sites of the asymmetric unit, and the centre each one is built around.
    let mut sites: Vec<Site> = Vec::new();
    let mut site_molecule: Vec<usize> = Vec::new();
    for (copy, (_, centre, rotation)) in chosen.iter().enumerate() {
        let centre_cart = lattice.to_cart(*centre);
        for (index, point) in molecule.points.iter().enumerate() {
            let rotated = apply_matrix(*rotation, *point);
            let cart = [
                centre_cart[0] + rotated[0],
                centre_cart[1] + rotated[1],
                centre_cart[2] + rotated[2],
            ];
            let frac = wrap(lattice.to_frac(cart));
            let element = molecule
                .elements
                .get(index)
                .cloned()
                .unwrap_or_else(|| "X".to_owned());
            let label = format!(
                "{}-{copy}",
                molecule
                    .labels
                    .get(index)
                    .cloned()
                    .unwrap_or_else(|| element.clone())
            );
            sites.push(Site::new(label, &element, frac)?);
            site_molecule.push(copy);
        }
    }
    // Expand manually to validate each unique image before committing.
    let mut images: Vec<(usize, Vec<[f64; 3]>)> = Vec::new();
    for (copy, (_, centre, _)) in chosen.iter().enumerate() {
        let base: Vec<[f64; 3]> = sites
            .iter()
            .enumerate()
            .filter(|(index, _)| site_molecule.get(*index) == Some(&copy))
            .map(|(_, site)| site.frac())
            .collect();
        let base_centre = *centre;
        let mut unique: Vec<(usize, [f64; 3], Vec<[f64; 3]>)> = Vec::new();
        for operation in &group.operations {
            let image_centre = wrap(operation.apply(base_centre));
            let image: Vec<[f64; 3]> = base.iter().map(|point| operation.apply(*point)).collect();
            if let Some(existing) = unique.iter().find(|(_, centre, _)| {
                let (_, distance) = lattice.minimum_image([
                    image_centre[0] - centre[0],
                    image_centre[1] - centre[1],
                    image_centre[2] - centre[2],
                ]);
                distance < SAME_POINT
            }) {
                if !same_molecule(&existing.2, &image, molecule, &lattice, options.tolerance) {
                    return Ok(None);
                }
            } else {
                unique.push((copy, image_centre, image));
            }
        }
        for (copy, centre, image) in unique {
            // Shift every image into its own molecule around the unique centre.
            let shift = image[0];
            let _ = shift;
            images.push((copy, image));
            let _ = centre;
        }
    }
    let mut positions_flat = Vec::new();
    let mut groups = Vec::new();
    for (group_index, (_, image)) in images.iter().enumerate() {
        for point in image {
            positions_flat.push(*point);
            groups.push(group_index);
        }
    }
    let distance = minimum_intergroup_distance(&positions_flat, &groups, &lattice);
    if distance < options.min_distance {
        return Ok(None);
    }
    let mut builder = StructureBuilder::new(format!("molecular-{}", group.number), lattice)
        .hall(group.hall_symbol.clone())
        .metadata("generator", json!("molecular"))
        .metadata("seed", json!(options.seed))
        .metadata("space_group", json!(group.number));
    for operation in &group.operations {
        builder = builder.operation(*operation);
    }
    for site in sites {
        builder = builder.site(site);
    }
    let structure = builder.build()?;
    let dataset = detect_symmetry(&structure, SymmetryOptions::default())?;
    if dataset.number != group.number {
        return Ok(None);
    }
    Ok(Some(MolecularCrystal {
        structure,
        space_group: group.number,
        wyckoffs: letters.iter().map(char::to_string).collect(),
        copies: options.copies,
        molecules_per_cell: images.len(),
        seed: options.seed,
        attempt: 1,
    }))
}

/// Whether two lists of same-order atoms describe the same molecule within
/// `tolerance`, allowing for a permutation of the atom order.
fn same_molecule(
    left: &[[f64; 3]],
    right: &[[f64; 3]],
    molecule: &Molecule,
    lattice: &xtalforge_core::Lattice,
    tolerance: f64,
) -> bool {
    if left.len() != right.len() {
        return false;
    }
    let mut used = vec![false; right.len()];
    for (index, first) in left.iter().enumerate() {
        let element = molecule.elements.get(index);
        let mut matched = false;
        for (other, second) in right.iter().enumerate() {
            if used.get(other).copied().unwrap_or(true) {
                continue;
            }
            if molecule.elements.get(other) != element {
                continue;
            }
            let (_, distance) = lattice.minimum_image([
                first[0] - second[0],
                first[1] - second[1],
                first[2] - second[2],
            ]);
            if distance <= tolerance {
                used[other] = true;
                matched = true;
                break;
            }
        }
        if !matched {
            return false;
        }
    }
    true
}

fn validate(options: &MolecularCrystalOptions) -> Result<(), Error> {
    if !(1..=230).contains(&options.space_group) {
        return Err(Error::InvalidInput(format!(
            "space-group number must be in 1..=230, got {}",
            options.space_group
        )));
    }
    if options.copies == 0 {
        return Err(Error::InvalidInput(
            "copies must be at least one".to_owned(),
        ));
    }
    if !options.volume_per_atom.is_finite() || options.volume_per_atom <= 0.0 {
        return Err(Error::InvalidInput(
            "volume per atom must be positive and finite".to_owned(),
        ));
    }
    if !options.min_distance.is_finite() || options.min_distance < 0.0 {
        return Err(Error::InvalidInput(
            "minimum distance must be finite and non-negative".to_owned(),
        ));
    }
    if !options.tolerance.is_finite() || options.tolerance <= 0.0 {
        return Err(Error::InvalidInput(
            "tolerance must be positive and finite".to_owned(),
        ));
    }
    if options.max_attempts == 0 {
        return Err(Error::InvalidInput(
            "max_attempts must be at least one".to_owned(),
        ));
    }
    Ok(())
}
