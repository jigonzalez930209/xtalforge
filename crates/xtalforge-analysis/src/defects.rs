#![allow(clippy::indexing_slicing, clippy::needless_range_loop)] // fixed-size numeric kernels

//! Point defects: symmetry-unique vacancies, substitutions and interstitials
//! (Task 9.5.1).

use serde::{Deserialize, Serialize};

use xtalforge_core::{Error, Lattice, Structure, add_sites, remove_sites, replace_site, voronoi};
use xtalforge_symmetry::{SymmetryOptions, detect_symmetry};

/// Kind of a point defect.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
#[serde(rename_all = "snake_case")]
pub enum DefectKind {
    /// One atom removed.
    Vacancy,
    /// One atom replaced by another species.
    Substitution,
    /// One atom inserted in a void.
    Interstitial,
}

/// One generated defect.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct PointDefect {
    /// Defect type.
    pub kind: DefectKind,
    /// P1 index in the input cell, absent for interstitials.
    pub site_index: Option<usize>,
    /// Original species, absent for interstitials.
    pub original: Option<String>,
    /// Species of the defect.
    pub species: String,
    /// Fractional position in the input cell.
    pub frac: [f64; 3],
    /// Supercell repeats used for the minimum-image distance.
    pub supercell: [u32; 3],
    /// Smallest distance between defect images in Angstrom.
    pub image_distance: f64,
    /// The defect structure.
    pub structure: Structure,
}

/// Symmetry-unique vacancies.
///
/// # Errors
///
/// Returns [`Error::InvalidInput`] for an empty structure or a bad distance,
/// or the errors of the symmetry search and the supercell construction.
pub fn vacancies(structure: &Structure, min_image: f64) -> Result<Vec<PointDefect>, Error> {
    let p1 = structure.p1()?;
    if p1.atoms.is_empty() {
        return Err(Error::InvalidInput(
            "a vacancy needs at least one atom".to_owned(),
        ));
    }
    let representatives = unique_sites(&p1)?;
    let mut defects = Vec::new();
    for index in representatives {
        let Some(atom) = p1.atoms.get(index) else {
            continue;
        };
        let removed = remove_sites(&p1_structure(&p1)?, &[index])?;
        let (expanded, repeats, distance) = expand_to_min_image(&removed, min_image)?;
        defects.push(PointDefect {
            kind: DefectKind::Vacancy,
            site_index: Some(index),
            original: Some(atom.element.clone()),
            species: "vacancy".to_owned(),
            frac: atom.frac,
            supercell: repeats,
            image_distance: distance,
            structure: expanded,
        });
    }
    Ok(defects)
}

/// Symmetry-unique substitutions.
///
/// # Errors
///
/// Returns [`Error::InvalidInput`] for an empty structure, a bad distance or
/// an unknown element, or the errors of the symmetry search.
pub fn substitutions(
    structure: &Structure,
    species: &str,
    min_image: f64,
) -> Result<Vec<PointDefect>, Error> {
    xtalforge_core::elements::element_by_symbol(species)?;
    let p1 = structure.p1()?;
    if p1.atoms.is_empty() {
        return Err(Error::InvalidInput(
            "a substitution needs at least one atom".to_owned(),
        ));
    }
    let representatives = unique_sites(&p1)?;
    let mut defects = Vec::new();
    for index in representatives {
        let Some(atom) = p1.atoms.get(index) else {
            continue;
        };
        if atom.element == species {
            continue;
        }
        let replaced = replace_site(&p1_structure(&p1)?, index, species)?;
        let (expanded, repeats, distance) = expand_to_min_image(&replaced, min_image)?;
        defects.push(PointDefect {
            kind: DefectKind::Substitution,
            site_index: Some(index),
            original: Some(atom.element.clone()),
            species: species.to_owned(),
            frac: atom.frac,
            supercell: repeats,
            image_distance: distance,
            structure: expanded,
        });
    }
    Ok(defects)
}

/// Symmetry-unique interstitial candidates from Voronoi vertices.
///
/// `cutoff` is the Voronoi neighbor cutoff; `minimum_distance` rejects a
/// candidate closer than that to an existing atom.
///
/// # Errors
///
/// Returns [`Error::InvalidInput`] for a bad distance or cutoff, an unknown
/// element, or the errors of the Voronoi tessellation.
pub fn interstitials(
    structure: &Structure,
    species: &str,
    min_image: f64,
    cutoff: f64,
    minimum_distance: f64,
) -> Result<Vec<PointDefect>, Error> {
    xtalforge_core::elements::element_by_symbol(species)?;
    if !minimum_distance.is_finite() || minimum_distance <= 0.0 {
        return Err(Error::InvalidInput(
            "the interstitial clearance must be positive and finite".to_owned(),
        ));
    }
    let p1 = structure.p1()?;
    if p1.atoms.is_empty() {
        return Err(Error::InvalidInput(
            "an interstitial needs at least one atom".to_owned(),
        ));
    }
    let report = voronoi(structure, cutoff)?;
    let operations = detect_symmetry(
        structure,
        SymmetryOptions {
            symprec: 0.05,
            ..SymmetryOptions::default()
        },
    )
    .map(|dataset| dataset.operations)
    .unwrap_or_default();
    let mut candidates: Vec<[f64; 3]> = Vec::new();
    for cell in &report.cells {
        let Some(atom) = p1.atoms.get(cell.index) else {
            continue;
        };
        let centre = p1.lattice.to_cart(atom.frac);
        for vertex in &cell.vertices {
            let frac = Lattice::wrap(p1.lattice.to_frac([
                centre[0] + vertex[0],
                centre[1] + vertex[1],
                centre[2] + vertex[2],
            ]));
            let clearance = p1
                .atoms
                .iter()
                .map(|other| p1.lattice.distance(other.frac, frac))
                .fold(f64::INFINITY, f64::min);
            if clearance < minimum_distance {
                continue;
            }
            if candidates
                .iter()
                .any(|existing| p1.lattice.distance(*existing, frac) < 0.5)
            {
                continue;
            }
            let equivalent = operations.iter().any(|operation| {
                let image = operation.apply(frac);
                candidates
                    .iter()
                    .any(|existing| p1.lattice.distance(*existing, image) < 0.5)
            });
            if equivalent {
                continue;
            }
            candidates.push(frac);
        }
    }
    let mut defects = Vec::new();
    for frac in candidates {
        let added = add_sites(&p1_structure(&p1)?, &[(species.to_owned(), frac)])?;
        let (expanded, repeats, distance) = expand_to_min_image(&added, min_image)?;
        defects.push(PointDefect {
            kind: DefectKind::Interstitial,
            site_index: None,
            original: None,
            species: species.to_owned(),
            frac,
            supercell: repeats,
            image_distance: distance,
            structure: expanded,
        });
    }
    Ok(defects)
}

fn unique_sites(p1: &xtalforge_core::P1Cell) -> Result<Vec<usize>, Error> {
    let structure = p1_structure(p1)?;
    let dataset = detect_symmetry(
        &structure,
        SymmetryOptions {
            symprec: 0.05,
            ..SymmetryOptions::default()
        },
    );
    match dataset {
        Ok(dataset) => Ok(dataset
            .orbits
            .iter()
            .enumerate()
            .filter(|(index, representative)| **representative == *index)
            .map(|(index, _)| index)
            .collect()),
        Err(_) => Ok((0..p1.atoms.len()).collect()),
    }
}

fn expand_to_min_image(
    structure: &Structure,
    min_image: f64,
) -> Result<(Structure, [u32; 3], f64), Error> {
    if !min_image.is_finite() || min_image <= 0.0 {
        return Err(Error::InvalidInput(
            "the minimum image distance must be positive and finite".to_owned(),
        ));
    }
    let spacing = [
        structure.lattice().d_spacing([1.0, 0.0, 0.0])?,
        structure.lattice().d_spacing([0.0, 1.0, 0.0])?,
        structure.lattice().d_spacing([0.0, 0.0, 1.0])?,
    ];
    let mut repeats = [1u32; 3];
    for (index, width) in spacing.iter().enumerate() {
        let count = (min_image / width.max(1e-9)).ceil().max(1.0) as u32;
        if let Some(slot) = repeats.get_mut(index) {
            *slot = count;
        }
    }
    let expanded = structure.supercell(repeats[0], repeats[1], repeats[2])?;
    let distance = repeats
        .iter()
        .zip(spacing)
        .map(|(repeat, width)| f64::from(*repeat) * width)
        .fold(f64::INFINITY, f64::min);
    Ok((expanded, repeats, distance))
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
