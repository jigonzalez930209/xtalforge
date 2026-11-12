#![allow(clippy::indexing_slicing, clippy::needless_range_loop)] // fixed-size numeric kernels

//! Structure matcher and equality classes (Task 9.1).
//!
//! Two structures match when a unimodular lattice mapping, an origin shift and
//! a species-preserving atom permutation bring every site within `stol`. The
//! search works on Niggli-reduced cells, so it does not depend on the input
//! basis. A supercell relationship is accepted when the volume ratio is an
//! integer and `allow_supercell` is set.

use std::collections::BTreeMap;

use serde::{Deserialize, Serialize};
use xtalforge_core::{Error, Lattice, Site, Structure, StructureBuilder};

use crate::matcher_lattice::find_integer_mappings;
use crate::reduction::{Reduction, reduce_structure};

/// How species are compared.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
#[serde(rename_all = "snake_case")]
pub enum SpeciesComparator {
    /// Exact element symbols.
    Exact,
    /// Canonical elements (isotopes collapse to their element).
    Element,
    /// Ignore species and compare geometry only.
    Ignore,
}

/// Options for [`match_structures`].
#[derive(Debug, Clone, Copy, PartialEq, Serialize, Deserialize)]
pub struct MatcherOptions {
    /// Fractional length tolerance for lattice vectors.
    pub ltol: f64,
    /// Site tolerance in Angstrom.
    pub stol: f64,
    /// Angle tolerance in degrees.
    pub angle_tol: f64,
    /// Species comparator.
    pub comparator: SpeciesComparator,
    /// Whether a supercell relationship is accepted.
    pub allow_supercell: bool,
}

impl Default for MatcherOptions {
    fn default() -> Self {
        Self {
            ltol: 0.2,
            stol: 0.3,
            angle_tol: 5.0,
            comparator: SpeciesComparator::Exact,
            allow_supercell: true,
        }
    }
}

/// The result of matching two structures.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct StructureMatch {
    /// Whether the structures match within the tolerances.
    pub matched: bool,
    /// Root mean square displacement in Angstrom.
    pub rms: f64,
    /// Mapping from first atoms to second atoms.
    pub mapping: Vec<Option<usize>>,
    /// Fractional origin shift applied to the second structure.
    pub translation: [f64; 3],
    /// Lattice mapping `new = M * second`.
    pub lattice_mapping: [[i32; 3]; 3],
    /// Supercell matrix applied to the second structure, when any.
    pub supercell: Option<[[i32; 3]; 3]>,
    /// Number of compared atoms.
    pub n_atoms: usize,
}

/// Decides whether two structures are the same up to basis, origin, mapping
/// and tolerance.
///
/// # Errors
///
/// Returns [`Error::InvalidInput`] for invalid tolerances, or the errors of
/// the primitive-cell reduction.
pub fn match_structures(
    first: &Structure,
    second: &Structure,
    options: &MatcherOptions,
) -> Result<StructureMatch, Error> {
    validate(options)?;
    let first = prepare(first, options)?;
    let mut second = prepare(second, options)?;
    if first.atoms.is_empty() || second.atoms.is_empty() {
        return Err(Error::InvalidInput(
            "the matcher needs at least one atom per structure".to_owned(),
        ));
    }
    if first.atoms.len() < second.atoms.len() && options.allow_supercell {
        let first_structure = p1_structure(&first)?;
        let second_structure = p1_structure(&second)?;
        let swapped = match_structures(&second_structure, &first_structure, options)?;
        if swapped.matched {
            return Ok(StructureMatch {
                matched: true,
                rms: swapped.rms,
                mapping: Vec::new(),
                translation: [0.0; 3],
                lattice_mapping: swapped.lattice_mapping,
                supercell: swapped.supercell,
                n_atoms: first.atoms.len(),
            });
        }
        return Ok(no_match());
    }
    let ratio = first.lattice.volume() / second.lattice.volume();
    let mut supercell = None;
    if first.atoms.len() != second.atoms.len() {
        if !options.allow_supercell {
            return Ok(no_match());
        }
        let count_ratio = first.atoms.len() as f64 / second.atoms.len() as f64;
        if count_ratio < 1.0 || (count_ratio - count_ratio.round()).abs() > 1e-6 {
            return Ok(no_match());
        }
        let factor = count_ratio.round() as u32;
        if factor == 0 || factor > 64 || (ratio - f64::from(factor)).abs() > 0.5 * f64::from(factor)
        {
            return Ok(no_match());
        }
        let mappings =
            find_integer_mappings(&first.lattice, &second.lattice, options, f64::from(factor))?;
        let Some(matrix) = mappings.first() else {
            return Ok(no_match());
        };
        let integer = matrix.map(|row| row.map(|value| value.round() as i32));
        let expanded = expand(&second, integer)?;
        if expanded.atoms.len() != first.atoms.len() {
            return Ok(no_match());
        }
        second = expanded;
        supercell = Some(integer);
    }
    let mappings = find_integer_mappings(&first.lattice, &second.lattice, options, 1.0)?;
    let mut best: Option<BestMatch> = None;
    for matrix in &mappings {
        let integer = matrix.map(|row| row.map(|value| value.round() as i32));
        let mapped = map_cell(&second, &integer)?;
        let Some((translation, mapping, rms)) = match_sites(&first, &mapped, options) else {
            continue;
        };
        if best
            .as_ref()
            .is_none_or(|(current, _, _, _)| rms < *current)
        {
            best = Some((rms, translation, mapping, integer));
        }
    }
    let Some((rms, translation, mapping, integer)) = best else {
        return Ok(no_match());
    };
    Ok(StructureMatch {
        matched: true,
        rms,
        mapping: mapping.into_iter().map(Some).collect(),
        translation,
        lattice_mapping: integer,
        supercell,
        n_atoms: first.atoms.len(),
    })
}

/// Groups structures into deterministic equality classes.
///
/// Classes are ordered by their first member; each structure joins the first
/// class whose representative matches it.
///
/// # Errors
///
/// Returns the errors of [`match_structures`].
pub fn group_structures(
    structures: &[Structure],
    options: &MatcherOptions,
) -> Result<Vec<Vec<usize>>, Error> {
    let mut groups: Vec<Vec<usize>> = Vec::new();
    let mut representatives: Vec<usize> = Vec::new();
    for (index, structure) in structures.iter().enumerate() {
        let mut found = None;
        for (group, representative) in representatives.iter().enumerate() {
            let Some(candidate) = structures.get(*representative) else {
                continue;
            };
            if match_structures(candidate, structure, options)?.matched {
                found = Some(group);
                break;
            }
        }
        match found {
            Some(group) => {
                if let Some(entries) = groups.get_mut(group) {
                    entries.push(index);
                }
            }
            None => {
                groups.push(vec![index]);
                representatives.push(index);
            }
        }
    }
    Ok(groups)
}

/// Interpolates between two matched structures.
///
/// # Errors
///
/// Returns the errors of [`match_structures`] when they do not match, or the
/// errors of the interpolation.
pub fn interpolate_matched(
    first: &Structure,
    second: &Structure,
    images: u32,
    options: &MatcherOptions,
) -> Result<Vec<Structure>, Error> {
    let outcome = match_structures(first, second, options)?;
    if !outcome.matched {
        return Err(Error::InvalidInput(
            "the structures do not match; interpolation needs a one-to-one mapping".to_owned(),
        ));
    }
    let start = prepare(first, options)?;
    let end = prepare(second, options)?;
    let mut ordered = Vec::with_capacity(start.atoms.len());
    for (index, atom) in start.atoms.iter().enumerate() {
        let Some(Some(target)) = outcome.mapping.get(index) else {
            return Err(Error::InvalidInput(
                "the match mapping is incomplete".to_owned(),
            ));
        };
        let Some(other) = end.atoms.get(*target) else {
            return Err(Error::InvalidInput(
                "the match mapping is out of range".to_owned(),
            ));
        };
        let cart = end.lattice.to_cart(other.frac);
        let frac = start.lattice.to_frac(cart);
        ordered.push(Site::new(
            atom.label.clone(),
            &atom.element,
            Lattice::wrap(frac),
        )?);
    }
    let mut builder = StructureBuilder::new(second.name(), start.lattice);
    for site in ordered {
        builder = builder.site(site);
    }
    let aligned = builder.build()?;
    let first_p1 = p1_structure(&start)?;
    xtalforge_core::interpolate(&first_p1, &aligned, images)
}

fn prepare(
    structure: &Structure,
    options: &MatcherOptions,
) -> Result<xtalforge_core::P1Cell, Error> {
    let _ = options;
    let reduced = reduce_structure(structure, Reduction::Niggli)?;
    reduced.structure.p1()
}

type BestMatch = (f64, [f64; 3], Vec<usize>, [[i32; 3]; 3]);

fn no_match() -> StructureMatch {
    StructureMatch {
        matched: false,
        rms: f64::INFINITY,
        mapping: Vec::new(),
        translation: [0.0; 3],
        lattice_mapping: [[1, 0, 0], [0, 1, 0], [0, 0, 1]],
        supercell: None,
        n_atoms: 0,
    }
}

fn validate(options: &MatcherOptions) -> Result<(), Error> {
    if !options.ltol.is_finite() || options.ltol <= 0.0 {
        return Err(Error::InvalidInput(
            "ltol must be positive and finite".to_owned(),
        ));
    }
    if !options.stol.is_finite() || options.stol <= 0.0 {
        return Err(Error::InvalidInput(
            "stol must be positive and finite".to_owned(),
        ));
    }
    if !options.angle_tol.is_finite() || options.angle_tol <= 0.0 {
        return Err(Error::InvalidInput(
            "angle_tol must be positive and finite".to_owned(),
        ));
    }
    Ok(())
}

fn expand(
    structure: &xtalforge_core::P1Cell,
    matrix: [[i32; 3]; 3],
) -> Result<xtalforge_core::P1Cell, Error> {
    let base = p1_structure(structure)?;
    let expanded = base.supercell_matrix(matrix.map(|row| row.map(f64::from)))?;
    expanded.p1()
}

fn map_cell(
    structure: &xtalforge_core::P1Cell,
    matrix: &[[i32; 3]; 3],
) -> Result<xtalforge_core::P1Cell, Error> {
    let rows = matrix.map(|row| row.map(f64::from));
    let lattice = Lattice::new(rows)?.transform(structure.lattice.matrix())?;
    let inverse = lattice.inverse();
    let site_lattice = lattice;
    let mut builder = StructureBuilder::new("mapped", site_lattice);
    for atom in &structure.atoms {
        let cart = structure.lattice.to_cart(atom.frac);
        let frac = [
            cart[0] * inverse[0][0] + cart[1] * inverse[1][0] + cart[2] * inverse[2][0],
            cart[0] * inverse[0][1] + cart[1] * inverse[1][1] + cart[2] * inverse[2][1],
            cart[0] * inverse[0][2] + cart[1] * inverse[1][2] + cart[2] * inverse[2][2],
        ];
        builder = builder.site(Site::new(
            atom.label.clone(),
            &atom.element,
            Lattice::wrap(frac),
        )?);
    }
    builder.build()?.p1()
}

fn match_sites(
    first: &xtalforge_core::P1Cell,
    second: &xtalforge_core::P1Cell,
    options: &MatcherOptions,
) -> Option<([f64; 3], Vec<usize>, f64)> {
    let mut best: Option<([f64; 3], Vec<usize>, f64)> = None;
    let anchors = first.atoms.len().min(4);
    for anchor in 0..anchors {
        let Some(atom) = first.atoms.get(anchor) else {
            continue;
        };
        for candidate in &second.atoms {
            if !compatible(&atom.element, &candidate.element, options.comparator) {
                continue;
            }
            let translation = [
                atom.frac[0] - candidate.frac[0],
                atom.frac[1] - candidate.frac[1],
                atom.frac[2] - candidate.frac[2],
            ];
            let Some((mapping, rms)) = assign(first, second, translation, options) else {
                continue;
            };
            if best.as_ref().is_none_or(|(_, _, current)| rms < *current) {
                best = Some((translation, mapping, rms));
            }
        }
    }
    best
}

fn assign(
    first: &xtalforge_core::P1Cell,
    second: &xtalforge_core::P1Cell,
    translation: [f64; 3],
    options: &MatcherOptions,
) -> Option<(Vec<usize>, f64)> {
    let mut used = vec![false; second.atoms.len()];
    let mut mapping = Vec::with_capacity(first.atoms.len());
    let mut squared = 0.0;
    for atom in &first.atoms {
        let target_frac = atom.frac;
        let mut best: Option<(f64, usize)> = None;
        for (index, candidate) in second.atoms.iter().enumerate() {
            if used.get(index).copied().unwrap_or(true) {
                continue;
            }
            if !compatible(&atom.element, &candidate.element, options.comparator) {
                continue;
            }
            let delta = [
                target_frac[0] - (candidate.frac[0] + translation[0]),
                target_frac[1] - (candidate.frac[1] + translation[1]),
                target_frac[2] - (candidate.frac[2] + translation[2]),
            ];
            let (_, distance) = first.lattice.minimum_image(delta);
            if distance > options.stol {
                continue;
            }
            if best.is_none_or(|(current, _)| distance < current) {
                best = Some((distance, index));
            }
        }
        let (distance, index) = best?;
        if let Some(slot) = used.get_mut(index) {
            *slot = true;
        }
        mapping.push(index);
        squared += distance * distance;
    }
    let rms = (squared / mapping.len() as f64).sqrt();
    Some((mapping, rms))
}

fn compatible(first: &str, second: &str, comparator: SpeciesComparator) -> bool {
    match comparator {
        SpeciesComparator::Exact => first == second,
        SpeciesComparator::Element => canonical(first) == canonical(second),
        SpeciesComparator::Ignore => true,
    }
}

fn canonical(element: &str) -> String {
    xtalforge_core::elements::element_by_symbol(element)
        .map(|info| info.symbol.clone())
        .unwrap_or_else(|_| element.to_owned())
}

fn p1_structure(cell: &xtalforge_core::P1Cell) -> Result<Structure, Error> {
    let mut builder = StructureBuilder::new("p1", cell.lattice);
    for atom in &cell.atoms {
        builder = builder.site(
            Site::new(atom.label.clone(), &atom.element, atom.frac)?
                .with_occupancy(atom.occupancy)?,
        );
    }
    builder.build()
}

/// Convenience wrapper returning the class of each structure.
///
/// # Errors
///
/// Returns the errors of [`group_structures`].
pub fn classify_structures(
    structures: &[Structure],
    options: &MatcherOptions,
) -> Result<BTreeMap<usize, Vec<usize>>, Error> {
    let groups = group_structures(structures, options)?;
    Ok(groups.into_iter().enumerate().collect())
}
