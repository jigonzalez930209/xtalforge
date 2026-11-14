//! Adsorption sites on a slab (Tasks 9.4.4 and 17.3.2).
//!
//! The default search covers the uppermost atomic layer. With
//! [`TerminationOptions::site_layers`] it also enumerates the layers below it.

use serde::{Deserialize, Serialize};

use crate::slab::TerminationOptions;
use crate::{Error, Site, Structure, StructureBuilder};

/// Type of an adsorption site.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
#[serde(rename_all = "snake_case")]
pub enum AdsorptionKind {
    /// On top of one surface atom.
    Top,
    /// Above the midpoint of two bonded surface atoms.
    Bridge,
    /// Above the centroid of three or four bonded surface atoms.
    Hollow,
}

/// One adsorption site on a slab.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct AdsorptionSite {
    /// Site type.
    pub kind: AdsorptionKind,
    /// Cartesian position on the surface.
    pub position: [f64; 3],
    /// Fractional coordinates in the slab cell.
    pub frac: [f64; 3],
    /// P1 indices of the surface atoms that define the site.
    pub atoms: Vec<usize>,
}

/// Enumerates top, bridge and hollow sites of the uppermost slab layer.
///
/// # Errors
///
/// Returns [`Error::InvalidInput`] for an empty slab or a non-positive
/// tolerance, or the errors of the P1 expansion.
pub fn adsorption_sites(slab: &Structure, tolerance: f64) -> Result<Vec<AdsorptionSite>, Error> {
    adsorption_sites_with_options(slab, tolerance, &TerminationOptions::default())
}

/// Enumerates adsorption sites down to the requested number of layers.
///
/// [`TerminationOptions::site_layers`] selects how many atomic layers below
/// the top are searched; the default of `1` reproduces [`adsorption_sites`].
/// Every site sits at the height of the layer that defines it, so lower-layer
/// sites keep the fractional `z` of their own layer.
///
/// # Errors
///
/// Returns [`Error::InvalidInput`] for an empty slab, a non-positive
/// tolerance or invalid options, or the errors of the P1 expansion.
pub fn adsorption_sites_with_options(
    slab: &Structure,
    tolerance: f64,
    options: &TerminationOptions,
) -> Result<Vec<AdsorptionSite>, Error> {
    if !tolerance.is_finite() || tolerance <= 0.0 {
        return Err(Error::InvalidInput(
            "adsorption tolerance must be positive and finite".to_owned(),
        ));
    }
    crate::slab::validate_termination_options(options)?;
    let p1 = slab.p1()?;
    if p1.atoms.is_empty() {
        return Err(Error::InvalidInput("the slab has no atoms".to_owned()));
    }
    let mut sites: Vec<AdsorptionSite> = Vec::new();
    let mut remaining: Vec<usize> = (0..p1.atoms.len()).collect();
    for _ in 0..options.site_layers {
        let Some(top) = remaining
            .iter()
            .filter_map(|index| p1.atoms.get(*index))
            .map(|atom| atom.frac[2])
            .reduce(f64::max)
        else {
            break;
        };
        let mut surface: Vec<usize> = Vec::new();
        let mut rest: Vec<usize> = Vec::new();
        for index in remaining {
            let inside = p1
                .atoms
                .get(index)
                .is_some_and(|atom| (atom.frac[2] - top).abs() <= tolerance);
            if inside {
                surface.push(index);
            } else {
                rest.push(index);
            }
        }
        layer_sites(&p1, top, &surface, &mut sites)?;
        remaining = rest;
    }
    deduplicate(&mut sites, tolerance);
    Ok(sites)
}

/// Places a molecule above an adsorption site.
///
/// The molecule centroid is moved to `site.position` plus `offset` Angstrom
/// along the slab normal (the `c` axis of a slab cell), and its atoms are
/// appended in fractional coordinates of the slab.
///
/// # Errors
///
/// Returns [`Error::InvalidInput`] for an empty molecule or a non-finite
/// offset, or a builder error.
pub fn place_molecule(
    slab: &Structure,
    molecule: &Structure,
    site: &AdsorptionSite,
    offset: f64,
) -> Result<Structure, Error> {
    if !offset.is_finite() {
        return Err(Error::InvalidInput(
            "the adsorption offset must be finite".to_owned(),
        ));
    }
    let molecule_p1 = molecule.p1()?;
    if molecule_p1.atoms.is_empty() {
        return Err(Error::InvalidInput("the molecule has no atoms".to_owned()));
    }
    let normal = slab.lattice().plane_normal([0.0, 0.0, 1.0])?;
    let target = [
        site.position[0] + normal[0] * offset,
        site.position[1] + normal[1] * offset,
        site.position[2] + normal[2] * offset,
    ];
    let mut centroid = [0.0; 3];
    for atom in &molecule_p1.atoms {
        let cart = molecule_p1.lattice.to_cart(atom.frac);
        centroid[0] += cart[0];
        centroid[1] += cart[1];
        centroid[2] += cart[2];
    }
    let count = molecule_p1.atoms.len() as f64;
    centroid = [
        centroid[0] / count,
        centroid[1] / count,
        centroid[2] / count,
    ];
    let mut builder = StructureBuilder::new(slab.name(), *slab.lattice()).hall(slab.hall());
    for operation in slab.operations() {
        builder = builder.operation(*operation);
    }
    for site in slab.sites() {
        builder = builder.site(site.clone());
    }
    for atom in &molecule_p1.atoms {
        let cart = molecule_p1.lattice.to_cart(atom.frac);
        let moved = [
            cart[0] - centroid[0] + target[0],
            cart[1] - centroid[1] + target[1],
            cart[2] - centroid[2] + target[2],
        ];
        let frac = slab.lattice().to_frac(moved);
        builder = builder.site(
            Site::new(atom.label.clone(), &atom.element, frac)?.with_occupancy(atom.occupancy)?,
        );
    }
    for (key, value) in slab.metadata() {
        builder = builder.metadata(key.clone(), value.clone());
    }
    builder.build()
}

/// Adds the top, bridge and hollow sites of one atomic layer.
fn layer_sites(
    p1: &crate::P1Cell,
    top: f64,
    surface: &[usize],
    sites: &mut Vec<AdsorptionSite>,
) -> Result<(), Error> {
    for &index in surface {
        let Some(atom) = p1.atoms.get(index) else {
            continue;
        };
        let frac = [atom.frac[0], atom.frac[1], top];
        sites.push(AdsorptionSite {
            kind: AdsorptionKind::Top,
            position: p1.lattice.to_cart(frac),
            frac,
            atoms: vec![index],
        });
    }
    let mut bonds = Vec::new();
    for (position, &first) in surface.iter().enumerate() {
        for &second in surface.iter().skip(position + 1) {
            let (Some(a), Some(b)) = (p1.atoms.get(first), p1.atoms.get(second)) else {
                continue;
            };
            let distance = p1.lattice.distance(a.frac, b.frac);
            let radius = crate::elements::element_by_symbol(&a.element)?.covalent_radius
                + crate::elements::element_by_symbol(&b.element)?.covalent_radius;
            if distance <= 1.25 * radius {
                bonds.push((first, second));
                let frac = [
                    (a.frac[0] + b.frac[0]) / 2.0,
                    (a.frac[1] + b.frac[1]) / 2.0,
                    top,
                ];
                sites.push(AdsorptionSite {
                    kind: AdsorptionKind::Bridge,
                    position: p1.lattice.to_cart(frac),
                    frac,
                    atoms: vec![first, second],
                });
            }
        }
    }
    for (position, &first) in surface.iter().enumerate() {
        for (other, &second) in surface.iter().enumerate().skip(position + 1) {
            for &third in surface.iter().skip(other + 1) {
                let has = |a: usize, b: usize| {
                    bonds.iter().any(|entry| {
                        (entry.0 == a && entry.1 == b) || (entry.0 == b && entry.1 == a)
                    })
                };
                if !(has(first, second) && has(second, third) && has(first, third)) {
                    continue;
                }
                let (Some(a), Some(b), Some(c)) = (
                    p1.atoms.get(first),
                    p1.atoms.get(second),
                    p1.atoms.get(third),
                ) else {
                    continue;
                };
                let frac = [
                    (a.frac[0] + b.frac[0] + c.frac[0]) / 3.0,
                    (a.frac[1] + b.frac[1] + c.frac[1]) / 3.0,
                    top,
                ];
                sites.push(AdsorptionSite {
                    kind: AdsorptionKind::Hollow,
                    position: p1.lattice.to_cart(frac),
                    frac,
                    atoms: vec![first, second, third],
                });
            }
        }
    }
    Ok(())
}

fn deduplicate(sites: &mut Vec<AdsorptionSite>, tolerance: f64) {
    let mut kept: Vec<AdsorptionSite> = Vec::new();
    for site in sites.drain(..) {
        let duplicate = kept.iter().any(|existing| {
            existing.kind == site.kind
                && existing.atoms.iter().all(|atom| site.atoms.contains(atom))
                && distance(existing.position, site.position) <= tolerance
        });
        if !duplicate {
            kept.push(site);
        }
    }
    kept.sort_by(|first, second| {
        kind_order(first.kind)
            .cmp(&kind_order(second.kind))
            .then_with(|| {
                first
                    .position
                    .partial_cmp(&second.position)
                    .unwrap_or(std::cmp::Ordering::Equal)
            })
    });
    *sites = kept;
}

fn kind_order(kind: AdsorptionKind) -> u8 {
    match kind {
        AdsorptionKind::Top => 0,
        AdsorptionKind::Bridge => 1,
        AdsorptionKind::Hollow => 2,
    }
}

fn distance(first: [f64; 3], second: [f64; 3]) -> f64 {
    let dx = first[0] - second[0];
    let dy = first[1] - second[1];
    let dz = first[2] - second[2];
    (dx * dx + dy * dy + dz * dz).sqrt()
}
