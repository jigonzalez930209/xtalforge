#![allow(clippy::indexing_slicing, clippy::needless_range_loop)] // fixed-size numeric kernels

//! Crystal surfaces: Wulff construction and adsorption sites (Task 9.4).
//!
//! The Wulff shape is the intersection of the half-spaces whose normals are
//! the reciprocal-lattice directions and whose offsets are the surface
//! energies. Adsorption sites are the top, bridge and hollow positions of the
//! uppermost layer of a slab.

use serde::{Deserialize, Serialize};

use crate::polyhedron::{HalfSpace, Polyhedron};
use crate::{Error, Lattice, Site, Structure, StructureBuilder};

/// Surface energy of one Miller direction.
#[derive(Debug, Clone, Copy, PartialEq, Serialize, Deserialize)]
pub struct SurfaceEnergy {
    /// Miller index.
    pub hkl: [i32; 3],
    /// Surface energy in J/m2; only ratios matter for the shape.
    pub energy: f64,
}

/// One facet of the Wulff shape.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct WulffFacet {
    /// Miller index of the family that produced the facet.
    pub hkl: [i32; 3],
    /// Outward unit normal in Cartesian coordinates.
    pub normal: [f64; 3],
    /// Facet area in square Angstrom.
    pub area: f64,
    /// Fraction of the total surface area.
    pub fraction: f64,
    /// Surface energy of the facet.
    pub energy: f64,
}

/// The Wulff polyhedron.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct WulffShape {
    /// Vertices in Cartesian coordinates.
    pub vertices: Vec<[f64; 3]>,
    /// Facets with their areas and fractions.
    pub facets: Vec<WulffFacet>,
    /// Volume in cubic Angstrom.
    pub volume: f64,
    /// Total surface area in square Angstrom.
    pub surface_area: f64,
}

/// Builds the Wulff shape from a lattice and a list of surface energies.
///
/// Directions are reciprocal-lattice vectors of `hkl`; include every symmetry
/// equivalent direction (the CLI expands them with the point group).
///
/// # Errors
///
/// Returns [`Error::InvalidInput`] for an empty list, a non-positive energy or
/// a degenerate intersection.
pub fn wulff(lattice: &Lattice, energies: &[SurfaceEnergy]) -> Result<WulffShape, Error> {
    if energies.is_empty() {
        return Err(Error::InvalidInput(
            "the Wulff construction needs at least one surface energy".to_owned(),
        ));
    }
    let reciprocal = lattice.reciprocal(false)?;
    let matrix = reciprocal.matrix();
    let mut planes = Vec::with_capacity(energies.len());
    let mut normals = Vec::with_capacity(energies.len());
    for surface in energies {
        if !surface.energy.is_finite() || surface.energy <= 0.0 {
            return Err(Error::InvalidInput(format!(
                "surface energy for {:?} must be positive and finite",
                surface.hkl
            )));
        }
        let [h, k, l] = surface.hkl.map(f64::from);
        let normal = [
            h * matrix[0][0] + k * matrix[1][0] + l * matrix[2][0],
            h * matrix[0][1] + k * matrix[1][1] + l * matrix[2][1],
            h * matrix[0][2] + k * matrix[1][2] + l * matrix[2][2],
        ];
        let length = (normal[0] * normal[0] + normal[1] * normal[1] + normal[2] * normal[2]).sqrt();
        if length < 1e-12 {
            return Err(Error::InvalidInput(format!(
                "direction {:?} has no reciprocal vector",
                surface.hkl
            )));
        }
        planes.push(HalfSpace::new(normal, surface.energy)?);
        normals.push((
            [normal[0] / length, normal[1] / length, normal[2] / length],
            surface.hkl,
            surface.energy,
        ));
    }
    let polyhedron = Polyhedron::from_half_spaces(&planes)?;
    let total = polyhedron.surface_area();
    if total <= 0.0 {
        return Err(Error::InvalidInput(
            "the Wulff shape has no surface".to_owned(),
        ));
    }
    let mut facets = Vec::new();
    for (face, _) in polyhedron.faces.iter().enumerate() {
        let Some(normal) = polyhedron.face_normal(face) else {
            continue;
        };
        let mut best: Option<(f64, [i32; 3], f64)> = None;
        for (candidate, hkl, energy) in &normals {
            let cosine =
                normal[0] * candidate[0] + normal[1] * candidate[1] + normal[2] * candidate[2];
            if cosine > best.map_or(0.99, |(value, _, _)| value) {
                best = Some((cosine, *hkl, *energy));
            }
        }
        let Some((_, hkl, energy)) = best else {
            continue;
        };
        let area = polyhedron.face_area(face);
        facets.push(WulffFacet {
            hkl,
            normal,
            area,
            fraction: area / total,
            energy,
        });
    }
    facets.sort_by(|first, second| {
        second
            .area
            .partial_cmp(&first.area)
            .unwrap_or(std::cmp::Ordering::Equal)
            .then_with(|| first.hkl.cmp(&second.hkl))
    });
    let volume = polyhedron.volume();
    Ok(WulffShape {
        vertices: polyhedron.vertices,
        facets,
        volume,
        surface_area: total,
    })
}

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
    if !tolerance.is_finite() || tolerance <= 0.0 {
        return Err(Error::InvalidInput(
            "adsorption tolerance must be positive and finite".to_owned(),
        ));
    }
    let p1 = slab.p1()?;
    if p1.atoms.is_empty() {
        return Err(Error::InvalidInput("the slab has no atoms".to_owned()));
    }
    let top = p1
        .atoms
        .iter()
        .map(|atom| atom.frac[2])
        .fold(f64::NEG_INFINITY, f64::max);
    let surface: Vec<usize> = p1
        .atoms
        .iter()
        .enumerate()
        .filter(|(_, atom)| (atom.frac[2] - top).abs() <= tolerance)
        .map(|(index, _)| index)
        .collect();
    let mut sites: Vec<AdsorptionSite> = Vec::new();
    for &index in &surface {
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
                bonds.push((first, second, distance));
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
