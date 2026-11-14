#![allow(clippy::indexing_slicing, clippy::needless_range_loop)] // fixed-size numeric kernels

//! Crystal surfaces: Wulff construction and adsorption sites (Task 9.4).
//!
//! The Wulff shape is the intersection of the half-spaces whose normals are
//! the reciprocal-lattice directions and whose offsets are the surface
//! energies. Adsorption sites are the top, bridge and hollow positions of the
//! uppermost layer of a slab, or of several layers with
//! [`TerminationOptions::site_layers`](crate::slab::TerminationOptions::site_layers).

use serde::{Deserialize, Serialize};

use crate::polyhedron::{HalfSpace, Polyhedron};
use crate::{Error, Lattice};

mod sites;

pub use sites::{
    AdsorptionKind, AdsorptionSite, adsorption_sites, adsorption_sites_with_options, place_molecule,
};

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
