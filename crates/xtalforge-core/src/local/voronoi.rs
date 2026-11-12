#![allow(clippy::indexing_slicing, clippy::needless_range_loop)] // fixed-size numeric kernels

//! Periodic Voronoi tessellation (Task 9.3.1).
//!
//! Every Voronoi cell is the intersection of the bisector half-spaces of the
//! atoms within a cutoff. The cell is clipped exactly with
//! [`crate::polyhedron`], and each face is assigned to the neighbor that
//! created it, which yields solid-angle weights for coordination analysis.

use serde::{Deserialize, Serialize};

use crate::polyhedron::{HalfSpace, Polyhedron};
use crate::{Error, Structure, neighbors};

/// One face of a Voronoi cell and the neighbor it belongs to.
#[derive(Debug, Clone, Copy, PartialEq, Serialize, Deserialize)]
pub struct VoronoiNeighbor {
    /// P1 index of the neighboring atom.
    pub index: usize,
    /// Lattice image of the neighbor relative to the central atom.
    pub image: [i32; 3],
    /// Distance in Angstrom.
    pub distance: f64,
    /// Face area in square Angstrom.
    pub area: f64,
    /// Solid angle of the face in steradians.
    pub solid_angle: f64,
}

/// The Voronoi cell of one atom.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct VoronoiCell {
    /// P1 index of the central atom.
    pub index: usize,
    /// Cell volume in cubic Angstrom.
    pub volume: f64,
    /// Faces ordered by decreasing solid angle.
    pub neighbors: Vec<VoronoiNeighbor>,
    /// Voronoi vertices relative to the central atom, in Cartesian
    /// coordinates.
    pub vertices: Vec<[f64; 3]>,
}

/// All Voronoi cells of a structure.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct VoronoiReport {
    /// One cell per P1 atom, in atom order.
    pub cells: Vec<VoronoiCell>,
    /// Sum of the cell volumes, which equals the cell volume when the cutoff
    /// captured every relevant face.
    pub total_volume: f64,
}

/// Computes the periodic Voronoi cell of every atom within `cutoff` Angstrom.
///
/// The cutoff must include the faces of the cell. A value of at least twice
/// the longest nearest-neighbor distance is safe and is documented on the
/// algorithm page.
///
/// # Errors
///
/// Returns [`Error::InvalidInput`] for a non-positive cutoff, a structure with
/// fewer than two atoms, or a neighbor whose Voronoi cell is empty.
pub fn voronoi(structure: &Structure, cutoff: f64) -> Result<VoronoiReport, Error> {
    if !cutoff.is_finite() || cutoff <= 0.0 {
        return Err(Error::InvalidInput(format!(
            "Voronoi cutoff must be positive and finite, got {cutoff}"
        )));
    }
    let p1 = structure.p1()?;
    if p1.atoms.len() < 2 {
        return Err(Error::InvalidInput(
            "Voronoi tessellation needs at least two atoms".to_owned(),
        ));
    }
    let matrix = p1.lattice.matrix();
    let cart: Vec<[f64; 3]> = p1
        .atoms
        .iter()
        .map(|atom| p1.lattice.to_cart(atom.frac))
        .collect();
    let pairs = neighbors(&p1_to_structure(&p1)?, cutoff)?;
    let mut adjacency: Vec<Vec<(usize, [i32; 3], f64)>> = vec![Vec::new(); p1.atoms.len()];
    for pair in pairs {
        let Some(first) = adjacency.get_mut(pair.i) else {
            continue;
        };
        first.push((pair.j, pair.image, pair.distance));
        let Some(second) = adjacency.get_mut(pair.j) else {
            continue;
        };
        second.push((
            pair.i,
            [-pair.image[0], -pair.image[1], -pair.image[2]],
            pair.distance,
        ));
    }
    let mut cells = Vec::with_capacity(p1.atoms.len());
    let mut total_volume = 0.0;
    for index in 0..p1.atoms.len() {
        let Some(&centre) = cart.get(index) else {
            continue;
        };
        let Some(entries) = adjacency.get(index) else {
            continue;
        };
        if entries.is_empty() {
            return Err(Error::InvalidInput(format!(
                "atom {index} has no neighbor within {cutoff} Angstrom"
            )));
        }
        let mut planes = Vec::with_capacity(entries.len());
        let mut vectors = Vec::with_capacity(entries.len());
        for &(other, image, distance) in entries {
            let Some(&target) = cart.get(other) else {
                continue;
            };
            let shifted = [
                target[0]
                    + f64::from(image[0]) * matrix[0][0]
                    + f64::from(image[1]) * matrix[1][0]
                    + f64::from(image[2]) * matrix[2][0],
                target[1]
                    + f64::from(image[0]) * matrix[0][1]
                    + f64::from(image[1]) * matrix[1][1]
                    + f64::from(image[2]) * matrix[2][1],
                target[2]
                    + f64::from(image[0]) * matrix[0][2]
                    + f64::from(image[1]) * matrix[1][2]
                    + f64::from(image[2]) * matrix[2][2],
            ];
            let normal = [
                shifted[0] - centre[0],
                shifted[1] - centre[1],
                shifted[2] - centre[2],
            ];
            let plane = HalfSpace::new(normal, 0.5 * distance * distance)?;
            planes.push(plane);
            vectors.push((
                HalfSpace::new(
                    [
                        normal[0] / distance,
                        normal[1] / distance,
                        normal[2] / distance,
                    ],
                    0.5 * distance,
                )?,
                other,
                image,
                distance,
            ));
        }
        let polyhedron = Polyhedron::from_half_spaces(&planes)?;
        let volume = polyhedron.volume();
        total_volume += volume;
        let mut faces = Vec::new();
        for (face, _) in polyhedron.faces.iter().enumerate() {
            let Some(normal) = polyhedron.face_normal(face) else {
                continue;
            };
            let mut best: Option<(f64, usize)> = None;
            for (position, (plane, _, _, _)) in vectors.iter().enumerate() {
                let cosine = normal[0] * plane.normal[0]
                    + normal[1] * plane.normal[1]
                    + normal[2] * plane.normal[2];
                if cosine > best.map_or(0.999, |(value, _)| value) {
                    best = Some((cosine, position));
                }
            }
            let Some((cosine, position)) = best else {
                continue;
            };
            if cosine < 0.9 {
                continue;
            }
            let Some((_, other, image, distance)) = vectors.get(position) else {
                continue;
            };
            faces.push(VoronoiNeighbor {
                index: *other,
                image: *image,
                distance: *distance,
                area: polyhedron.face_area(face),
                solid_angle: polyhedron.face_solid_angle(face, centre),
            });
        }
        faces.sort_by(|first, second| {
            second
                .solid_angle
                .partial_cmp(&first.solid_angle)
                .unwrap_or(std::cmp::Ordering::Equal)
                .then_with(|| (first.index, first.image).cmp(&(second.index, second.image)))
        });
        cells.push(VoronoiCell {
            index,
            volume,
            neighbors: faces,
            vertices: polyhedron.vertices.clone(),
        });
    }
    Ok(VoronoiReport {
        cells,
        total_volume,
    })
}

fn p1_to_structure(p1: &crate::P1Cell) -> Result<Structure, Error> {
    let mut builder = crate::StructureBuilder::new("p1", p1.lattice);
    for atom in &p1.atoms {
        builder = builder.site(
            crate::Site::new(atom.label.clone(), &atom.element, atom.frac)?
                .with_occupancy(atom.occupancy)?,
        );
    }
    builder.build()
}
