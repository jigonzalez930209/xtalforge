#![allow(clippy::indexing_slicing, clippy::needless_range_loop)] // fixed-size numeric kernels

//! Coordination numbers and continuous shape measures (Task 9.3.2).

use serde::{Deserialize, Serialize};

use crate::local::shapes::{ShapeMeasure, shape_measures};
use crate::local::voronoi::{VoronoiNeighbor, voronoi};
use crate::{Error, Structure};

/// How a coordination shell is selected.
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum CoordinationStrategy {
    /// Every neighbor within a distance cutoff.
    Cutoff(f64),
    /// Voronoi faces with a solid angle above the threshold.
    Voronoi {
        /// Cutoff in Angstrom.
        cutoff: f64,
        /// Minimum solid angle in steradians.
        solid_angle_min: f64,
    },
    /// A CrystalNN-like distance and solid-angle weighting.
    CrystalNn {
        /// Cutoff in Angstrom.
        cutoff: f64,
    },
}

/// One neighbor of a coordination shell.
#[derive(Debug, Clone, Copy, PartialEq, Serialize, Deserialize)]
pub struct CoordinationNeighbor {
    /// P1 index.
    pub index: usize,
    /// Lattice image.
    pub image: [i32; 3],
    /// Distance in Angstrom.
    pub distance: f64,
    /// Weight of the neighbor in the shell.
    pub weight: f64,
    /// Voronoi solid angle, when the strategy provides one.
    pub solid_angle: Option<f64>,
}

/// A coordination shell.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct CoordinationShell {
    /// P1 index of the central atom.
    pub index: usize,
    /// Integer coordination number.
    pub number: usize,
    /// Weighted (effective) coordination number.
    pub effective_number: f64,
    /// Neighbors ordered by decreasing weight.
    pub neighbors: Vec<CoordinationNeighbor>,
    /// Continuous shape measures sorted by increasing distance from ideal.
    pub shapes: Vec<ShapeMeasure>,
}

/// Computes the coordination shell and shape measures of one P1 atom.
///
/// # Errors
///
/// Returns [`Error::InvalidInput`] for an out-of-range index, an invalid
/// cutoff or an empty shell, and the errors of [`voronoi`] for the Voronoi
/// strategies.
pub fn coordination(
    structure: &Structure,
    index: usize,
    strategy: CoordinationStrategy,
) -> Result<CoordinationShell, Error> {
    let p1 = structure.p1()?;
    if index >= p1.atoms.len() {
        return Err(Error::InvalidInput(format!(
            "atom index {index} is out of range for {} atoms",
            p1.atoms.len()
        )));
    }
    let matrix = p1.lattice.matrix();
    let (neighbors, require_voronoi) = match strategy {
        CoordinationStrategy::Cutoff(cutoff) => {
            let pairs = crate::neighbors(&p1_structure(&p1)?, cutoff)?;
            let mut list = Vec::new();
            for pair in pairs {
                if pair.i == index {
                    list.push(CoordinationNeighbor {
                        index: pair.j,
                        image: pair.image,
                        distance: pair.distance,
                        weight: 1.0,
                        solid_angle: None,
                    });
                } else if pair.j == index {
                    list.push(CoordinationNeighbor {
                        index: pair.i,
                        image: [-pair.image[0], -pair.image[1], -pair.image[2]],
                        distance: pair.distance,
                        weight: 1.0,
                        solid_angle: None,
                    });
                }
            }
            (list, false)
        }
        CoordinationStrategy::Voronoi {
            cutoff,
            solid_angle_min,
        } => {
            let report = voronoi(structure, cutoff)?;
            let Some(cell) = report.cells.get(index) else {
                return Err(Error::InvalidInput(format!(
                    "no Voronoi cell for atom {index}"
                )));
            };
            let list = cell
                .neighbors
                .iter()
                .filter(|neighbor| neighbor.solid_angle >= solid_angle_min)
                .map(to_neighbor)
                .collect();
            (list, true)
        }
        CoordinationStrategy::CrystalNn { cutoff } => {
            let report = voronoi(structure, cutoff)?;
            let Some(cell) = report.cells.get(index) else {
                return Err(Error::InvalidInput(format!(
                    "no Voronoi cell for atom {index}"
                )));
            };
            let minimum = cell
                .neighbors
                .iter()
                .map(|neighbor| neighbor.distance)
                .fold(f64::INFINITY, f64::min);
            let list = cell
                .neighbors
                .iter()
                .map(|neighbor| {
                    let distance_weight = (-((neighbor.distance - minimum) / 0.4).powi(2)).exp();
                    let solid_weight = (neighbor.solid_angle / std::f64::consts::PI).sqrt();
                    let mut converted = to_neighbor(neighbor);
                    converted.weight = distance_weight * solid_weight;
                    converted
                })
                .filter(|neighbor| neighbor.weight > 0.05)
                .collect();
            (list, true)
        }
    };
    let mut neighbors: Vec<CoordinationNeighbor> = neighbors;
    neighbors.sort_by(|first, second| {
        second
            .weight
            .partial_cmp(&first.weight)
            .unwrap_or(std::cmp::Ordering::Equal)
            .then_with(|| {
                first
                    .distance
                    .partial_cmp(&second.distance)
                    .unwrap_or(std::cmp::Ordering::Equal)
            })
            .then_with(|| (first.index, first.image).cmp(&(second.index, second.image)))
    });
    let effective_number = neighbors.iter().map(|neighbor| neighbor.weight).sum();
    let centre_cart = p1.lattice.to_cart(p1.atoms[index].frac);
    let points: Vec<[f64; 3]> = neighbors
        .iter()
        .map(|neighbor| {
            let Some(target) = p1.atoms.get(neighbor.index) else {
                return centre_cart;
            };
            let base = p1.lattice.to_cart(target.frac);
            [
                base[0]
                    + f64::from(neighbor.image[0]) * matrix[0][0]
                    + f64::from(neighbor.image[1]) * matrix[1][0]
                    + f64::from(neighbor.image[2]) * matrix[2][0]
                    - centre_cart[0],
                base[1]
                    + f64::from(neighbor.image[0]) * matrix[0][1]
                    + f64::from(neighbor.image[1]) * matrix[1][1]
                    + f64::from(neighbor.image[2]) * matrix[2][1]
                    - centre_cart[1],
                base[2]
                    + f64::from(neighbor.image[0]) * matrix[0][2]
                    + f64::from(neighbor.image[1]) * matrix[1][2]
                    + f64::from(neighbor.image[2]) * matrix[2][2]
                    - centre_cart[2],
            ]
        })
        .collect();
    let shapes = if require_voronoi || neighbors.len() >= 2 {
        shape_measures(&points)
    } else {
        Vec::new()
    };
    Ok(CoordinationShell {
        index,
        number: neighbors.len(),
        effective_number,
        neighbors,
        shapes,
    })
}

fn to_neighbor(neighbor: &VoronoiNeighbor) -> CoordinationNeighbor {
    CoordinationNeighbor {
        index: neighbor.index,
        image: neighbor.image,
        distance: neighbor.distance,
        weight: neighbor.solid_angle,
        solid_angle: Some(neighbor.solid_angle),
    }
}

fn p1_structure(p1: &crate::P1Cell) -> Result<Structure, Error> {
    let mut builder = crate::StructureBuilder::new("p1", p1.lattice);
    for atom in &p1.atoms {
        builder = builder.site(
            crate::Site::new(atom.label.clone(), &atom.element, atom.frac)?
                .with_occupancy(atom.occupancy)?,
        );
    }
    builder.build()
}
