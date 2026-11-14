#![allow(clippy::indexing_slicing, clippy::needless_range_loop)] // fixed-size numeric kernels

//! Coordination numbers and continuous shape measures (Task 9.3.2).
//!
//! # CrystalNN-like weighting
//!
//! [`CoordinationStrategy::CrystalNn`] is an independent, pure-Rust
//! approximation of the weighting idea behind pymatgen's `CrystalNN`; it does
//! not call pymatgen and its tests are ordinary Rust integration tests. For
//! every Voronoi face the reported weight is
//!
//! ```text
//! w = exp(-(d - d_min)^2 / (2 * sigma^2)) * Omega / (4 * pi) * c / c_max
//! ```
//!
//! with `sigma = 0.3` Angstrom, `d_min` the shortest neighbor distance in the
//! cell, `Omega = 4 * A / d^2` the planar solid-angle estimate of a face of
//! area `A` at the bisector distance `d / 2`, and `c` the chemical factor
//! below. The factors `c` are divided by the largest one in the shell, so a
//! weight never exceeds its solid-angle fraction. Faces with `w <= 0.05` are
//! discarded and the rest are sorted by weight, distance, index and image,
//! which makes the result deterministic.
//!
//! The chemical factor is
//!
//! ```text
//! c = 1 + 3 * sqrt(|X_center - X_neighbor| / 3.3)
//! ```
//!
//! when the two atoms lie on opposite sides of the arithmetic mean Pauling
//! electronegativity of the P1 cell, and `c = 1` otherwise. Atoms below the
//! mean are treated as cation-like and atoms at or above it as anion-like.
//! Missing electronegativities disable the correction for that pair.
//!
//! pymatgen's `cation_anion` option restricts Voronoi targets to sites with
//! opposite (or zero) formal charge, which needs oxidation states that
//! xtalforge P1 atoms do not carry; the mean-electronegativity classification
//! is therefore an approximation. The chemical factor borrows the form of
//! pymatgen's `x_diff_weight` term. The distance term is also deliberate:
//! pymatgen uses a cosine cutoff based on covalent radii, whereas this model
//! uses a Gaussian around the nearest-neighbor distance.

use serde::{Deserialize, Serialize};

use crate::elements::element_by_symbol;
use crate::local::shapes::{ShapeMeasure, shape_measures};
use crate::local::voronoi::{VoronoiNeighbor, voronoi};
use crate::{Error, P1Atom, Structure};

/// Width of the Gaussian distance term, in Angstrom.
const CRYSTAL_NN_SIGMA: f64 = 0.3;
/// Preference for cation-anion pairs over same-class pairs.
const CRYSTAL_NN_X_DIFF_WEIGHT: f64 = 3.0;
/// Largest Pauling electronegativity difference used by pymatgen.
const CRYSTAL_NN_MAX_X_DIFF: f64 = 3.3;
/// Smallest weight retained in a CrystalNN-like shell.
const CRYSTAL_NN_WEIGHT_MIN: f64 = 0.05;

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
    /// A CrystalNN-like distance, solid-angle and electronegativity weighting.
    ///
    /// See the [module documentation](self) for the exact formula.
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
            let mean = mean_electronegativity(&p1.atoms);
            let centre = p1
                .atoms
                .get(index)
                .and_then(|atom| electronegativity(&atom.element));
            let mut strongest_chemical = 1.0_f64;
            let mut weighted = Vec::with_capacity(cell.neighbors.len());
            for neighbor in &cell.neighbors {
                let neighbor_x = p1
                    .atoms
                    .get(neighbor.index)
                    .and_then(|atom| electronegativity(&atom.element));
                let chemical = cation_anion_factor(centre, neighbor_x, mean);
                strongest_chemical = strongest_chemical.max(chemical);
                let distance_weight = (-((neighbor.distance - minimum).powi(2))
                    / (2.0 * CRYSTAL_NN_SIGMA * CRYSTAL_NN_SIGMA))
                    .exp();
                // The stored Voronoi solid angle is evaluated with an absolute
                // origin, which is only correct for atoms at the cell origin.
                // The planar estimate from the face area and the bisector
                // distance is symmetric for equivalent faces and independent
                // of the central atom position.
                let solid_angle = if neighbor.distance > 0.0 {
                    4.0 * neighbor.area / (neighbor.distance * neighbor.distance)
                } else {
                    neighbor.solid_angle
                };
                let solid_weight = solid_angle / (4.0 * std::f64::consts::PI);
                let mut converted = to_neighbor(neighbor);
                converted.solid_angle = Some(solid_angle);
                converted.weight = distance_weight * solid_weight * chemical;
                weighted.push(converted);
            }
            let list = weighted
                .into_iter()
                .map(|mut neighbor| {
                    neighbor.weight /= strongest_chemical;
                    neighbor
                })
                .filter(|neighbor| neighbor.weight > CRYSTAL_NN_WEIGHT_MIN)
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

/// Pauling electronegativity of an element symbol, when the table defines it.
fn electronegativity(symbol: &str) -> Option<f64> {
    element_by_symbol(symbol)
        .ok()
        .and_then(|element| element.electronegativity)
}

/// Arithmetic mean over the P1 atoms whose electronegativity is defined.
fn mean_electronegativity(atoms: &[P1Atom]) -> Option<f64> {
    let mut total = 0.0;
    let mut count = 0usize;
    for atom in atoms {
        if let Some(value) = electronegativity(&atom.element) {
            total += value;
            count += 1;
        }
    }
    if count == 0 {
        None
    } else {
        Some(total / count as f64)
    }
}

/// Chemical preference of a pair, relative to the cell mean electronegativity.
///
/// Opposite classes (one below the mean, one at or above it) receive the
/// pymatgen-style `1 + 3 * sqrt(delta / 3.3)` boost; same-class pairs and
/// pairs with missing data receive `1`.
fn cation_anion_factor(centre: Option<f64>, neighbor: Option<f64>, mean: Option<f64>) -> f64 {
    let (Some(centre), Some(neighbor), Some(mean)) = (centre, neighbor, mean) else {
        return 1.0;
    };
    let centre_is_cation = centre < mean;
    let neighbor_is_cation = neighbor < mean;
    if centre_is_cation == neighbor_is_cation {
        return 1.0;
    }
    let difference = (centre - neighbor).abs();
    1.0 + CRYSTAL_NN_X_DIFF_WEIGHT * (difference / CRYSTAL_NN_MAX_X_DIFF).sqrt()
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
