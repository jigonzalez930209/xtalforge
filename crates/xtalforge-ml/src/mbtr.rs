//! MBTR descriptors (Task 21.1.3).
//!
//! The many-body tensor uses Gaussian broadening on a fixed grid:
//!
//! * `k=1` — pair distances, weight `1/r`;
//! * `k=2` — angles at every bonded triple, weight `1/(r_ij r_jk)`.
//!
//! The feature vector concatenates the element-pair distance blocks and the
//! element-triple angle blocks in sorted element order. The descriptor is
//! deterministic and invariant to rotation, translation and atom order.

#![allow(clippy::indexing_slicing)] // fixed grid and element blocks

use std::collections::BTreeSet;

use serde::{Deserialize, Serialize};
use xtalforge_core::{Error, Structure, neighbors};

/// Options for [`mbtr`].
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct MbtrOptions {
    /// Element symbols; empty means every element in the structure.
    pub elements: Vec<String>,
    /// Smallest distance in Angstrom.
    pub r_min: f64,
    /// Largest distance in Angstrom.
    pub r_max: f64,
    /// Number of distance and angle grid points.
    pub n_grid: usize,
    /// Gaussian standard deviation.
    pub sigma: f64,
    /// Neighbour cutoff for the k=2 triples.
    pub cutoff: f64,
    /// Include the angular block.
    pub angular: bool,
}

impl Default for MbtrOptions {
    fn default() -> Self {
        Self {
            elements: Vec::new(),
            r_min: 0.0,
            r_max: 8.0,
            n_grid: 64,
            sigma: 0.1,
            cutoff: 5.0,
            angular: true,
        }
    }
}

/// An MBTR feature vector.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct MbtrDescriptor {
    /// Feature values.
    pub features: Vec<f64>,
    /// Element order.
    pub elements: Vec<String>,
    /// Grid values for the distance blocks.
    pub distance_grid: Vec<f64>,
    /// Grid values for the angular blocks.
    pub angle_grid: Vec<f64>,
}

impl MbtrDescriptor {
    /// Serializes the descriptor to JSON.
    ///
    /// # Errors
    ///
    /// Returns [`Error::InvalidInput`] when serialization fails.
    pub fn to_json(&self) -> Result<String, Error> {
        serde_json::to_string(self)
            .map_err(|error| Error::InvalidInput(format!("could not serialize MBTR: {error}")))
    }
}

/// Computes an MBTR descriptor.
///
/// # Errors
///
/// Returns [`Error::InvalidInput`] for invalid options and the neighbour
/// errors.
pub fn mbtr(structure: &Structure, options: &MbtrOptions) -> Result<MbtrDescriptor, Error> {
    validate(options)?;
    let p1 = structure.p1()?;
    if p1.atoms.is_empty() {
        return Err(Error::InvalidInput(
            "MBTR needs at least one atom".to_owned(),
        ));
    }
    let mut elements: Vec<String> = options.elements.clone();
    if elements.is_empty() {
        elements = p1
            .atoms
            .iter()
            .map(|atom| atom.element.clone())
            .collect::<BTreeSet<_>>()
            .into_iter()
            .collect();
    }
    let mut element_index = std::collections::BTreeMap::new();
    for (index, element) in elements.iter().enumerate() {
        element_index.insert(element.clone(), index);
    }
    let n_elements = elements.len();
    let distance_grid: Vec<f64> = (0..options.n_grid)
        .map(|index| {
            options.r_min
                + (options.r_max - options.r_min) * index as f64
                    / (options.n_grid - 1).max(1) as f64
        })
        .collect();
    let angle_grid: Vec<f64> = (0..options.n_grid)
        .map(|index| std::f64::consts::PI * index as f64 / (options.n_grid - 1).max(1) as f64)
        .collect();
    let pair_count = n_elements * (n_elements + 1) / 2;
    let triple_count = n_elements * n_elements * n_elements;
    let mut distances = vec![0.0; pair_count * options.n_grid];
    let mut angles = vec![0.0; triple_count * options.n_grid];
    let raw = neighbors(structure, options.cutoff)?;
    let mut neighbour_lists: Vec<Vec<(usize, [f64; 3], f64)>> = vec![Vec::new(); p1.atoms.len()];
    for neighbor in &raw {
        if neighbor.i == neighbor.j {
            continue;
        }
        let from = p1.lattice.to_cart(p1.atoms[neighbor.i].frac);
        let to = p1.lattice.to_cart([
            p1.atoms[neighbor.j].frac[0] + neighbor.image[0] as f64,
            p1.atoms[neighbor.j].frac[1] + neighbor.image[1] as f64,
            p1.atoms[neighbor.j].frac[2] + neighbor.image[2] as f64,
        ]);
        let delta = [to[0] - from[0], to[1] - from[1], to[2] - from[2]];
        neighbour_lists[neighbor.i].push((neighbor.j, delta, neighbor.distance));
        neighbour_lists[neighbor.j].push((
            neighbor.i,
            [-delta[0], -delta[1], -delta[2]],
            neighbor.distance,
        ));
        let first = element_index
            .get(&p1.atoms[neighbor.i].element)
            .copied()
            .unwrap_or(0);
        let second = element_index
            .get(&p1.atoms[neighbor.j].element)
            .copied()
            .unwrap_or(0);
        let (low, high) = if first <= second {
            (first, second)
        } else {
            (second, first)
        };
        let block = low * n_elements - low * (low + 1) / 2 + high;
        let weight = 1.0 / neighbor.distance.max(1e-6);
        for (grid_index, value) in distance_grid.iter().enumerate() {
            let z = (neighbor.distance - value) / options.sigma;
            distances[block * options.n_grid + grid_index] += weight * (-0.5 * z * z).exp();
        }
    }
    if options.angular {
        for (centre, list) in neighbour_lists.iter().enumerate() {
            let centre_element = element_index
                .get(&p1.atoms[centre].element)
                .copied()
                .unwrap_or(0);
            for (first_index, first) in list.iter().enumerate() {
                for second in list.iter().skip(first_index + 1) {
                    let first_element = element_index
                        .get(&p1.atoms[first.0].element)
                        .copied()
                        .unwrap_or(0);
                    let second_element = element_index
                        .get(&p1.atoms[second.0].element)
                        .copied()
                        .unwrap_or(0);
                    let dot = first.1[0] * second.1[0]
                        + first.1[1] * second.1[1]
                        + first.1[2] * second.1[2];
                    let cosine = (dot / (first.2 * second.2)).clamp(-1.0, 1.0);
                    let angle = cosine.acos();
                    let block =
                        (centre_element * n_elements + first_element) * n_elements + second_element;
                    let weight = 1.0 / (first.2 * second.2).max(1e-12);
                    for (grid_index, value) in angle_grid.iter().enumerate() {
                        let z = (angle - value) / options.sigma;
                        angles[block * options.n_grid + grid_index] +=
                            weight * (-0.5 * z * z).exp();
                    }
                }
            }
        }
    }
    let mut features = distances;
    if options.angular {
        features.extend(angles);
    }
    Ok(MbtrDescriptor {
        features,
        elements,
        distance_grid,
        angle_grid,
    })
}

fn validate(options: &MbtrOptions) -> Result<(), Error> {
    if !options.r_min.is_finite() || !options.r_max.is_finite() || options.r_max <= options.r_min {
        return Err(Error::InvalidInput(
            "MBTR needs a finite r_max greater than r_min".to_owned(),
        ));
    }
    if options.n_grid < 2 || options.n_grid > 4096 {
        return Err(Error::InvalidInput(
            "MBTR n_grid must be in 2..=4096".to_owned(),
        ));
    }
    if !options.sigma.is_finite() || options.sigma <= 0.0 {
        return Err(Error::InvalidInput(
            "MBTR sigma must be positive and finite".to_owned(),
        ));
    }
    if !options.cutoff.is_finite() || options.cutoff <= 0.0 {
        return Err(Error::InvalidInput(
            "MBTR cutoff must be positive and finite".to_owned(),
        ));
    }
    Ok(())
}
