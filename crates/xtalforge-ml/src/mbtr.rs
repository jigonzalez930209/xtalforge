#![allow(clippy::indexing_slicing, clippy::needless_range_loop)] // fixed grid and element blocks

//! Many-body tensor representation (Task 21.1.3).
//!
//! The implementation follows DScribe's MBTR conventions:
//!
//! * `k=2` — pair distances with weight one;
//! * `k=3` — angles at every unordered triple, in **degrees**, with weight one;
//! * Gaussian broadening is integrated over each bin (the difference of the
//!   normal CDF at the bin edges divided by the bin width), which preserves the
//!   norm on a coarse grid.
//!
//! The feature vector concatenates the enabled `k` blocks in element-pair and
//! element-triple order. A single `MbtrGrid` sets `min`, `max`, `n` and
//! `sigma` for one block; DScribe allows a different grid per k and so does
//! this API. The recorded DScribe fixtures for both blocks are matched within
//! `5e-6`.

use std::collections::{BTreeMap, BTreeSet};

use serde::{Deserialize, Serialize};
use xtalforge_core::{Error, Structure};

/// One MBTR grid.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct MbtrGrid {
    /// Smallest grid value.
    pub min: f64,
    /// Largest grid value.
    pub max: f64,
    /// Number of bins.
    pub n: usize,
    /// Gaussian standard deviation.
    pub sigma: f64,
}

impl MbtrGrid {
    /// A grid with the documented defaults.
    #[must_use]
    pub fn new(min: f64, max: f64, n: usize, sigma: f64) -> Self {
        Self { min, max, n, sigma }
    }

    fn validate(&self, name: &str) -> Result<(), Error> {
        if !self.min.is_finite() || !self.max.is_finite() || self.max <= self.min {
            return Err(Error::InvalidInput(format!(
                "MBTR {name} grid needs a finite max greater than min"
            )));
        }
        if self.n < 2 || self.n > 4096 {
            return Err(Error::InvalidInput(format!(
                "MBTR {name} grid n must be in 2..=4096"
            )));
        }
        if !self.sigma.is_finite() || self.sigma <= 0.0 {
            return Err(Error::InvalidInput(format!(
                "MBTR {name} sigma must be positive and finite"
            )));
        }
        Ok(())
    }

    fn values(&self) -> Vec<f64> {
        let dx = (self.max - self.min) / (self.n - 1) as f64;
        (0..self.n)
            .map(|index| self.min + index as f64 * dx)
            .collect()
    }
}

/// Options for [`mbtr`].
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct MbtrOptions {
    /// Element symbols; empty means every element in the structure.
    pub elements: Vec<String>,
    /// Distance block (`k=2`); `None` disables it.
    pub k2: Option<MbtrGrid>,
    /// Angle block (`k=3`) in degrees; `None` disables it.
    pub k3: Option<MbtrGrid>,
    /// Neighbour cutoff for the angle triples in Angstrom.
    pub cutoff: f64,
}

impl Default for MbtrOptions {
    fn default() -> Self {
        Self {
            elements: Vec::new(),
            k2: Some(MbtrGrid::new(0.0, 8.0, 64, 0.1)),
            k3: Some(MbtrGrid::new(0.0, 180.0, 64, 0.1)),
            cutoff: 5.0,
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
    /// Distance grid values, empty when `k2` is disabled.
    pub distance_grid: Vec<f64>,
    /// Angle grid values in degrees, empty when `k3` is disabled.
    pub angle_grid: Vec<f64>,
    /// Feature names in output order.
    pub names: Vec<String>,
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
/// Returns [`Error::InvalidInput`] for invalid options and an empty structure.
pub fn mbtr(structure: &Structure, options: &MbtrOptions) -> Result<MbtrDescriptor, Error> {
    if options.k2.is_none() && options.k3.is_none() {
        return Err(Error::InvalidInput(
            "MBTR needs at least one of k2 or k3".to_owned(),
        ));
    }
    if !options.cutoff.is_finite() || options.cutoff <= 0.0 {
        return Err(Error::InvalidInput(
            "MBTR cutoff must be positive and finite".to_owned(),
        ));
    }
    if let Some(grid) = &options.k2 {
        grid.validate("distance")?;
    }
    if let Some(grid) = &options.k3 {
        grid.validate("angle")?;
    }
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
    let mut element_index: BTreeMap<String, usize> = BTreeMap::new();
    for (index, element) in elements.iter().enumerate() {
        element_index.insert(element.clone(), index);
    }
    let n_elements = elements.len();
    let positions: Vec<[f64; 3]> = p1
        .atoms
        .iter()
        .map(|atom| p1.lattice.to_cart(atom.frac))
        .collect();
    let mut distances = vec![vec![0.0; p1.atoms.len()]; p1.atoms.len()];
    for i in 0..p1.atoms.len() {
        for j in 0..p1.atoms.len() {
            if i == j {
                continue;
            }
            let (_, distance) = p1.lattice.minimum_image([
                p1.atoms[j].frac[0] - p1.atoms[i].frac[0],
                p1.atoms[j].frac[1] - p1.atoms[i].frac[1],
                p1.atoms[j].frac[2] - p1.atoms[i].frac[2],
            ]);
            distances[i][j] = distance;
        }
    }
    let mut features = Vec::new();
    let mut names = Vec::new();
    let distance_grid = options
        .k2
        .as_ref()
        .map(MbtrGrid::values)
        .unwrap_or_default();
    let angle_grid = options
        .k3
        .as_ref()
        .map(MbtrGrid::values)
        .unwrap_or_default();
    if let Some(grid) = &options.k2 {
        let pair_count = n_elements * (n_elements + 1) / 2;
        let mut blocks = vec![vec![0.0; grid.n]; pair_count];
        for i in 0..p1.atoms.len() {
            for j in (i + 1)..p1.atoms.len() {
                let first = element_index
                    .get(&p1.atoms[i].element)
                    .copied()
                    .unwrap_or(0);
                let second = element_index
                    .get(&p1.atoms[j].element)
                    .copied()
                    .unwrap_or(0);
                let (low, high) = if first <= second {
                    (first, second)
                } else {
                    (second, first)
                };
                let block = low * n_elements - low * (low + 1) / 2 + high;
                let pdf = integrated_gaussian(grid, distances[i][j]);
                for (index, value) in pdf.iter().enumerate() {
                    blocks[block][index] += value;
                }
            }
        }
        for (block, values) in blocks.iter().enumerate() {
            features.extend(values);
            for (index, _) in values.iter().enumerate() {
                names.push(format!("MBTR_k2_{block}_{index}"));
            }
        }
    }
    if let Some(grid) = &options.k3 {
        let triple_count = n_elements * n_elements * (n_elements + 1) / 2;
        let mut blocks = vec![vec![0.0; grid.n]; triple_count];
        for centre in 0..p1.atoms.len() {
            for first in 0..p1.atoms.len() {
                if first == centre {
                    continue;
                }
                if distances[centre][first] > options.cutoff {
                    continue;
                }
                for second in (first + 1)..p1.atoms.len() {
                    if second == centre {
                        continue;
                    }
                    if distances[centre][second] > options.cutoff
                        || distances[first][second] > options.cutoff
                    {
                        continue;
                    }
                    let first_element = element_index
                        .get(&p1.atoms[first].element)
                        .copied()
                        .unwrap_or(0);
                    let second_element = element_index
                        .get(&p1.atoms[second].element)
                        .copied()
                        .unwrap_or(0);
                    let centre_element = element_index
                        .get(&p1.atoms[centre].element)
                        .copied()
                        .unwrap_or(0);
                    let key = element_triple_key(
                        centre_element,
                        first_element,
                        second_element,
                        n_elements,
                    );
                    let a = [
                        positions[first][0] - positions[centre][0],
                        positions[first][1] - positions[centre][1],
                        positions[first][2] - positions[centre][2],
                    ];
                    let b = [
                        positions[second][0] - positions[centre][0],
                        positions[second][1] - positions[centre][1],
                        positions[second][2] - positions[centre][2],
                    ];
                    let dot = a[0] * b[0] + a[1] * b[1] + a[2] * b[2];
                    let cosine = (dot / (distances[centre][first] * distances[centre][second]))
                        .clamp(-1.0, 1.0);
                    let angle = cosine.acos() * 180.0 / std::f64::consts::PI;
                    let pdf = integrated_gaussian(grid, angle);
                    for (index, value) in pdf.iter().enumerate() {
                        blocks[key][index] += value;
                    }
                }
            }
        }
        for (block, values) in blocks.iter().enumerate() {
            features.extend(values);
            for (index, _) in values.iter().enumerate() {
                names.push(format!("MBTR_k3_{block}_{index}"));
            }
        }
    }
    Ok(MbtrDescriptor {
        features,
        elements,
        distance_grid,
        angle_grid,
        names,
    })
}

fn element_triple_key(centre: usize, first: usize, second: usize, n_elements: usize) -> usize {
    // DScribe enumerates the three-dimensional element array with k >= i and
    // loops j, i, k: m = j E (E+1)/2 + k + i E - i(i+1)/2.
    let (low, high) = if first <= second {
        (first, second)
    } else {
        (second, first)
    };
    centre * n_elements * (n_elements + 1) / 2 + high + low * n_elements - low * (low + 1) / 2
}

fn integrated_gaussian(grid: &MbtrGrid, center: f64) -> Vec<f64> {
    let dx = (grid.max - grid.min) / (grid.n - 1) as f64;
    let start = grid.min - dx / 2.0;
    let scale = grid.sigma * std::f64::consts::SQRT_2;
    let mut cdf = Vec::with_capacity(grid.n + 1);
    for index in 0..=grid.n {
        let x = start + index as f64 * dx;
        cdf.push(0.5 * (1.0 + libm::erf((x - center) / scale)));
    }
    (0..grid.n)
        .map(|index| (cdf[index + 1] - cdf[index]) / dx)
        .collect()
}

/// Re-export the grid type in the crate root.
pub use MbtrGrid as Grid;
