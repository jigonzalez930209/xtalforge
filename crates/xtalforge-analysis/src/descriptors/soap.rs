#![allow(clippy::indexing_slicing, clippy::needless_range_loop)] // dense descriptor kernel

//! SOAP-like local descriptor (Task 16.3.3).

use std::collections::BTreeSet;

use serde::{Deserialize, Serialize};
use xtalforge_core::{Error, Structure};

use super::{canonical_cell, canonical_snapshot, to_json};

/// Largest Legendre order accepted by [`SoapLikeOptions`].
pub const MAX_SOAP_L: usize = 8;
/// Largest radial basis size accepted by [`SoapLikeOptions`].
pub const MAX_SOAP_RADIAL: usize = 32;
/// Largest neighbour cap accepted by [`SoapLikeOptions`].
pub const MAX_SOAP_NEIGHBORS: usize = 4096;
/// Largest feature vector accepted per centre.
pub const MAX_SOAP_FEATURES: usize = 100_000;

/// Options for [`soap_like_descriptor`].
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct SoapLikeOptions {
    /// Neighbour cutoff in Angstrom.
    pub cutoff: f64,
    /// Number of Gaussian radial basis functions.
    pub n_radial: usize,
    /// Gaussian width in Angstrom.
    pub sigma: f64,
    /// Largest Legendre polynomial order.
    pub l_max: usize,
    /// Maximum neighbours kept per centre, nearest first.
    pub max_neighbors: usize,
}

impl Default for SoapLikeOptions {
    fn default() -> Self {
        Self {
            cutoff: 5.0,
            n_radial: 4,
            sigma: 0.5,
            l_max: 2,
            max_neighbors: 128,
        }
    }
}

/// A SOAP-like descriptor around every atom plus its structure mean.
///
/// The layout of every per-centre vector is
///
/// ```text
/// [radial(s=0,n=0) ... radial(s=0,n=N-1), radial(s=1,n=0) ...,
///  angular(s=0,t=0,n,n',l) ..., angular(s=0,t=1,n,n',l) ...,
///  angular(s=1,t=1,n,n',l) ...]
/// ```
///
/// with species in the alphabetical order of `species`, `n` and `n'` varying
/// fastest before `l`, and `t >= s` for the unordered species pairs. Species
/// pairs are enumerated as `(0, 0), (0, 1), ..., (1, 1), ...` in lexicographic
/// order. The formulas are documented in the
/// [`descriptors`](crate::descriptors) module.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct SoapLikeDescriptor {
    /// Element symbols that span the species blocks, alphabetically sorted.
    pub species: Vec<String>,
    /// Radial Gaussian centres in Angstrom.
    pub radial_centers: Vec<f64>,
    /// Number of radial basis functions.
    pub n_radial: usize,
    /// Largest Legendre order.
    pub l_max: usize,
    /// Element of each centre in canonical order.
    pub center_elements: Vec<String>,
    /// Fractional position of each centre in canonical order.
    pub center_positions: Vec<[f64; 3]>,
    /// One feature vector per centre.
    pub per_center: Vec<Vec<f64>>,
    /// Arithmetic mean of `per_center`, suitable as a structure-level vector.
    pub mean: Vec<f64>,
}

impl SoapLikeDescriptor {
    /// Number of features per centre.
    #[must_use]
    pub fn feature_len(&self) -> usize {
        self.mean.len()
    }

    /// Number of centres.
    #[must_use]
    pub fn n_centers(&self) -> usize {
        self.per_center.len()
    }

    /// Structure-level feature vector (the mean over centres).
    #[must_use]
    pub fn to_vector(&self) -> Vec<f64> {
        self.mean.clone()
    }

    /// Serializes the descriptor as pretty JSON.
    ///
    /// # Errors
    ///
    /// Returns [`Error::InvalidInput`] when the serializer fails.
    pub fn to_json(&self) -> Result<String, Error> {
        to_json(self, "SOAP-like descriptor")
    }
}

/// Computes a SOAP-like local descriptor.
///
/// # Errors
///
/// Returns [`Error::InvalidInput`] for invalid options, an empty structure, an
/// unknown element, a feature vector above [`MAX_SOAP_FEATURES`], or a
/// neighbour-search or P1-expansion failure.
pub fn soap_like_descriptor(
    structure: &Structure,
    options: &SoapLikeOptions,
) -> Result<SoapLikeDescriptor, Error> {
    validate(options)?;
    let (lattice, sites) = canonical_cell(structure)?;
    let snapshot = canonical_snapshot(lattice, &sites)?;
    let mut species_set = BTreeSet::new();
    for site in &sites {
        species_set.insert(site.element.clone());
    }
    let species: Vec<String> = species_set.into_iter().collect();
    let species_index: std::collections::BTreeMap<&str, usize> = species
        .iter()
        .enumerate()
        .map(|(index, symbol)| (symbol.as_str(), index))
        .collect();
    let atom_species: Vec<usize> = sites
        .iter()
        .map(|site| {
            species_index
                .get(site.element.as_str())
                .copied()
                .unwrap_or_default()
        })
        .collect();

    let pairs = xtalforge_core::neighbors(&snapshot, options.cutoff)?;
    let mut lists: Vec<Vec<LocalNeighbor>> = (0..sites.len()).map(|_| Vec::new()).collect();
    for pair in pairs {
        let (Some(first), Some(second)) = (sites.get(pair.i), sites.get(pair.j)) else {
            continue;
        };
        if !pair.distance.is_finite() || pair.distance <= 1e-8 {
            continue;
        }
        let delta = [
            second.frac[0] + f64::from(pair.image[0]) - first.frac[0],
            second.frac[1] + f64::from(pair.image[1]) - first.frac[1],
            second.frac[2] + f64::from(pair.image[2]) - first.frac[2],
        ];
        let vector = lattice.to_cart(delta);
        let first_species = atom_species.get(pair.i).copied().unwrap_or_default();
        let second_species = atom_species.get(pair.j).copied().unwrap_or_default();
        if pair.i == pair.j {
            if let Some(list) = lists.get_mut(pair.i) {
                list.push(LocalNeighbor {
                    species: first_species,
                    distance: pair.distance,
                    vector,
                    weights: Vec::new(),
                });
            }
            continue;
        }
        if let Some(list) = lists.get_mut(pair.i) {
            list.push(LocalNeighbor {
                species: second_species,
                distance: pair.distance,
                vector,
                weights: Vec::new(),
            });
        }
        if let Some(list) = lists.get_mut(pair.j) {
            list.push(LocalNeighbor {
                species: first_species,
                distance: pair.distance,
                vector: [-vector[0], -vector[1], -vector[2]],
                weights: Vec::new(),
            });
        }
    }

    let radial_centers: Vec<f64> = (0..options.n_radial)
        .map(|index| options.cutoff * (index as f64 + 1.0) / (options.n_radial as f64 + 1.0))
        .collect();
    for list in &mut lists {
        list.sort_by(|left, right| {
            left.distance
                .total_cmp(&right.distance)
                .then_with(|| left.species.cmp(&right.species))
                .then_with(|| left.vector[0].total_cmp(&right.vector[0]))
                .then_with(|| left.vector[1].total_cmp(&right.vector[1]))
                .then_with(|| left.vector[2].total_cmp(&right.vector[2]))
        });
        list.truncate(options.max_neighbors);
        for neighbor in list.iter_mut() {
            neighbor.weights = radial_centers
                .iter()
                .map(|center| gaussian(neighbor.distance, *center, options.sigma))
                .collect();
        }
    }

    let n_species = species.len();
    let n_pairs = n_species * (n_species + 1) / 2;
    let feature_len = n_species * options.n_radial
        + n_pairs * options.n_radial * options.n_radial * (options.l_max + 1);
    if feature_len > MAX_SOAP_FEATURES {
        return Err(Error::InvalidInput(format!(
            "SOAP-like descriptor would need {feature_len} features per centre, \
             above the limit of {MAX_SOAP_FEATURES}; lower n_radial or l_max"
        )));
    }

    let mut per_center = Vec::with_capacity(lists.len());
    for list in &lists {
        let mut groups: Vec<Vec<usize>> = vec![Vec::new(); n_species];
        for (index, neighbor) in list.iter().enumerate() {
            if let Some(group) = groups.get_mut(neighbor.species) {
                group.push(index);
            }
        }
        let mut features = Vec::with_capacity(feature_len);
        for group in &groups {
            let mut sum = vec![0.0; options.n_radial];
            for index in group {
                let Some(neighbor) = list.get(*index) else {
                    continue;
                };
                for (slot, weight) in sum.iter_mut().zip(neighbor.weights.iter()) {
                    *slot += weight;
                }
            }
            let norm = if group.is_empty() {
                1.0
            } else {
                group.len() as f64
            };
            features.extend(sum.into_iter().map(|value| value / norm));
        }
        for first_species in 0..n_species {
            let first_group = groups.get(first_species).map(Vec::as_slice).unwrap_or(&[]);
            for second_species in first_species..n_species {
                let second_group = groups.get(second_species).map(Vec::as_slice).unwrap_or(&[]);
                let denominator = if first_species == second_species {
                    first_group.len() * first_group.len().saturating_sub(1)
                } else {
                    first_group.len() * second_group.len()
                };
                let norm = if denominator == 0 {
                    1.0
                } else {
                    denominator as f64
                };
                for radial_first in 0..options.n_radial {
                    for radial_second in 0..options.n_radial {
                        for order in 0..=options.l_max {
                            let mut sum = 0.0;
                            for (first_position, first_index) in first_group.iter().enumerate() {
                                let Some(left) = list.get(*first_index) else {
                                    continue;
                                };
                                for (second_position, second_index) in
                                    second_group.iter().enumerate()
                                {
                                    if first_species == second_species
                                        && first_position == second_position
                                    {
                                        continue;
                                    }
                                    let Some(right) = list.get(*second_index) else {
                                        continue;
                                    };
                                    let cosine = ((left.vector[0] * right.vector[0]
                                        + left.vector[1] * right.vector[1]
                                        + left.vector[2] * right.vector[2])
                                        / (left.distance * right.distance))
                                        .clamp(-1.0, 1.0);
                                    let (Some(left_weight), Some(right_weight)) = (
                                        left.weights.get(radial_first),
                                        right.weights.get(radial_second),
                                    ) else {
                                        continue;
                                    };
                                    sum += left_weight * right_weight * legendre(order, cosine);
                                }
                            }
                            features.push(sum / norm);
                        }
                    }
                }
            }
        }
        per_center.push(features);
    }

    let mut mean = vec![0.0; feature_len];
    for values in &per_center {
        for (slot, value) in mean.iter_mut().zip(values.iter()) {
            *slot += value;
        }
    }
    let count = per_center.len() as f64;
    for value in &mut mean {
        *value /= count;
    }
    Ok(SoapLikeDescriptor {
        species,
        radial_centers,
        n_radial: options.n_radial,
        l_max: options.l_max,
        center_elements: sites.iter().map(|site| site.element.clone()).collect(),
        center_positions: sites.iter().map(|site| site.frac).collect(),
        per_center,
        mean,
    })
}

struct LocalNeighbor {
    species: usize,
    distance: f64,
    vector: [f64; 3],
    weights: Vec<f64>,
}

fn validate(options: &SoapLikeOptions) -> Result<(), Error> {
    if !options.cutoff.is_finite() || options.cutoff <= 0.0 {
        return Err(Error::InvalidInput(
            "SOAP-like cutoff must be positive and finite".to_owned(),
        ));
    }
    if options.n_radial == 0 || options.n_radial > MAX_SOAP_RADIAL {
        return Err(Error::InvalidInput(format!(
            "SOAP-like n_radial must be in 1..={MAX_SOAP_RADIAL}, got {}",
            options.n_radial
        )));
    }
    if !options.sigma.is_finite() || options.sigma <= 0.0 {
        return Err(Error::InvalidInput(
            "SOAP-like sigma must be positive and finite".to_owned(),
        ));
    }
    if options.l_max > MAX_SOAP_L {
        return Err(Error::InvalidInput(format!(
            "SOAP-like l_max must not exceed {MAX_SOAP_L}, got {}",
            options.l_max
        )));
    }
    if options.max_neighbors == 0 || options.max_neighbors > MAX_SOAP_NEIGHBORS {
        return Err(Error::InvalidInput(format!(
            "SOAP-like max_neighbors must be in 1..={MAX_SOAP_NEIGHBORS}, got {}",
            options.max_neighbors
        )));
    }
    Ok(())
}

fn gaussian(distance: f64, center: f64, sigma: f64) -> f64 {
    let scaled = (distance - center) / sigma;
    (-0.5 * scaled * scaled).exp()
}

fn legendre(order: usize, cosine: f64) -> f64 {
    if order == 0 {
        return 1.0;
    }
    if order == 1 {
        return cosine;
    }
    let mut previous = 1.0;
    let mut current = cosine;
    for value in 2..=order {
        let next = ((2 * value - 1) as f64 * cosine * current - (value - 1) as f64 * previous)
            / value as f64;
        previous = current;
        current = next;
    }
    current
}
