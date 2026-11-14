//! Atom-centered symmetry functions (Task 21.1.3).
//!
//! The G2 radial and G4 angular functions follow Behler and Parrinello:
//!
//! ```text
//! G2_i = sum_j exp(-eta (r_ij - r_s)^2) fc(r_ij)
//! G4_i = 2^(1-zeta) sum_{j,k != i} (1 + lambda cos theta_jik)^zeta
//!        exp(-eta (r_ij^2 + r_ik^2 + r_jk^2)) fc(r_ij) fc(r_ik) fc(r_jk)
//! fc(r) = 0.5 (cos(pi r / r_cut) + 1) for r < r_cut, 0 otherwise
//! ```
//!
//! Every element pair and triple gets its own block, in the sorted element
//! order, so the feature vector has a stable layout. The descriptor is
//! rotation and translation invariant and depends only on the P1 cell.

#![allow(clippy::indexing_slicing)] // fixed-size element and angle blocks

use std::collections::BTreeMap;

use serde::{Deserialize, Serialize};
use xtalforge_core::{Error, Structure, neighbors};

/// Options for [`acsf`].
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct AcsfOptions {
    /// Cutoff in Angstrom.
    pub cutoff: f64,
    /// Element symbols; empty means every element in the structure.
    pub elements: Vec<String>,
    /// G2 exponent values.
    pub eta: Vec<f64>,
    /// G2 shift values.
    pub r_shift: Vec<f64>,
    /// G4 exponent values.
    pub eta_angular: Vec<f64>,
    /// G4 exponent `zeta`.
    pub zeta: Vec<f64>,
    /// G4 sign `lambda`.
    pub lambda: Vec<f64>,
    /// Average over centres.
    pub average: bool,
}

impl Default for AcsfOptions {
    fn default() -> Self {
        Self {
            cutoff: 5.0,
            elements: Vec::new(),
            eta: vec![0.05, 0.3],
            r_shift: vec![0.0, 1.5],
            eta_angular: vec![0.005],
            zeta: vec![1.0, 4.0],
            lambda: vec![-1.0, 1.0],
            average: false,
        }
    }
}

/// An ACSF feature matrix.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct AcsfDescriptor {
    /// One row per centre, or one averaged row.
    pub rows: Vec<Vec<f64>>,
    /// Element order used by the blocks.
    pub elements: Vec<String>,
    /// Number of features per row.
    pub n_features: usize,
}

impl AcsfDescriptor {
    /// Serializes the descriptor to JSON.
    ///
    /// # Errors
    ///
    /// Returns [`Error::InvalidInput`] when serialization fails.
    pub fn to_json(&self) -> Result<String, Error> {
        serde_json::to_string(self)
            .map_err(|error| Error::InvalidInput(format!("could not serialize ACSF: {error}")))
    }
}

/// Computes ACSF features for a structure.
///
/// # Errors
///
/// Returns [`Error::InvalidInput`] for invalid options or an empty structure.
pub fn acsf(structure: &Structure, options: &AcsfOptions) -> Result<AcsfDescriptor, Error> {
    validate(options)?;
    let p1 = structure.p1()?;
    if p1.atoms.is_empty() {
        return Err(Error::InvalidInput(
            "ACSF needs at least one atom".to_owned(),
        ));
    }
    let mut elements: Vec<String> = options.elements.clone();
    if elements.is_empty() {
        elements = p1
            .atoms
            .iter()
            .map(|atom| atom.element.clone())
            .collect::<std::collections::BTreeSet<_>>()
            .into_iter()
            .collect();
    }
    let mut element_index: BTreeMap<String, usize> = BTreeMap::new();
    for (index, element) in elements.iter().enumerate() {
        element_index.insert(element.clone(), index);
    }
    let raw = neighbors(structure, options.cutoff)?;
    let mut neighbours: Vec<Vec<(usize, [f64; 3], f64)>> = vec![Vec::new(); p1.atoms.len()];
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
        neighbours[neighbor.i].push((neighbor.j, delta, neighbor.distance));
        neighbours[neighbor.j].push((
            neighbor.i,
            [-delta[0], -delta[1], -delta[2]],
            neighbor.distance,
        ));
    }
    let n_g2 = options.eta.len() * options.r_shift.len() * elements.len();
    let n_g4 = options.eta_angular.len()
        * options.zeta.len()
        * options.lambda.len()
        * elements.len()
        * (elements.len() + 1)
        / 2;
    let n_features = n_g2 + n_g4;
    let mut rows = Vec::with_capacity(p1.atoms.len());
    for (centre, atom) in p1.atoms.iter().enumerate() {
        let centre_element = element_index.get(&atom.element).copied().unwrap_or(0);
        let mut row = vec![0.0; n_features];
        let list = &neighbours[centre];
        let mut offset = 0;
        for (target_index, target) in elements.iter().enumerate() {
            let mut radial = vec![0.0; options.eta.len() * options.r_shift.len()];
            for (j, _, distance) in list {
                if p1.atoms[*j].element != *target {
                    continue;
                }
                let fc = cutoff_function(*distance, options.cutoff);
                if fc == 0.0 {
                    continue;
                }
                for (eta_index, eta) in options.eta.iter().enumerate() {
                    for (shift_index, shift) in options.r_shift.iter().enumerate() {
                        let delta = distance - shift;
                        radial[eta_index * options.r_shift.len() + shift_index] +=
                            (-eta * delta * delta).exp() * fc;
                    }
                }
            }
            for (value_index, value) in radial.iter().enumerate() {
                row[offset + target_index * radial.len() + value_index] = *value;
            }
        }
        offset += n_g2;
        for first in 0..elements.len() {
            for second in first..elements.len() {
                let mut angular =
                    vec![
                        0.0;
                        options.eta_angular.len() * options.zeta.len() * options.lambda.len()
                    ];
                for (j, first_delta, first_distance) in list {
                    if element_index.get(&p1.atoms[*j].element).copied() != Some(first) {
                        continue;
                    }
                    for (k, second_delta, second_distance) in list {
                        if j == k {
                            continue;
                        }
                        if element_index.get(&p1.atoms[*k].element).copied() != Some(second) {
                            continue;
                        }
                        let dot = first_delta[0] * second_delta[0]
                            + first_delta[1] * second_delta[1]
                            + first_delta[2] * second_delta[2];
                        let cosine = (dot / (first_distance * second_distance)).clamp(-1.0, 1.0);
                        let third = distance_between(first_delta, second_delta);
                        let fc = cutoff_function(*first_distance, options.cutoff)
                            * cutoff_function(*second_distance, options.cutoff)
                            * cutoff_function(third, options.cutoff);
                        if fc == 0.0 {
                            continue;
                        }
                        for (eta_index, eta) in options.eta_angular.iter().enumerate() {
                            let decay = (-eta
                                * (first_distance * first_distance
                                    + second_distance * second_distance
                                    + third * third))
                                .exp();
                            for (zeta_index, zeta) in options.zeta.iter().enumerate() {
                                for (lambda_index, lambda) in options.lambda.iter().enumerate() {
                                    let index =
                                        eta_index * options.zeta.len() * options.lambda.len()
                                            + zeta_index * options.lambda.len()
                                            + lambda_index;
                                    angular[index] += 2.0_f64.powf(1.0 - zeta)
                                        * (1.0 + lambda * cosine).powf(*zeta)
                                        * decay
                                        * fc;
                                }
                            }
                        }
                    }
                }
                let width = angular.len();
                let block = first * elements.len() - first * (first + 1) / 2 + second;
                for (value_index, value) in angular.iter().enumerate() {
                    row[offset + block * width + value_index] = *value;
                }
            }
        }
        let _ = centre_element;
        rows.push(row);
    }
    if options.average {
        let mut mean = vec![0.0; n_features];
        let count = rows.len().max(1) as f64;
        for row in &rows {
            for (index, value) in row.iter().enumerate() {
                mean[index] += value / count;
            }
        }
        rows = vec![mean];
    }
    Ok(AcsfDescriptor {
        rows,
        elements,
        n_features,
    })
}

fn distance_between(first: &[f64; 3], second: &[f64; 3]) -> f64 {
    let delta = [
        first[0] - second[0],
        first[1] - second[1],
        first[2] - second[2],
    ];
    (delta[0] * delta[0] + delta[1] * delta[1] + delta[2] * delta[2]).sqrt()
}

fn cutoff_function(distance: f64, cutoff: f64) -> f64 {
    if distance <= 0.0 || distance >= cutoff {
        0.0
    } else {
        0.5 * ((std::f64::consts::PI * distance / cutoff).cos() + 1.0)
    }
}

fn validate(options: &AcsfOptions) -> Result<(), Error> {
    if !options.cutoff.is_finite() || options.cutoff <= 0.0 {
        return Err(Error::InvalidInput(
            "the ACSF cutoff must be positive and finite".to_owned(),
        ));
    }
    if options.eta.is_empty()
        || options.r_shift.is_empty()
        || options.eta_angular.is_empty()
        || options.zeta.is_empty()
        || options.lambda.is_empty()
    {
        return Err(Error::InvalidInput(
            "every ACSF parameter list must be non-empty".to_owned(),
        ));
    }
    if options
        .eta
        .iter()
        .any(|value| !value.is_finite() || *value < 0.0)
        || options
            .eta_angular
            .iter()
            .any(|value| !value.is_finite() || *value < 0.0)
    {
        return Err(Error::InvalidInput(
            "ACSF exponents must be finite and non-negative".to_owned(),
        ));
    }
    if options
        .zeta
        .iter()
        .any(|value| !value.is_finite() || *value <= 0.0)
    {
        return Err(Error::InvalidInput(
            "ACSF zeta must be positive and finite".to_owned(),
        ));
    }
    Ok(())
}
