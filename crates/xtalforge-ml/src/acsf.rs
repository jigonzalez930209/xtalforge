#![allow(clippy::indexing_slicing)] // fixed-size element and angle blocks

//! Atom-centered symmetry functions (Task 21.1.3).
//!
//! The G1, G2, G3 and G4 functions follow DScribe's `ACSF` implementation:
//!
//! ```text
//! fc(r) = 0.5 (cos(pi r / r_cut) + 1)
//! G1_i  = sum_j fc(r_ij)
//! G2_i  = sum_j exp(-eta (r_ij - r_s)^2) fc(r_ij)
//! G3_i  = sum_j cos(kappa r_ij) fc(r_ij)
//! G4_i  = 2 sum_{j<k} (0.5 (1 + lambda cos theta_jik))^zeta
//!         exp(-eta (r_ij^2 + r_ik^2 + r_jk^2)) fc(r_ij) fc(r_ik) fc(r_jk)
//! ```
//!
//! The neighbour loop excludes the centre, as DScribe does. The feature
//! layout is one `(G1, G2*, G3*)` block per element followed by one `G4*`
//! block per unordered element pair, matching DScribe's `acsf.cpp`. The
//! recorded DScribe fixture is matched within `5e-6`.

use std::collections::{BTreeMap, BTreeSet};

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
    /// Include the G1 term.
    pub include_g1: bool,
    /// G2 `(eta, r_s)` pairs.
    pub g2: Vec<[f64; 2]>,
    /// G3 `kappa` values.
    pub g3: Vec<f64>,
    /// G4 `(eta, zeta, lambda)` triples.
    pub g4: Vec<[f64; 3]>,
    /// Average over centres.
    pub average: bool,
}

impl Default for AcsfOptions {
    fn default() -> Self {
        Self {
            cutoff: 5.0,
            elements: Vec::new(),
            include_g1: true,
            g2: vec![[0.05, 0.0], [0.3, 1.5]],
            g3: Vec::new(),
            g4: vec![[0.005, 1.0, -1.0], [0.005, 4.0, 1.0]],
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
    /// Feature names in output order.
    pub names: Vec<String>,
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

#[derive(Clone)]
struct Neighbor {
    index: usize,
    distance: f64,
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
            .collect::<BTreeSet<_>>()
            .into_iter()
            .collect();
    }
    let mut element_index: BTreeMap<String, usize> = BTreeMap::new();
    for (index, element) in elements.iter().enumerate() {
        element_index.insert(element.clone(), index);
    }
    let positions: Vec<[f64; 3]> = p1
        .atoms
        .iter()
        .map(|atom| p1.lattice.to_cart(atom.frac))
        .collect();
    let raw = neighbors(structure, options.cutoff)?;
    let mut lists: Vec<Vec<Neighbor>> = vec![Vec::new(); p1.atoms.len()];
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
        let _ = (to, from);
        lists[neighbor.i].push(Neighbor {
            index: neighbor.j,
            distance: neighbor.distance,
        });
        lists[neighbor.j].push(Neighbor {
            index: neighbor.i,
            distance: neighbor.distance,
        });
    }
    let n_elements = elements.len();
    let n_g2 = options.g2.len();
    let n_g3 = options.g3.len();
    let n_g4 = options.g4.len();
    let pair_count = n_elements * (n_elements + 1) / 2;
    let n_features =
        (usize::from(options.include_g1) + n_g2 + n_g3) * n_elements + n_g4 * pair_count;
    let mut names = Vec::with_capacity(n_features);
    for element in &elements {
        if options.include_g1 {
            names.push(format!("G1_{element}"));
        }
        for (index, _) in options.g2.iter().enumerate() {
            names.push(format!("G2_{element}_{index}"));
        }
        for (index, _) in options.g3.iter().enumerate() {
            names.push(format!("G3_{element}_{index}"));
        }
    }
    for first in 0..n_elements {
        for second in first..n_elements {
            for (index, _) in options.g4.iter().enumerate() {
                names.push(format!(
                    "G4_{}_{}_{}",
                    elements[first], elements[second], index
                ));
            }
        }
    }
    let mut rows = Vec::with_capacity(p1.atoms.len());
    for (centre, atom) in p1.atoms.iter().enumerate() {
        let mut row = vec![0.0; n_features];
        let list = &lists[centre];
        let mut offset = 0usize;
        let centre_element = element_index.get(&atom.element).copied().unwrap_or(0);
        let _ = centre_element;
        for (target_index, target) in elements.iter().enumerate() {
            let mut g1 = 0.0;
            let mut g2 = vec![0.0; n_g2];
            let mut g3 = vec![0.0; n_g3];
            for neighbor in list {
                if p1.atoms[neighbor.index].element != *target {
                    continue;
                }
                let fc = cutoff_function(neighbor.distance, options.cutoff);
                if fc == 0.0 {
                    continue;
                }
                g1 += fc;
                for (index, [eta, shift]) in options.g2.iter().enumerate() {
                    let delta = neighbor.distance - shift;
                    g2[index] += (-eta * delta * delta).exp() * fc;
                }
                for (index, kappa) in options.g3.iter().enumerate() {
                    g3[index] += (kappa * neighbor.distance).cos() * fc;
                }
            }
            if options.include_g1 {
                row[offset] = g1;
                offset += 1;
            }
            for value in g2 {
                row[offset] = value;
                offset += 1;
            }
            for value in g3 {
                row[offset] = value;
                offset += 1;
            }
            let _ = target_index;
        }
        for first in 0..n_elements {
            for second in first..n_elements {
                let mut g4 = vec![0.0; n_g4];
                for neighbor_j in list {
                    if element_index
                        .get(&p1.atoms[neighbor_j.index].element)
                        .copied()
                        != Some(first)
                    {
                        continue;
                    }
                    for neighbor_k in list {
                        if neighbor_k.index >= neighbor_j.index {
                            continue;
                        }
                        if element_index
                            .get(&p1.atoms[neighbor_k.index].element)
                            .copied()
                            != Some(second)
                        {
                            continue;
                        }
                        let jk = [
                            positions[neighbor_k.index][0] - positions[neighbor_j.index][0],
                            positions[neighbor_k.index][1] - positions[neighbor_j.index][1],
                            positions[neighbor_k.index][2] - positions[neighbor_j.index][2],
                        ];
                        let r_jk = (jk[0] * jk[0] + jk[1] * jk[1] + jk[2] * jk[2]).sqrt();
                        if r_jk > options.cutoff {
                            continue;
                        }
                        let costheta = 0.5 / (neighbor_j.distance * neighbor_k.distance)
                            * (neighbor_j.distance * neighbor_j.distance
                                + neighbor_k.distance * neighbor_k.distance
                                - r_jk * r_jk);
                        let fc4 = cutoff_function(neighbor_j.distance, options.cutoff)
                            * cutoff_function(neighbor_k.distance, options.cutoff)
                            * cutoff_function(r_jk, options.cutoff);
                        if fc4 == 0.0 {
                            continue;
                        }
                        for (index, [eta, zeta, lambda]) in options.g4.iter().enumerate() {
                            g4[index] += 2.0
                                * (0.5 * (1.0 + lambda * costheta)).powf(*zeta)
                                * (-eta
                                    * (neighbor_j.distance * neighbor_j.distance
                                        + neighbor_k.distance * neighbor_k.distance
                                        + r_jk * r_jk))
                                    .exp()
                                * fc4;
                        }
                    }
                }
                let block = first * n_elements - first * (first + 1) / 2 + second;
                for (index, value) in g4.iter().enumerate() {
                    row[offset + block * n_g4 + index] = *value;
                }
            }
        }
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
        names,
    })
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
    if options
        .g2
        .iter()
        .any(|[eta, shift]| !eta.is_finite() || *eta < 0.0 || !shift.is_finite())
    {
        return Err(Error::InvalidInput(
            "ACSF G2 parameters must be finite and non-negative eta".to_owned(),
        ));
    }
    if options.g3.iter().any(|value| !value.is_finite()) {
        return Err(Error::InvalidInput(
            "ACSF G3 parameters must be finite".to_owned(),
        ));
    }
    if options.g4.iter().any(|[eta, zeta, lambda]| {
        !eta.is_finite() || !zeta.is_finite() || *zeta <= 0.0 || !lambda.is_finite()
    }) {
        return Err(Error::InvalidInput(
            "ACSF G4 parameters must be finite with positive zeta".to_owned(),
        ));
    }
    Ok(())
}
