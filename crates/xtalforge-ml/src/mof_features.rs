//! MOF-specific feature vectors (Task 21.1.4).
//!
//! Two blocks are concatenated:
//!
//! * revised autocorrelations (RACs) for atomic number, Pauling
//!   electronegativity, covalent radius and mass at bond-path lengths one to
//!   three. For property `p` and length `k`,
//!   `RAC(p,k) = sum_{i<j, d_ij=k} (p_i - <p>)(p_j - <p>) / N_pairs`.
//! * the pore geometry vector: PLD, LCD, accessible surface areas, accessible
//!   volume, void fraction, channel dimension, density and probe radius.
//!
//! The bond graph is the Phase 1 perception; the pore numbers come from the
//! Phase 10 distance grid. The result is a feature table for MOF screening
//! models, not a new descriptor theory.

#![allow(clippy::indexing_slicing)] // fixed property and descriptor arrays

use std::collections::{BTreeMap, VecDeque};

use serde::{Deserialize, Serialize};
use xtalforge_core::{Error, Structure};
use xtalforge_mof::platform::porosity::{PoreDescriptors, PoreOptions, pore_descriptors};

/// Options for [`mof_features`].
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct MofFeaturesOptions {
    /// Longest bond-path length for the RACs.
    pub max_path: usize,
    /// Pore grid spacing in Angstrom.
    pub pore_spacing: f64,
    /// Probe radius in Angstrom.
    pub probe_radius: f64,
    /// Whether to include the pore block.
    pub pores: bool,
}

impl Default for MofFeaturesOptions {
    fn default() -> Self {
        Self {
            max_path: 3,
            pore_spacing: 0.8,
            probe_radius: 1.2,
            pores: true,
        }
    }
}

/// The MOF feature vector.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct MofFeatures {
    /// Revised autocorrelations, property-major then path length.
    pub racs: Vec<f64>,
    /// Pore geometry vector, empty when disabled.
    pub pore_vector: Vec<f64>,
    /// Full concatenated vector.
    pub features: Vec<f64>,
    /// Property names in RAC order.
    pub properties: Vec<String>,
    /// Pore descriptor names in vector order.
    pub pore_names: Vec<String>,
}

impl MofFeatures {
    /// Serializes the feature vector to JSON.
    ///
    /// # Errors
    ///
    /// Returns [`Error::InvalidInput`] when serialization fails.
    pub fn to_json(&self) -> Result<String, Error> {
        serde_json::to_string(self).map_err(|error| {
            Error::InvalidInput(format!("could not serialize MOF features: {error}"))
        })
    }
}

/// Computes the MOF feature vector.
///
/// # Errors
///
/// Returns [`Error::InvalidInput`] for invalid options and the bond or pore
/// errors.
pub fn mof_features(
    structure: &Structure,
    options: &MofFeaturesOptions,
) -> Result<MofFeatures, Error> {
    if options.max_path == 0 || options.max_path > 6 {
        return Err(Error::InvalidInput("max_path must be in 1..=6".to_owned()));
    }
    if !options.pore_spacing.is_finite() || options.pore_spacing <= 0.0 {
        return Err(Error::InvalidInput(
            "pore spacing must be positive and finite".to_owned(),
        ));
    }
    let (racs, properties) = revised_autocorrelations(structure, options.max_path)?;
    let (pore_vector, pore_names) = if options.pores {
        let report = pore_descriptors(
            structure,
            &PoreOptions {
                spacing: options.pore_spacing,
                probe_radius: options.probe_radius,
                channel_probe: options.probe_radius,
            },
        )?;
        pore_geometry_vector(&report.descriptors)
    } else {
        (Vec::new(), Vec::new())
    };
    let mut features = racs.clone();
    features.extend(&pore_vector);
    Ok(MofFeatures {
        racs,
        pore_vector,
        features,
        properties,
        pore_names,
    })
}

/// Revised autocorrelations for the four element properties.
///
/// # Errors
///
/// Returns the bond graph errors.
pub fn revised_autocorrelations(
    structure: &Structure,
    max_path: usize,
) -> Result<(Vec<f64>, Vec<String>), Error> {
    let p1 = structure.p1()?;
    let canonical = crate::graph::canonical_p1(structure)?;
    let graph = xtalforge_core::bonding::perceive_bonds(
        &canonical,
        xtalforge_core::bonding::BondRules::default(),
    )?;
    let mut properties: Vec<(&str, Vec<f64>)> = Vec::new();
    let mut z = Vec::with_capacity(p1.atoms.len());
    let mut en = Vec::with_capacity(p1.atoms.len());
    let mut radius = Vec::with_capacity(p1.atoms.len());
    let mut mass = Vec::with_capacity(p1.atoms.len());
    for atom in &p1.atoms {
        let element = xtalforge_core::elements::element_by_symbol(&atom.element)?;
        z.push(f64::from(element.z));
        en.push(element.electronegativity.unwrap_or(0.0));
        radius.push(element.covalent_radius);
        mass.push(element.mass);
    }
    properties.push(("Z", z));
    properties.push(("EN", en));
    properties.push(("cov_radius", radius));
    properties.push(("mass", mass));
    let n = p1.atoms.len();
    let mut features = Vec::with_capacity(properties.len() * max_path);
    let mut names = Vec::new();
    for (name, values) in &properties {
        let mean = values.iter().sum::<f64>() / n.max(1) as f64;
        let centered: Vec<f64> = values.iter().map(|value| value - mean).collect();
        for path in 1..=max_path {
            let mut sum = 0.0;
            let mut pairs = 0usize;
            for start in 0..n {
                let mut distances = vec![usize::MAX; n];
                distances[start] = 0;
                let mut queue = VecDeque::from([start]);
                while let Some(current) = queue.pop_front() {
                    let distance = distances[current];
                    if distance >= path {
                        continue;
                    }
                    for next in graph.bonds_of(current) {
                        if *next >= n {
                            continue;
                        }
                        if distances[*next] == usize::MAX {
                            distances[*next] = distance + 1;
                            queue.push_back(*next);
                        }
                    }
                }
                for target in (start + 1)..n {
                    if distances[target] == path {
                        sum += centered[start] * centered[target];
                        pairs += 1;
                    }
                }
            }
            features.push(sum / pairs.max(1) as f64);
            names.push(format!("RAC_{name}_{path}"));
        }
    }
    Ok((features, names))
}

/// The pore geometry vector and its names.
#[must_use]
pub fn pore_geometry_vector(descriptors: &PoreDescriptors) -> (Vec<f64>, Vec<String>) {
    let values = vec![
        descriptors.pld,
        descriptors.lcd,
        descriptors.asa_m2_g,
        descriptors.asa_m2_cm3,
        descriptors.accessible_volume,
        descriptors.void_fraction,
        descriptors.channel_dimension as f64,
        descriptors.density,
        descriptors.probe_radius,
    ];
    let names = vec![
        "pld",
        "lcd",
        "asa_m2_g",
        "asa_m2_cm3",
        "accessible_volume",
        "void_fraction",
        "channel_dimension",
        "density",
        "probe_radius",
    ];
    (values, names.into_iter().map(str::to_owned).collect())
}

/// Groups RACs by property for reporting.
#[must_use]
pub fn rac_groups(features: &MofFeatures) -> BTreeMap<String, Vec<f64>> {
    let mut grouped: BTreeMap<String, Vec<f64>> = BTreeMap::new();
    for (index, name) in features.properties.iter().enumerate() {
        grouped.insert(name.clone(), vec![features.racs[index]]);
    }
    grouped
}
