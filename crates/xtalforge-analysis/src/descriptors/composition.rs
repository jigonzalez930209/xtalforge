//! Composition features (Task 16.3.3).

use std::collections::BTreeMap;

use serde::{Deserialize, Serialize};
use xtalforge_core::{Error, Structure};

use super::{canonical_cell, to_json};

/// Composition of a structure as a feature vector.
///
/// Fractions are atom fractions, not mass fractions. Occupancies are treated
/// as one. Electronegativity statistics use only the elements that define a
/// Pauling value.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct CompositionFeatures {
    /// Element symbols in alphabetical order.
    pub elements: Vec<String>,
    /// Atom count per element, in the same order as `elements`.
    pub counts: Vec<usize>,
    /// Atom fractions that sum to one, in the same order as `elements`.
    pub fractions: Vec<f64>,
    /// Total number of atoms in the P1 cell.
    pub n_atoms: usize,
    /// Mean atomic number over all atoms.
    pub mean_atomic_number: f64,
    /// Mean standard atomic mass in unified atomic mass units.
    pub mean_mass: f64,
    /// Mean Pauling electronegativity over the atoms that define one.
    pub mean_electronegativity: Option<f64>,
    /// Fraction of atoms whose element is a metal.
    pub metal_fraction: f64,
}

impl CompositionFeatures {
    /// Atom fraction of one element; zero when the element is absent.
    #[must_use]
    pub fn fraction(&self, element: &str) -> f64 {
        self.elements
            .iter()
            .position(|entry| entry == element)
            .and_then(|index| self.fractions.get(index).copied())
            .unwrap_or(0.0)
    }

    /// Atom count of one element; zero when the element is absent.
    #[must_use]
    pub fn count(&self, element: &str) -> usize {
        self.elements
            .iter()
            .position(|entry| entry == element)
            .and_then(|index| self.counts.get(index).copied())
            .unwrap_or(0)
    }

    /// Feature vector of atom fractions in `elements` order.
    #[must_use]
    pub fn to_vector(&self) -> Vec<f64> {
        self.fractions.clone()
    }

    /// Serializes the features as pretty JSON.
    ///
    /// # Errors
    ///
    /// Returns [`Error::InvalidInput`] when the serializer fails.
    pub fn to_json(&self) -> Result<String, Error> {
        to_json(self, "composition features")
    }
}

/// Computes the composition features of a structure.
///
/// # Errors
///
/// Returns [`Error::InvalidInput`] for an empty structure or an unknown
/// element, or the errors of the P1 expansion.
pub fn composition_features(structure: &Structure) -> Result<CompositionFeatures, Error> {
    let (_, sites) = canonical_cell(structure)?;
    let n_atoms = sites.len();
    let mut entries: BTreeMap<String, (usize, f64, f64, usize, usize)> = BTreeMap::new();
    let mut z_sum = 0.0;
    for site in &sites {
        let entry = entries
            .entry(site.element.clone())
            .or_insert((0, 0.0, 0.0, 0, 0));
        entry.0 += 1;
        entry.1 += site.mass;
        if let Some(value) = site.electronegativity {
            entry.2 += value;
            entry.3 += 1;
        }
        if site.is_metal {
            entry.4 += 1;
        }
        z_sum += f64::from(site.z);
    }
    let total = n_atoms as f64;
    let mut elements = Vec::with_capacity(entries.len());
    let mut counts = Vec::with_capacity(entries.len());
    let mut fractions = Vec::with_capacity(entries.len());
    let mut mass_sum = 0.0;
    let mut electronegativity_sum = 0.0;
    let mut electronegativity_count = 0;
    let mut metal_count = 0;
    for (element, (count, mass, electronegativity, en_count, metal)) in entries {
        elements.push(element);
        counts.push(count);
        fractions.push(count as f64 / total);
        mass_sum += mass;
        electronegativity_sum += electronegativity;
        electronegativity_count += en_count;
        metal_count += metal;
    }
    Ok(CompositionFeatures {
        elements,
        counts,
        fractions,
        n_atoms,
        mean_atomic_number: z_sum / total,
        mean_mass: mass_sum / total,
        mean_electronegativity: if electronegativity_count == 0 {
            None
        } else {
            Some(electronegativity_sum / electronegativity_count as f64)
        },
        metal_fraction: metal_count as f64 / total,
    })
}
