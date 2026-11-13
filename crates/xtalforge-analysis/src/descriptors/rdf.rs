//! Radial-distribution descriptor (Task 16.3.3).

use std::collections::BTreeMap;

use serde::{Deserialize, Serialize};
use xtalforge_core::{Error, RadialDistribution, RdfOptions, Structure, radial_distribution};

use super::to_json;

/// Total and partial radial distribution function of a structure.
///
/// The values come from [`radial_distribution`](xtalforge_core::radial_distribution)
/// (Task 9.3.4). Bin centres are in Angstrom; `total` is the total `g(r)` and
/// `partial` holds one curve per unordered element pair, keyed `A-B`.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct RdfDescriptor {
    /// Bin centres in Angstrom.
    pub r: Vec<f64>,
    /// Total `g(r)`.
    pub total: Vec<f64>,
    /// Partial `g_ab(r)` per unordered element pair.
    pub partial: BTreeMap<String, Vec<f64>>,
    /// Running coordination number from the total `g(r)`.
    pub coordination: Vec<f64>,
    /// Number density in atoms per cubic Angstrom.
    pub density: f64,
}

impl RdfDescriptor {
    /// Names of the partial curves in sorted order.
    #[must_use]
    pub fn partial_names(&self) -> Vec<String> {
        self.partial.keys().cloned().collect()
    }

    /// Feature vector with the total `g(r)` only.
    #[must_use]
    pub fn to_vector(&self) -> Vec<f64> {
        self.total.clone()
    }

    /// Feature vector with the total `g(r)` followed by every partial curve in
    /// sorted key order.
    #[must_use]
    pub fn to_vector_with_partials(&self) -> Vec<f64> {
        let mut vector = self.total.clone();
        for values in self.partial.values() {
            vector.extend_from_slice(values);
        }
        vector
    }

    /// Serializes the descriptor as pretty JSON.
    ///
    /// # Errors
    ///
    /// Returns [`Error::InvalidInput`] when the serializer fails.
    pub fn to_json(&self) -> Result<String, Error> {
        to_json(self, "RDF descriptor")
    }
}

impl From<RadialDistribution> for RdfDescriptor {
    fn from(distribution: RadialDistribution) -> Self {
        Self {
            r: distribution.r,
            total: distribution.total,
            partial: distribution.partial,
            coordination: distribution.coordination,
            density: distribution.density,
        }
    }
}

/// Computes the radial-distribution descriptor of a structure.
///
/// # Errors
///
/// Returns the errors of [`radial_distribution`](xtalforge_core::radial_distribution):
/// a non-positive range or bin width, an empty structure, or a neighbour-search
/// failure.
pub fn rdf_descriptor(structure: &Structure, options: &RdfOptions) -> Result<RdfDescriptor, Error> {
    Ok(RdfDescriptor::from(radial_distribution(
        structure, options,
    )?))
}
