#![allow(clippy::indexing_slicing)] // bounded indices from validated loops and fixed-size vectors

//! Ensemble predictions and uncertainty (Task 21.3.3).

use serde::{Deserialize, Serialize};
use xtalforge_core::Error;

use crate::property::{Predictor, PropertyModel};

/// One prediction with its committee disagreement.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct Prediction {
    /// Mean prediction.
    pub mean: f64,
    /// Standard deviation across members.
    pub std: f64,
    /// Individual member predictions.
    pub members: Vec<f64>,
}

/// An ensemble of linear models.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct Ensemble {
    /// Member models.
    pub members: Vec<PropertyModel>,
}

/// An ensemble prediction table.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct EnsemblePrediction {
    /// Record identifiers.
    pub ids: Vec<String>,
    /// Mean, deviation and member values per record.
    pub predictions: Vec<Prediction>,
}

impl Ensemble {
    /// Builds an ensemble; every member must have the same feature width.
    ///
    /// # Errors
    ///
    /// Returns [`Error::InvalidInput`] for an empty ensemble or mixed widths.
    pub fn new(members: Vec<PropertyModel>) -> Result<Self, Error> {
        if members.is_empty() {
            return Err(Error::InvalidInput(
                "an ensemble needs at least one member".to_owned(),
            ));
        }
        let width = members[0].weights.len();
        if members.iter().any(|member| member.weights.len() != width) {
            return Err(Error::InvalidInput(
                "ensemble members must have the same feature width".to_owned(),
            ));
        }
        Ok(Self { members })
    }

    /// Predicts one vector.
    ///
    /// # Errors
    ///
    /// Returns the member errors.
    pub fn predict(&self, features: &[f64]) -> Result<Prediction, Error> {
        let mut values = Vec::with_capacity(self.members.len());
        for member in &self.members {
            values.push(member.predict(features)?);
        }
        Ok(uncertainty(&values))
    }
}

/// Mean and standard deviation of a committee.
#[must_use]
pub fn uncertainty(values: &[f64]) -> Prediction {
    if values.is_empty() {
        return Prediction {
            mean: 0.0,
            std: 0.0,
            members: Vec::new(),
        };
    }
    let mean = values.iter().sum::<f64>() / values.len() as f64;
    let variance = values
        .iter()
        .map(|value| (value - mean) * (value - mean))
        .sum::<f64>()
        / values.len() as f64;
    Prediction {
        mean,
        std: variance.sqrt(),
        members: values.to_vec(),
    }
}
