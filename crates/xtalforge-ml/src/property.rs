#![allow(clippy::indexing_slicing)] // bounded indices from validated loops and fixed-size vectors

//! Generic property models and dataset prediction (Tasks 21.3.2 and 21.3.3).

use serde::{Deserialize, Serialize};
use xtalforge_analysis::dataset::{Dataset, PropertyValue};
use xtalforge_core::{Error, Structure};

/// Anything that maps a descriptor vector to one scalar.
pub trait Predictor: Send + Sync {
    /// Model name.
    fn name(&self) -> &str;
    /// Expected input width.
    fn n_features(&self) -> usize;
    /// Predicts one scalar.
    ///
    /// # Errors
    ///
    /// Returns [`Error::InvalidInput`] for a wrong input width or an
    /// inference failure.
    fn predict(&self, features: &[f64]) -> Result<f64, Error>;
}

/// A linear model with optional input centring and scaling.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct PropertyModel {
    /// Model name reported with predictions.
    pub name: String,
    /// Feature weights.
    pub weights: Vec<f64>,
    /// Bias.
    pub bias: f64,
    /// Optional per-feature means subtracted before the dot product.
    pub mean: Vec<f64>,
    /// Optional per-feature scales divided before the dot product.
    pub scale: Vec<f64>,
}

impl PropertyModel {
    /// Creates a model from weights and a bias.
    ///
    /// # Errors
    ///
    /// Returns [`Error::InvalidInput`] for an empty or non-finite weight list.
    pub fn new(name: impl Into<String>, weights: Vec<f64>, bias: f64) -> Result<Self, Error> {
        if weights.is_empty() || !weights.iter().all(|value| value.is_finite()) {
            return Err(Error::InvalidInput(
                "property weights must be finite and non-empty".to_owned(),
            ));
        }
        if !bias.is_finite() {
            return Err(Error::InvalidInput(
                "the property bias must be finite".to_owned(),
            ));
        }
        Ok(Self {
            name: name.into(),
            weights,
            bias,
            mean: Vec::new(),
            scale: Vec::new(),
        })
    }

    /// Adds input centring and scaling.
    ///
    /// # Errors
    ///
    /// Returns [`Error::InvalidInput`] when the lengths do not match the
    /// weights or a scale is zero.
    pub fn with_normalization(mut self, mean: Vec<f64>, scale: Vec<f64>) -> Result<Self, Error> {
        if mean.len() != self.weights.len() || scale.len() != self.weights.len() {
            return Err(Error::InvalidInput(
                "normalization vectors must match the weight length".to_owned(),
            ));
        }
        if scale
            .iter()
            .any(|value| !value.is_finite() || *value == 0.0)
        {
            return Err(Error::InvalidInput(
                "normalization scales must be finite and non-zero".to_owned(),
            ));
        }
        self.mean = mean;
        self.scale = scale;
        Ok(self)
    }

    /// Serializes the model to JSON.
    ///
    /// # Errors
    ///
    /// Returns [`Error::InvalidInput`] when serialization fails.
    pub fn to_json(&self) -> Result<String, Error> {
        serde_json::to_string(self)
            .map_err(|error| Error::InvalidInput(format!("could not serialize model: {error}")))
    }

    /// Parses a model from JSON.
    ///
    /// # Errors
    ///
    /// Returns [`Error::InvalidInput`] for invalid JSON or weights.
    pub fn from_json(text: &str) -> Result<Self, Error> {
        let model: Self = serde_json::from_str(text)
            .map_err(|error| Error::InvalidInput(format!("could not parse model: {error}")))?;
        if model.weights.is_empty() || !model.weights.iter().all(|value| value.is_finite()) {
            return Err(Error::InvalidInput(
                "property weights must be finite and non-empty".to_owned(),
            ));
        }
        Ok(model)
    }
}

impl Predictor for PropertyModel {
    fn name(&self) -> &str {
        &self.name
    }

    fn n_features(&self) -> usize {
        self.weights.len()
    }

    fn predict(&self, features: &[f64]) -> Result<f64, Error> {
        if features.len() != self.weights.len() {
            return Err(Error::InvalidInput(format!(
                "model `{}` expects {} features, got {}",
                self.name,
                self.weights.len(),
                features.len()
            )));
        }
        let mut value = self.bias;
        for (index, feature) in features.iter().enumerate() {
            let centered = if self.mean.is_empty() {
                *feature
            } else {
                feature - self.mean[index]
            };
            let scaled = if self.scale.is_empty() {
                centered
            } else {
                centered / self.scale[index]
            };
            value += self.weights[index] * scaled;
        }
        Ok(value)
    }
}

/// One dataset prediction table.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct WorkflowPrediction {
    /// Model name.
    pub model: String,
    /// Record identifiers.
    pub ids: Vec<String>,
    /// Predicted values.
    pub predictions: Vec<f64>,
}

impl WorkflowPrediction {
    /// Serializes the prediction to JSON.
    ///
    /// # Errors
    ///
    /// Returns [`Error::InvalidInput`] when serialization fails.
    pub fn to_json(&self) -> Result<String, Error> {
        serde_json::to_string(self).map_err(|error| {
            Error::InvalidInput(format!("could not serialize prediction: {error}"))
        })
    }
}

/// Predicts one scalar for every dataset record.
///
/// # Errors
///
/// Returns the featurizer and predictor errors.
pub fn predict_dataset<F>(
    dataset: &Dataset,
    model: &dyn Predictor,
    mut featurizer: F,
) -> Result<WorkflowPrediction, Error>
where
    F: FnMut(&Structure) -> Result<Vec<f64>, Error>,
{
    let mut ids = Vec::with_capacity(dataset.len());
    let mut predictions = Vec::with_capacity(dataset.len());
    for record in dataset.iter() {
        let features = featurizer(&record.structure)?;
        ids.push(record.id.clone());
        predictions.push(model.predict(&features)?);
    }
    Ok(WorkflowPrediction {
        model: model.name().to_owned(),
        ids,
        predictions,
    })
}

/// Appends predictions as a dataset property.
///
/// # Errors
///
/// Returns [`Error::InvalidInput`] for a length mismatch or an invalid value.
pub fn append_predictions(
    dataset: &Dataset,
    prediction: &WorkflowPrediction,
    property: &str,
) -> Result<Dataset, Error> {
    if prediction.predictions.len() != dataset.len() {
        return Err(Error::InvalidInput(
            "prediction and dataset lengths differ".to_owned(),
        ));
    }
    let mut output = Dataset::new(dataset.metadata().name.clone());
    *output.metadata_mut() = dataset.metadata().clone();
    for (record, value) in dataset.iter().zip(&prediction.predictions) {
        let mut record = record.clone();
        record = record.with_property(property.to_owned(), PropertyValue::Number(*value))?;
        output.push(record)?;
    }
    Ok(output)
}
