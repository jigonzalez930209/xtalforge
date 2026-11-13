//! Record and property types of the JSON Lines dataset container.

use std::collections::BTreeMap;

use serde::{Deserialize, Serialize};
use serde_json::Value;
use xtalforge_core::{Error, Structure};

/// One typed property of a dataset record.
///
/// The JSON encoding is untagged, so a property is written as a bare JSON
/// scalar or array. Floating-point properties must be finite; the constructors
/// and the reader reject `NaN` and infinities.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
#[serde(untagged)]
pub enum PropertyValue {
    /// Boolean property.
    Boolean(bool),
    /// Signed integer property.
    Integer(i64),
    /// Floating-point property.
    Number(f64),
    /// Text property.
    Text(String),
    /// Vector of finite floating-point values, for example a descriptor.
    Vector(Vec<f64>),
}

impl PropertyValue {
    /// Stable name of the property type (`boolean`, `integer`, `number`,
    /// `text` or `vector`).
    #[must_use]
    pub fn kind(&self) -> &'static str {
        match self {
            Self::Boolean(_) => "boolean",
            Self::Integer(_) => "integer",
            Self::Number(_) => "number",
            Self::Text(_) => "text",
            Self::Vector(_) => "vector",
        }
    }

    /// Checks that every stored number is finite.
    ///
    /// # Errors
    ///
    /// Returns [`Error::InvalidInput`] for a non-finite number or vector entry.
    pub fn validate(&self) -> Result<(), Error> {
        match self {
            Self::Boolean(_) | Self::Integer(_) | Self::Text(_) => Ok(()),
            Self::Number(value) if value.is_finite() => Ok(()),
            Self::Number(_) => Err(Error::InvalidInput(
                "dataset properties must be finite numbers".to_owned(),
            )),
            Self::Vector(values) if values.iter().all(|value| value.is_finite()) => Ok(()),
            Self::Vector(_) => Err(Error::InvalidInput(
                "dataset property vectors must contain finite numbers".to_owned(),
            )),
        }
    }

    /// Numeric value, with integers converted to `f64`.
    #[must_use]
    pub fn as_f64(&self) -> Option<f64> {
        match self {
            Self::Integer(value) => Some(*value as f64),
            Self::Number(value) => Some(*value),
            Self::Boolean(_) | Self::Text(_) | Self::Vector(_) => None,
        }
    }

    /// Text value.
    #[must_use]
    pub fn as_str(&self) -> Option<&str> {
        match self {
            Self::Text(value) => Some(value),
            Self::Boolean(_) | Self::Integer(_) | Self::Number(_) | Self::Vector(_) => None,
        }
    }

    /// Boolean value.
    #[must_use]
    pub fn as_bool(&self) -> Option<bool> {
        match self {
            Self::Boolean(value) => Some(*value),
            Self::Integer(_) | Self::Number(_) | Self::Text(_) | Self::Vector(_) => None,
        }
    }

    /// Vector value.
    #[must_use]
    pub fn as_vector(&self) -> Option<&[f64]> {
        match self {
            Self::Vector(values) => Some(values),
            Self::Boolean(_) | Self::Integer(_) | Self::Number(_) | Self::Text(_) => None,
        }
    }
}

impl From<bool> for PropertyValue {
    fn from(value: bool) -> Self {
        Self::Boolean(value)
    }
}

impl From<i32> for PropertyValue {
    fn from(value: i32) -> Self {
        Self::Integer(i64::from(value))
    }
}

impl From<u32> for PropertyValue {
    fn from(value: u32) -> Self {
        Self::Integer(i64::from(value))
    }
}

impl From<i64> for PropertyValue {
    fn from(value: i64) -> Self {
        Self::Integer(value)
    }
}

impl From<f64> for PropertyValue {
    fn from(value: f64) -> Self {
        Self::Number(value)
    }
}

impl From<&str> for PropertyValue {
    fn from(value: &str) -> Self {
        Self::Text(value.to_owned())
    }
}

impl From<String> for PropertyValue {
    fn from(value: String) -> Self {
        Self::Text(value)
    }
}

impl From<Vec<f64>> for PropertyValue {
    fn from(values: Vec<f64>) -> Self {
        Self::Vector(values)
    }
}

/// One structure and its properties inside a [`Dataset`](crate::dataset::Dataset).
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct DatasetRecord {
    /// Unique record identifier inside the dataset.
    pub id: String,
    /// The structure carried by this record.
    pub structure: Structure,
    /// Named properties; keys are sorted on write.
    #[serde(default)]
    pub properties: BTreeMap<String, PropertyValue>,
}

impl DatasetRecord {
    /// Creates a record without properties.
    ///
    /// # Errors
    ///
    /// Returns [`Error::InvalidInput`] for an empty identifier.
    pub fn new(id: impl Into<String>, structure: Structure) -> Result<Self, Error> {
        let id = id.into();
        if id.trim().is_empty() {
            return Err(Error::InvalidInput(
                "dataset record identifiers must not be empty".to_owned(),
            ));
        }
        Ok(Self {
            id,
            structure,
            properties: BTreeMap::new(),
        })
    }

    /// Adds one validated property.
    ///
    /// # Errors
    ///
    /// Returns [`Error::InvalidInput`] for an empty key or a non-finite value.
    pub fn with_property(
        mut self,
        key: impl Into<String>,
        value: impl Into<PropertyValue>,
    ) -> Result<Self, Error> {
        let key = key.into();
        if key.trim().is_empty() {
            return Err(Error::InvalidInput(
                "dataset property keys must not be empty".to_owned(),
            ));
        }
        let value = value.into();
        value.validate()?;
        self.properties.insert(key, value);
        Ok(self)
    }

    /// Property by name.
    #[must_use]
    pub fn property(&self, key: &str) -> Option<&PropertyValue> {
        self.properties.get(key)
    }

    /// Numeric property by name; integers are converted to `f64`.
    #[must_use]
    pub fn number(&self, key: &str) -> Option<f64> {
        self.property(key).and_then(PropertyValue::as_f64)
    }

    /// Checks the identifier, the property keys and every property value.
    ///
    /// # Errors
    ///
    /// Returns [`Error::InvalidInput`] for an empty identifier or property key,
    /// or a non-finite number.
    pub fn validate(&self) -> Result<(), Error> {
        if self.id.trim().is_empty() {
            return Err(Error::InvalidInput(
                "dataset record identifiers must not be empty".to_owned(),
            ));
        }
        for (key, value) in &self.properties {
            if key.trim().is_empty() {
                return Err(Error::InvalidInput(
                    "dataset property keys must not be empty".to_owned(),
                ));
            }
            value.validate()?;
        }
        Ok(())
    }
}

/// Free-form metadata of a dataset.
#[derive(Debug, Clone, Default, PartialEq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct DatasetMetadata {
    /// Short dataset name.
    pub name: String,
    /// Human-readable description.
    #[serde(default)]
    pub description: String,
    /// Where the records come from.
    #[serde(default)]
    pub source: String,
    /// Arbitrary JSON labels; keys are sorted on write.
    #[serde(default)]
    pub labels: BTreeMap<String, Value>,
}

impl DatasetMetadata {
    /// Adds one arbitrary JSON label.
    ///
    /// # Errors
    ///
    /// Returns [`Error::InvalidInput`] for an empty key.
    pub fn with_label(
        mut self,
        key: impl Into<String>,
        value: impl Into<Value>,
    ) -> Result<Self, Error> {
        let key = key.into();
        if key.trim().is_empty() {
            return Err(Error::InvalidInput(
                "dataset metadata keys must not be empty".to_owned(),
            ));
        }
        self.labels.insert(key, value.into());
        Ok(self)
    }
}
