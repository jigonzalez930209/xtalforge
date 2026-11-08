//! Typed parameter descriptions shared by modules, engines and builders
//! (Task 3.8.1).

use std::collections::BTreeMap;

use serde::{Deserialize, Serialize};
use serde_json::Value;

use crate::Error;

/// The type and constraints of a parameter.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
#[serde(tag = "type", rename_all = "snake_case")]
pub enum ParamType {
    /// A whole number with optional bounds.
    #[serde(rename = "integer")]
    Integer {
        /// Minimum, inclusive.
        min: Option<i64>,
        /// Maximum, inclusive.
        max: Option<i64>,
    },
    /// A floating-point number with optional bounds and unit.
    #[serde(rename = "number")]
    Float {
        /// Minimum, inclusive.
        min: Option<f64>,
        /// Maximum, inclusive.
        max: Option<f64>,
        /// Unit shown to the user, for example `Angstrom`.
        unit: Option<String>,
    },
    /// A boolean switch.
    #[serde(rename = "boolean")]
    Boolean,
    /// One of a fixed set of strings.
    Choice {
        /// Allowed values.
        choices: Vec<String>,
    },
    /// Free text.
    #[serde(rename = "string")]
    Text,
    /// Three floating-point numbers.
    #[serde(rename = "array")]
    Vector3 {
        /// Unit shown to the user.
        unit: Option<String>,
    },
    /// A filesystem path (a string on the browser).
    Path {
        /// Whether the path must already exist.
        must_exist: bool,
    },
}

/// One described parameter.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct ParamSpec {
    /// Stable identifier, `snake_case`.
    pub name: String,
    /// Type and constraints.
    pub kind: ParamType,
    /// Default value, already valid for the type.
    pub default: Value,
    /// Help text.
    pub help: String,
    /// Whether the caller must provide it.
    pub required: bool,
    /// Optional grouping label for a user interface.
    pub group: Option<String>,
}

impl ParamSpec {
    /// Validates a value against this specification.
    ///
    /// # Errors
    ///
    /// Returns [`Error::InvalidInput`] with a message naming the parameter.
    pub fn validate(&self, value: &Value) -> Result<(), Error> {
        let fail = |reason: &str| {
            Err(Error::InvalidInput(format!(
                "parameter {}: {reason}",
                self.name
            )))
        };
        match &self.kind {
            ParamType::Integer { min, max } => {
                let Some(number) = value.as_i64() else {
                    return fail("expected an integer");
                };
                if min.is_some_and(|min| number < min) {
                    return fail("below the minimum");
                }
                if max.is_some_and(|max| number > max) {
                    return fail("above the maximum");
                }
            }
            ParamType::Float { min, max, .. } => {
                let Some(number) = value.as_f64() else {
                    return fail("expected a number");
                };
                if min.is_some_and(|min| number < min) {
                    return fail("below the minimum");
                }
                if max.is_some_and(|max| number > max) {
                    return fail("above the maximum");
                }
            }
            ParamType::Boolean => {
                if !value.is_boolean() {
                    return fail("expected true or false");
                }
            }
            ParamType::Choice { choices } => {
                let Some(text) = value.as_str() else {
                    return fail("expected a choice string");
                };
                if !choices.iter().any(|choice| choice == text) {
                    return fail("not one of the allowed choices");
                }
            }
            ParamType::Text => {
                if !value.is_string() {
                    return fail("expected text");
                }
            }
            ParamType::Vector3 { .. } => {
                let Some(values) = value.as_array() else {
                    return fail("expected three numbers");
                };
                if values.len() != 3 || !values.iter().all(Value::is_number) {
                    return fail("expected three numbers");
                }
            }
            ParamType::Path { .. } => {
                if !value.is_string() {
                    return fail("expected a path string");
                }
            }
        }
        Ok(())
    }
}

/// A set of parameter descriptions and their current values.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct ParamSet {
    /// Descriptions in declaration order.
    pub specs: Vec<ParamSpec>,
    /// Current values; missing entries use the default.
    pub values: BTreeMap<String, Value>,
}

impl ParamSet {
    /// Creates an empty set.
    #[must_use]
    pub fn new() -> Self {
        Self {
            specs: Vec::new(),
            values: BTreeMap::new(),
        }
    }

    /// Adds a specification and its default value.
    ///
    /// # Errors
    ///
    /// Returns [`Error::InvalidInput`] for a duplicate name or a default that
    /// violates the specification.
    pub fn with_param(mut self, spec: ParamSpec) -> Result<Self, Error> {
        if self.specs.iter().any(|existing| existing.name == spec.name) {
            return Err(Error::InvalidInput(format!(
                "duplicate parameter {}",
                spec.name
            )));
        }
        spec.validate(&spec.default)?;
        self.values.insert(spec.name.clone(), spec.default.clone());
        self.specs.push(spec);
        Ok(self)
    }

    /// Sets a value after validation.
    ///
    /// # Errors
    ///
    /// Returns [`Error::InvalidInput`] for an unknown parameter or an invalid
    /// value.
    pub fn set(&mut self, name: &str, value: Value) -> Result<(), Error> {
        let spec = self
            .specs
            .iter()
            .find(|spec| spec.name == name)
            .ok_or_else(|| Error::InvalidInput(format!("unknown parameter {name}")))?;
        spec.validate(&value)?;
        self.values.insert(name.to_owned(), value);
        Ok(())
    }

    /// Gets a value.
    #[must_use]
    pub fn get(&self, name: &str) -> Option<&Value> {
        self.values.get(name)
    }

    /// Checks that every required parameter has a value.
    ///
    /// # Errors
    ///
    /// Returns [`Error::InvalidInput`] naming the missing parameter.
    pub fn check_required(&self) -> Result<(), Error> {
        for spec in &self.specs {
            if spec.required && !self.values.contains_key(&spec.name) {
                return Err(Error::InvalidInput(format!(
                    "missing required parameter {}",
                    spec.name
                )));
            }
        }
        Ok(())
    }

    /// A JSON schema-like object for the agent layer.
    #[must_use]
    pub fn schema(&self) -> Value {
        let properties: BTreeMap<&str, Value> = self
            .specs
            .iter()
            .map(|spec| (spec.name.as_str(), spec.kind.clone().into_value(spec)))
            .collect();
        serde_json::json!({
            "type": "object",
            "properties": properties,
            "required": self.specs.iter().filter(|spec| spec.required).map(|spec| spec.name.as_str()).collect::<Vec<_>>(),
        })
    }
}

impl Default for ParamSet {
    fn default() -> Self {
        Self::new()
    }
}

impl ParamType {
    fn into_value(self, spec: &ParamSpec) -> Value {
        let mut object = serde_json::to_value(&self).unwrap_or(Value::Null);
        if let Some(map) = object.as_object_mut() {
            map.insert("default".to_owned(), spec.default.clone());
            map.insert("description".to_owned(), Value::String(spec.help.clone()));
            if let Some(group) = &spec.group {
                map.insert("group".to_owned(), Value::String(group.clone()));
            }
        }
        object
    }
}

impl ParamSet {
    /// Serializes the descriptions and current values.
    ///
    /// # Errors
    ///
    /// Returns [`Error::InvalidInput`] when serialization fails.
    pub fn to_json(&self) -> Result<String, Error> {
        serde_json::to_string_pretty(self).map_err(|error| {
            Error::InvalidInput(format!("could not serialize parameters: {error}"))
        })
    }

    /// Parses a stored parameter set and validates every value.
    ///
    /// # Errors
    ///
    /// Returns [`Error::InvalidInput`] for invalid JSON, unknown fields or a
    /// value that violates its specification.
    pub fn from_json(text: &str) -> Result<Self, Error> {
        let set: Self = serde_json::from_str(text)
            .map_err(|error| Error::InvalidInput(format!("invalid parameters: {error}")))?;
        for (name, value) in &set.values {
            let spec = set
                .specs
                .iter()
                .find(|spec| spec.name == *name)
                .ok_or_else(|| Error::InvalidInput(format!("unknown parameter {name}")))?;
            spec.validate(value)?;
        }
        Ok(set)
    }
}
