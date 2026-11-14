//! Typed access to merged preset settings.

use std::collections::BTreeMap;

use serde_json::{Map, Value};

use crate::Error;

/// Merged settings with the preset name used in diagnostics.
#[derive(Debug, Clone)]
pub(crate) struct Settings {
    preset: String,
    values: BTreeMap<String, Value>,
}

impl Settings {
    pub(crate) fn new(preset: String, values: BTreeMap<String, Value>) -> Self {
        Self { preset, values }
    }

    pub(crate) fn raw(&self, path: &str) -> Option<&Value> {
        self.value(path)
    }

    fn value(&self, path: &str) -> Option<&Value> {
        let mut parts = path.split('.');
        let first = parts.next()?;
        let mut current = self.values.get(first)?;
        for part in parts {
            current = current.as_object()?.get(part)?;
        }
        Some(current)
    }

    pub(crate) fn text(&self, path: &str) -> Result<&str, Error> {
        match self.value(path) {
            Some(Value::String(text)) => Ok(text),
            Some(_) => Err(self.invalid(path, "a string")),
            None => Err(self.missing(path)),
        }
    }

    pub(crate) fn optional_text(&self, path: &str) -> Result<Option<&str>, Error> {
        match self.value(path) {
            None | Some(Value::Null) => Ok(None),
            Some(Value::String(text)) if text.is_empty() => Ok(None),
            Some(Value::String(text)) => Ok(Some(text)),
            Some(_) => Err(self.invalid(path, "a string or null")),
        }
    }

    pub(crate) fn float(&self, path: &str) -> Result<f64, Error> {
        match self.value(path) {
            Some(Value::Number(number)) => number
                .as_f64()
                .ok_or_else(|| self.invalid(path, "a finite number")),
            Some(_) => Err(self.invalid(path, "a number")),
            None => Err(self.missing(path)),
        }
    }

    pub(crate) fn flag(&self, path: &str) -> Result<bool, Error> {
        match self.value(path) {
            Some(Value::Bool(value)) => Ok(*value),
            Some(_) => Err(self.invalid(path, "a boolean")),
            None => Err(self.missing(path)),
        }
    }

    pub(crate) fn integer(&self, path: &str) -> Result<i64, Error> {
        match self.value(path) {
            Some(Value::Number(number)) => {
                if let Some(value) = number.as_i64() {
                    return Ok(value);
                }
                let value = number
                    .as_f64()
                    .filter(|value| value.fract().abs() < 1e-9)
                    .ok_or_else(|| self.invalid(path, "an integer"))?;
                Ok(value as i64)
            }
            Some(_) => Err(self.invalid(path, "an integer")),
            None => Err(self.missing(path)),
        }
    }

    pub(crate) fn usize(&self, path: &str) -> Result<usize, Error> {
        let value = self.integer(path)?;
        usize::try_from(value).map_err(|_| self.invalid(path, "a non-negative integer"))
    }

    pub(crate) fn table(&self, path: &str) -> Result<&Map<String, Value>, Error> {
        match self.value(path) {
            Some(Value::Object(table)) => Ok(table),
            Some(_) => Err(self.invalid(path, "an object")),
            None => Err(self.missing(path)),
        }
    }

    pub(crate) fn array(&self, path: &str) -> Result<&[Value], Error> {
        match self.value(path) {
            Some(Value::Array(values)) => Ok(values),
            Some(_) => Err(self.invalid(path, "an array")),
            None => Err(self.missing(path)),
        }
    }

    pub(crate) fn integer_array(&self, path: &str) -> Result<Vec<i64>, Error> {
        self.array(path)?
            .iter()
            .map(|value| {
                value
                    .as_i64()
                    .ok_or_else(|| self.invalid(path, "an array of integers"))
            })
            .collect()
    }

    fn missing(&self, path: &str) -> Error {
        Error::Job(format!("input set {}: missing setting {path}", self.preset))
    }

    fn invalid(&self, path: &str, expected: &str) -> Error {
        Error::Job(format!(
            "input set {}: setting {path} must be {expected}",
            self.preset
        ))
    }
}
