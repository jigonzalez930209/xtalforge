//! Deterministic nested override merging and unknown-key validation.

use std::collections::BTreeMap;

use serde_json::{Map, Value};

use crate::Error;

/// Merges `overrides` into `base` and validates every key.
///
/// The reserved key `allow_unknown` must be a boolean. Every other key is
/// checked against the base tree and the preset's open branches. The result is
/// a [`BTreeMap`], so output order never depends on hash-map iteration.
pub(crate) fn merge(
    base: &Map<String, Value>,
    overrides: &Map<String, Value>,
    open: &[String],
) -> Result<BTreeMap<String, Value>, Error> {
    let allow_unknown = match overrides.get("allow_unknown") {
        None => false,
        Some(Value::Bool(value)) => *value,
        Some(_) => {
            return Err(Error::Job(
                "override allow_unknown must be a boolean".to_owned(),
            ));
        }
    };
    validate(base, overrides, "", open, allow_unknown)?;
    let mut merged = BTreeMap::new();
    for (key, value) in base {
        merged.insert(key.clone(), value.clone());
    }
    for (key, value) in overrides {
        if key == "allow_unknown" {
            continue;
        }
        match (merged.get(key), value) {
            (Some(Value::Object(base_child)), Value::Object(override_child)) => {
                let child = deep_merge(base_child, override_child);
                merged.insert(key.clone(), Value::Object(child.into_iter().collect()));
            }
            _ => {
                merged.insert(key.clone(), value.clone());
            }
        }
    }
    Ok(merged)
}

/// Deep-merges two objects without validation.
pub(crate) fn deep_merge(
    base: &Map<String, Value>,
    overrides: &Map<String, Value>,
) -> BTreeMap<String, Value> {
    let mut merged = BTreeMap::new();
    for (key, value) in base {
        merged.insert(key.clone(), value.clone());
    }
    for (key, value) in overrides {
        match (merged.get(key), value) {
            (Some(Value::Object(base_child)), Value::Object(override_child)) => {
                let child = deep_merge(base_child, override_child);
                merged.insert(key.clone(), Value::Object(child.into_iter().collect()));
            }
            _ => {
                merged.insert(key.clone(), value.clone());
            }
        }
    }
    merged
}

fn validate(
    base: &Map<String, Value>,
    overrides: &Map<String, Value>,
    prefix: &str,
    open: &[String],
    allow_unknown: bool,
) -> Result<(), Error> {
    for (key, value) in overrides {
        if key == "allow_unknown" && prefix.is_empty() {
            continue;
        }
        let path = if prefix.is_empty() {
            key.clone()
        } else {
            format!("{prefix}.{key}")
        };
        let known = base.contains_key(key) || is_open(&path, open);
        if !known && !allow_unknown {
            return Err(Error::UnknownOverride(path));
        }
        if let (Some(Value::Object(base_child)), Value::Object(override_child)) =
            (base.get(key), value)
        {
            validate(base_child, override_child, &path, open, allow_unknown)?;
        }
    }
    Ok(())
}

fn is_open(path: &str, open: &[String]) -> bool {
    open.iter().any(|branch| {
        path == branch
            || path
                .strip_prefix(branch.as_str())
                .is_some_and(|rest| rest.starts_with('.'))
    })
}
