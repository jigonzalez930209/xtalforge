//! Provenance records for reproducible results (Task 16.2.1).
//!
//! Every operation can produce a [`Provenance`] record naming the producing
//! version, the operation and its options, the SHA-256 digest of the inputs,
//! the data bundles used, the random seed and an optional timestamp. The
//! record serializes to deterministic JSON: fields appear in declaration
//! order and every map is a [`BTreeMap`], so two records with the same values
//! always produce the same text regardless of insertion order.
//!
//! The schema has its own version, [`PROVENANCE_SCHEMA_VERSION`], independent
//! of the crate version (see `AGENTS.md`, section 5.2). Readers reject a
//! different version instead of guessing. Records are attached to structures
//! through the [`PROVENANCE_METADATA_KEY`] metadata entry, which holds the
//! canonical JSON string; the format readers in `xtalforge-io` embed that same
//! string in CIF, extended XYZ and project files.

use std::collections::BTreeMap;

use serde::{Deserialize, Serialize};
use serde_json::Value;
use sha2::{Digest, Sha256};

use crate::{Error, Structure};

/// Version of the provenance schema, independent of the crate version.
pub const PROVENANCE_SCHEMA_VERSION: u32 = 1;

/// Structure metadata key that carries a serialized [`Provenance`].
pub const PROVENANCE_METADATA_KEY: &str = "provenance";

/// Lowercase hexadecimal SHA-256 digest of `bytes`.
#[must_use]
pub fn sha256_hex(bytes: &[u8]) -> String {
    let digest = Sha256::digest(bytes);
    let mut text = String::with_capacity(64);
    for byte in digest {
        text.push_str(&format!("{byte:02x}"));
    }
    text
}

/// One input of an operation and its SHA-256 digest.
#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct ProvenanceInput {
    /// Optional logical name: a path, URL or dataset identifier.
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub name: Option<String>,
    /// Lowercase hexadecimal SHA-256 digest of the input bytes.
    pub sha256: String,
}

impl ProvenanceInput {
    /// Builds an input from a digest, rejecting a malformed digest.
    ///
    /// # Errors
    ///
    /// Returns [`Error::InvalidInput`] when `sha256` is not 64 hexadecimal
    /// characters.
    pub fn new(name: Option<String>, sha256: impl Into<String>) -> Result<Self, Error> {
        let input = Self {
            name,
            sha256: sha256.into(),
        };
        input.validate()?;
        Ok(input)
    }

    /// Builds an input by hashing `bytes`.
    #[must_use]
    pub fn from_bytes(name: Option<String>, bytes: &[u8]) -> Self {
        Self {
            name,
            sha256: sha256_hex(bytes),
        }
    }

    /// Checks that the digest is lowercase or uppercase hexadecimal.
    ///
    /// # Errors
    ///
    /// Returns [`Error::InvalidInput`] when the digest has the wrong length or
    /// contains a non-hexadecimal character.
    pub fn validate(&self) -> Result<(), Error> {
        if is_sha256_hex(&self.sha256) {
            Ok(())
        } else {
            Err(Error::InvalidInput(format!(
                "provenance: input digest must be 64 hexadecimal characters, got {:?}",
                self.sha256
            )))
        }
    }
}

/// Record of how a result was produced.
///
/// The JSON produced by [`Provenance::to_json`] is the interchange form used
/// by the CLI, Python, WASM, CIF, extended XYZ and project files.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct Provenance {
    /// Schema version, independent of the crate version.
    pub schema_version: u32,
    /// Version of xtalforge that produced the result.
    pub xtalforge_version: String,
    /// Operation name: a verb, command or workflow step.
    pub operation: String,
    /// Operation options, keyed by name.
    #[serde(default, skip_serializing_if = "BTreeMap::is_empty")]
    pub options: BTreeMap<String, Value>,
    /// Inputs and their SHA-256 digests, in declaration order.
    #[serde(default, skip_serializing_if = "Vec::is_empty")]
    pub inputs: Vec<ProvenanceInput>,
    /// Data bundles used, mapped to their version or digest.
    #[serde(default, skip_serializing_if = "BTreeMap::is_empty")]
    pub bundles: BTreeMap<String, String>,
    /// Random seed, when the operation is stochastic.
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub seed: Option<u64>,
    /// ISO 8601 timestamp, when the caller records one.
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub timestamp: Option<String>,
}

impl Provenance {
    /// Starts a record for `operation` with the current xtalforge version.
    #[must_use]
    pub fn new(operation: impl Into<String>) -> Self {
        Self {
            schema_version: PROVENANCE_SCHEMA_VERSION,
            xtalforge_version: crate::version().to_owned(),
            operation: operation.into(),
            options: BTreeMap::new(),
            inputs: Vec::new(),
            bundles: BTreeMap::new(),
            seed: None,
            timestamp: None,
        }
    }

    /// Adds or replaces one operation option.
    #[must_use]
    pub fn with_option(mut self, key: impl Into<String>, value: impl Into<Value>) -> Self {
        self.options.insert(key.into(), value.into());
        self
    }

    /// Adds every option from `options`.
    #[must_use]
    pub fn with_options(mut self, options: BTreeMap<String, Value>) -> Self {
        self.options.extend(options);
        self
    }

    /// Adds one input digest.
    #[must_use]
    pub fn with_input(mut self, input: ProvenanceInput) -> Self {
        self.inputs.push(input);
        self
    }

    /// Adds one input by hashing its bytes.
    #[must_use]
    pub fn with_input_bytes(mut self, name: Option<String>, bytes: &[u8]) -> Self {
        self.inputs.push(ProvenanceInput::from_bytes(name, bytes));
        self
    }

    /// Adds one input digest after validating it.
    ///
    /// # Errors
    ///
    /// Returns [`Error::InvalidInput`] when the digest is malformed.
    pub fn try_with_input_sha256(
        self,
        name: Option<String>,
        sha256: impl Into<String>,
    ) -> Result<Self, Error> {
        Ok(self.with_input(ProvenanceInput::new(name, sha256)?))
    }

    /// Records a data bundle and its version or digest.
    #[must_use]
    pub fn with_bundle(mut self, name: impl Into<String>, version: impl Into<String>) -> Self {
        self.bundles.insert(name.into(), version.into());
        self
    }

    /// Records the random seed.
    #[must_use]
    pub fn with_seed(mut self, seed: u64) -> Self {
        self.seed = Some(seed);
        self
    }

    /// Records a timestamp. The caller owns the format; ISO 8601 UTC is
    /// recommended so records sort chronologically.
    #[must_use]
    pub fn with_timestamp(mut self, timestamp: impl Into<String>) -> Self {
        self.timestamp = Some(timestamp.into());
        self
    }

    /// Checks the schema version, required strings, digests and bundle keys.
    ///
    /// # Errors
    ///
    /// Returns [`Error::InvalidInput`] describing the first invalid field.
    pub fn validate(&self) -> Result<(), Error> {
        if self.schema_version != PROVENANCE_SCHEMA_VERSION {
            return Err(Error::InvalidInput(format!(
                "provenance: unsupported schema_version {}, expected {PROVENANCE_SCHEMA_VERSION}",
                self.schema_version
            )));
        }
        if self.xtalforge_version.trim().is_empty() {
            return Err(Error::InvalidInput(
                "provenance: xtalforge_version must not be empty".to_owned(),
            ));
        }
        if self.operation.trim().is_empty() {
            return Err(Error::InvalidInput(
                "provenance: operation must not be empty".to_owned(),
            ));
        }
        for input in &self.inputs {
            input.validate()?;
        }
        for (name, version) in &self.bundles {
            if name.trim().is_empty() || version.trim().is_empty() {
                return Err(Error::InvalidInput(
                    "provenance: bundle names and versions must not be empty".to_owned(),
                ));
            }
        }
        if self
            .timestamp
            .as_deref()
            .is_some_and(|timestamp| timestamp.trim().is_empty())
        {
            return Err(Error::InvalidInput(
                "provenance: timestamp must not be empty when present".to_owned(),
            ));
        }
        Ok(())
    }

    /// Serializes the record to deterministic JSON.
    ///
    /// Fields appear in declaration order, maps are sorted by key and absent
    /// optional fields are omitted.
    ///
    /// # Errors
    ///
    /// Returns [`Error::InvalidInput`] when validation or serialization fails.
    pub fn to_json(&self) -> Result<String, Error> {
        self.validate()?;
        serde_json::to_string(self).map_err(|error| {
            Error::InvalidInput(format!("provenance: could not serialize record: {error}"))
        })
    }

    /// Parses a record and checks its schema version.
    ///
    /// # Errors
    ///
    /// Returns [`Error::InvalidInput`] for invalid JSON, unknown fields, a
    /// newer schema version or an invalid field value.
    pub fn from_json(text: &str) -> Result<Self, Error> {
        let record: Self = serde_json::from_str(text).map_err(|error| {
            Error::InvalidInput(format!("provenance: invalid JSON record: {error}"))
        })?;
        record.validate()?;
        Ok(record)
    }

    /// Parses a record stored as a JSON string or a JSON object.
    ///
    /// # Errors
    ///
    /// Returns [`Error::InvalidInput`] for any other JSON value or an invalid
    /// record.
    pub fn from_value(value: &Value) -> Result<Self, Error> {
        match value {
            Value::String(text) => Self::from_json(text),
            Value::Object(_) => {
                let record: Self = serde_json::from_value(value.clone()).map_err(|error| {
                    Error::InvalidInput(format!("provenance: invalid JSON record: {error}"))
                })?;
                record.validate()?;
                Ok(record)
            }
            _ => Err(Error::InvalidInput(
                "provenance: metadata must be a JSON string or object".to_owned(),
            )),
        }
    }

    /// Attaches this record to a copy of `structure`.
    ///
    /// The canonical JSON replaces any previous provenance entry in
    /// [`PROVENANCE_METADATA_KEY`].
    ///
    /// # Errors
    ///
    /// Returns the errors of [`Provenance::to_json`].
    pub fn attach(&self, structure: &Structure) -> Result<Structure, Error> {
        let json = self.to_json()?;
        Ok(structure
            .clone()
            .with_metadata(PROVENANCE_METADATA_KEY, Value::String(json)))
    }

    /// Reads the record attached to `structure`, if any.
    ///
    /// # Errors
    ///
    /// Returns [`Error::InvalidInput`] when the metadata entry exists but is
    /// not a valid provenance record.
    pub fn from_structure(structure: &Structure) -> Result<Option<Self>, Error> {
        structure
            .metadata()
            .get(PROVENANCE_METADATA_KEY)
            .map(Self::from_value)
            .transpose()
    }

    /// Finds the first input whose name equals `name`; `None` matches an
    /// unnamed input.
    #[must_use]
    pub fn input(&self, name: Option<&str>) -> Option<&ProvenanceInput> {
        self.inputs
            .iter()
            .find(|input| input.name.as_deref() == name)
    }

    /// Recomputes the SHA-256 of `bytes` and compares it with a recorded
    /// input.
    ///
    /// # Errors
    ///
    /// Returns [`Error::InvalidInput`] when no input matches `name` or the
    /// digest differs, which detects a modified input.
    pub fn verify_input(&self, name: Option<&str>, bytes: &[u8]) -> Result<(), Error> {
        let input = self.input(name).ok_or_else(|| {
            Error::InvalidInput(format!(
                "provenance: no recorded input named {}",
                name.unwrap_or("<unnamed>")
            ))
        })?;
        input.validate()?;
        let actual = sha256_hex(bytes);
        if actual.eq_ignore_ascii_case(&input.sha256) {
            Ok(())
        } else {
            Err(Error::InvalidInput(format!(
                "provenance: input {} does not match the recorded digest",
                name.unwrap_or("<unnamed>")
            )))
        }
    }
}

fn is_sha256_hex(text: &str) -> bool {
    text.len() == 64 && text.bytes().all(|byte| byte.is_ascii_hexdigit())
}
