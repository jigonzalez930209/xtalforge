//! JSON Lines dataset container (Task 16.3.1).
//!
//! A [`Dataset`] holds many structures and their properties in one file. The
//! file is JSON Lines: the first line is a versioned header and every other
//! line is one [`DatasetRecord`].
//!
//! ```text
//! {"dataset_version":1,"n_records":2,"metadata":{...}}
//! {"id":"frame-1","structure":{...},"properties":{"energy":-1.5}}
//! {"id":"frame-2","structure":{...},"properties":{"energy":-2.0}}
//! ```
//!
//! The writer preserves the record order and sorts every map key, so the same
//! dataset always produces byte-identical output. The reader validates the
//! header, the declared record count, the property types and the uniqueness of
//! the record identifiers.
//!
//! JSON numbers are parsed by `serde_json` with its default float parser, which
//! can differ from the written value in the last bit. Writing an in-memory
//! dataset is deterministic; a read/write round trip is bit-exact only when
//! the `serde_json` dependency enables its `float_roundtrip` feature.

use std::collections::{BTreeMap, BTreeSet};
#[cfg(not(target_arch = "wasm32"))]
use std::path::Path;

use serde::{Deserialize, Serialize};
use serde_json::Value;
use xtalforge_core::Error;

mod record;

pub use record::{DatasetMetadata, DatasetRecord, PropertyValue};

/// Version of the JSON Lines container written by this crate.
pub const DATASET_VERSION: u32 = 1;

#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
struct DatasetHeader {
    dataset_version: u32,
    n_records: usize,
    metadata: DatasetMetadata,
}

/// A deterministic container of structures and properties.
#[derive(Debug, Clone, PartialEq)]
pub struct Dataset {
    metadata: DatasetMetadata,
    records: Vec<DatasetRecord>,
    index: BTreeMap<String, usize>,
}

impl Dataset {
    /// Creates an empty dataset with a name.
    #[must_use]
    pub fn new(name: impl Into<String>) -> Self {
        Self {
            metadata: DatasetMetadata {
                name: name.into(),
                ..DatasetMetadata::default()
            },
            records: Vec::new(),
            index: BTreeMap::new(),
        }
    }

    /// Builds a dataset from ordered records, rejecting duplicate identifiers.
    ///
    /// # Errors
    ///
    /// Returns [`Error::InvalidInput`] for an invalid or duplicate record.
    pub fn from_records(
        name: impl Into<String>,
        records: Vec<DatasetRecord>,
    ) -> Result<Self, Error> {
        let mut dataset = Self::new(name);
        for record in records {
            dataset.push(record)?;
        }
        Ok(dataset)
    }

    /// Dataset metadata.
    #[must_use]
    pub fn metadata(&self) -> &DatasetMetadata {
        &self.metadata
    }

    /// Mutable dataset metadata. Changing the name does not change the
    /// records.
    pub fn metadata_mut(&mut self) -> &mut DatasetMetadata {
        &mut self.metadata
    }

    /// Sets the description.
    #[must_use]
    pub fn with_description(mut self, description: impl Into<String>) -> Self {
        self.metadata.description = description.into();
        self
    }

    /// Sets the source.
    #[must_use]
    pub fn with_source(mut self, source: impl Into<String>) -> Self {
        self.metadata.source = source.into();
        self
    }

    /// Adds one metadata label.
    ///
    /// # Errors
    ///
    /// Returns [`Error::InvalidInput`] for an empty key.
    pub fn with_label(
        mut self,
        key: impl Into<String>,
        value: impl Into<Value>,
    ) -> Result<Self, Error> {
        self.metadata = self.metadata.with_label(key, value)?;
        Ok(self)
    }

    /// Number of records.
    #[must_use]
    pub fn len(&self) -> usize {
        self.records.len()
    }

    /// Whether the dataset has no records.
    #[must_use]
    pub fn is_empty(&self) -> bool {
        self.records.is_empty()
    }

    /// All records in insertion order.
    #[must_use]
    pub fn records(&self) -> &[DatasetRecord] {
        &self.records
    }

    /// Record at an index.
    #[must_use]
    pub fn get(&self, index: usize) -> Option<&DatasetRecord> {
        self.records.get(index)
    }

    /// Record by identifier.
    #[must_use]
    pub fn by_id(&self, id: &str) -> Option<&DatasetRecord> {
        self.index
            .get(id)
            .and_then(|index| self.records.get(*index))
    }

    /// Iterates over the records in insertion order.
    pub fn iter(&self) -> std::slice::Iter<'_, DatasetRecord> {
        self.records.iter()
    }

    /// Appends a record.
    ///
    /// # Errors
    ///
    /// Returns [`Error::InvalidInput`] for an invalid record or a duplicate
    /// identifier.
    pub fn push(&mut self, record: DatasetRecord) -> Result<(), Error> {
        record.validate()?;
        if self.index.contains_key(&record.id) {
            return Err(Error::InvalidInput(format!(
                "duplicate dataset record identifier {:?}",
                record.id
            )));
        }
        let index = self.records.len();
        self.index.insert(record.id.clone(), index);
        self.records.push(record);
        Ok(())
    }

    /// Appends a record and returns the dataset.
    ///
    /// # Errors
    ///
    /// Returns the errors of [`Dataset::push`].
    pub fn with_record(mut self, record: DatasetRecord) -> Result<Self, Error> {
        self.push(record)?;
        Ok(self)
    }

    /// Sorted union of every property key in the dataset.
    #[must_use]
    pub fn property_names(&self) -> Vec<String> {
        let mut names = BTreeSet::new();
        for record in &self.records {
            names.extend(record.properties.keys().cloned());
        }
        names.into_iter().collect()
    }

    /// One property column in record order; missing properties are `None`.
    #[must_use]
    pub fn column(&self, name: &str) -> Vec<Option<&PropertyValue>> {
        self.records
            .iter()
            .map(|record| record.property(name))
            .collect()
    }

    /// Serializes the dataset as JSON Lines with a trailing newline.
    ///
    /// The header carries the declared record count, so a truncated file is
    /// detected on read.
    ///
    /// # Errors
    ///
    /// Returns [`Error::InvalidInput`] for an invalid record or a serialization
    /// failure.
    pub fn to_json_lines(&self) -> Result<String, Error> {
        for record in &self.records {
            record.validate()?;
        }
        let header = DatasetHeader {
            dataset_version: DATASET_VERSION,
            n_records: self.records.len(),
            metadata: self.metadata.clone(),
        };
        let mut text = serde_json::to_string(&header)
            .map_err(|error| Error::InvalidInput(format!("dataset header encoding: {error}")))?;
        text.push('\n');
        for record in &self.records {
            let line = serde_json::to_string(record).map_err(|error| {
                Error::InvalidInput(format!("dataset record encoding: {error}"))
            })?;
            text.push_str(&line);
            text.push('\n');
        }
        Ok(text)
    }

    /// Parses a JSON Lines dataset.
    ///
    /// Blank lines are ignored. Errors report the physical line number.
    ///
    /// # Errors
    ///
    /// Returns [`Error::InvalidInput`] for an empty input, a malformed header
    /// or record, an unsupported version, a duplicate identifier, or a record
    /// count that differs from the header.
    pub fn from_json_lines(text: &str) -> Result<Self, Error> {
        let mut header: Option<DatasetHeader> = None;
        let mut dataset = Self::new("");
        let mut declared = 0;
        for (line_index, line) in text.lines().enumerate() {
            let line_number = line_index + 1;
            if line.trim().is_empty() {
                continue;
            }
            if header.is_none() {
                let parsed: DatasetHeader = serde_json::from_str(line).map_err(|error| {
                    Error::InvalidInput(format!(
                        "dataset line {line_number}: invalid header: {error}"
                    ))
                })?;
                if parsed.dataset_version != DATASET_VERSION {
                    return Err(Error::InvalidInput(format!(
                        "dataset line {line_number}: unsupported dataset_version {}, expected {DATASET_VERSION}",
                        parsed.dataset_version
                    )));
                }
                declared = parsed.n_records;
                dataset.metadata = parsed.metadata.clone();
                header = Some(parsed);
                continue;
            }
            let record: DatasetRecord = serde_json::from_str(line).map_err(|error| {
                Error::InvalidInput(format!(
                    "dataset line {line_number}: invalid record: {error}"
                ))
            })?;
            dataset.push(record).map_err(|error| {
                Error::InvalidInput(format!("dataset line {line_number}: {error}"))
            })?;
        }
        if header.is_none() {
            return Err(Error::InvalidInput(
                "empty dataset: expected a JSON Lines header".to_owned(),
            ));
        }
        if dataset.records.len() != declared {
            return Err(Error::InvalidInput(format!(
                "dataset header declares {declared} records but the file contains {}",
                dataset.records.len()
            )));
        }
        Ok(dataset)
    }

    /// Writes the dataset to a JSON Lines file.
    ///
    /// # Errors
    ///
    /// Returns [`Error::InvalidInput`] when the dataset cannot be encoded or
    /// the file cannot be written.
    #[cfg(not(target_arch = "wasm32"))]
    pub fn save(&self, path: impl AsRef<Path>) -> Result<(), Error> {
        let path = path.as_ref();
        std::fs::write(path, self.to_json_lines()?).map_err(|error| {
            Error::InvalidInput(format!("could not write {}: {error}", path.display()))
        })
    }

    /// Reads a dataset from a JSON Lines file.
    ///
    /// # Errors
    ///
    /// Returns [`Error::InvalidInput`] when the file cannot be read or parsed.
    #[cfg(not(target_arch = "wasm32"))]
    pub fn load(path: impl AsRef<Path>) -> Result<Self, Error> {
        let path = path.as_ref();
        let text = std::fs::read_to_string(path).map_err(|error| {
            Error::InvalidInput(format!("could not read {}: {error}", path.display()))
        })?;
        Self::from_json_lines(&text)
    }
}

impl<'a> IntoIterator for &'a Dataset {
    type Item = &'a DatasetRecord;
    type IntoIter = std::slice::Iter<'a, DatasetRecord>;

    fn into_iter(self) -> Self::IntoIter {
        self.records.iter()
    }
}

impl IntoIterator for Dataset {
    type Item = DatasetRecord;
    type IntoIter = std::vec::IntoIter<DatasetRecord>;

    fn into_iter(self) -> Self::IntoIter {
        self.records.into_iter()
    }
}
