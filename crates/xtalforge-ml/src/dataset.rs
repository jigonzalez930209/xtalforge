//! Feature tables and columnar Parquet containers (Tasks 21.1 and 21.2.1).

use serde::{Deserialize, Serialize};
#[cfg(feature = "parquet")]
use std::collections::BTreeMap;
use xtalforge_analysis::dataset::Dataset;
#[cfg(feature = "parquet")]
use xtalforge_core::StructureDocument;
use xtalforge_core::{Error, Structure};

/// A dense feature table ready for a model.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct FeatureTable {
    /// Record identifiers, in dataset order.
    pub ids: Vec<String>,
    /// One feature row per record.
    pub features: Vec<Vec<f64>>,
    /// Optional target row per record.
    pub targets: Vec<Vec<f64>>,
    /// Feature names, when the featurizer supplies them.
    pub feature_names: Vec<String>,
    /// Target names, when the target extractor supplies them.
    pub target_names: Vec<String>,
}

impl FeatureTable {
    /// Number of rows.
    #[must_use]
    pub fn len(&self) -> usize {
        self.ids.len()
    }

    /// Whether the table is empty.
    #[must_use]
    pub fn is_empty(&self) -> bool {
        self.ids.is_empty()
    }

    /// Serializes the table to JSON.
    ///
    /// # Errors
    ///
    /// Returns [`Error::InvalidInput`] when serialization fails.
    pub fn to_json(&self) -> Result<String, Error> {
        serde_json::to_string(self)
            .map_err(|error| Error::InvalidInput(format!("could not serialize table: {error}")))
    }
}

/// Builds a feature table from a dataset and a featurizer.
///
/// `featurizer` returns the feature vector for one structure; `targets` names
/// dataset properties to collect in order.
///
/// # Errors
///
/// Returns the featurizer errors and [`Error::InvalidInput`] when a requested
/// property is not numeric or dimensions differ.
pub fn feature_table<F>(
    dataset: &Dataset,
    targets: &[String],
    mut featurizer: F,
) -> Result<FeatureTable, Error>
where
    F: FnMut(&Structure) -> Result<Vec<f64>, Error>,
{
    let mut ids = Vec::with_capacity(dataset.len());
    let mut features = Vec::with_capacity(dataset.len());
    let mut target_rows = Vec::with_capacity(dataset.len());
    let mut width: Option<usize> = None;
    for record in dataset.iter() {
        let row = featurizer(&record.structure)?;
        if let Some(width) = width
            && row.len() != width
        {
            return Err(Error::InvalidInput(format!(
                "record `{}` has {} features, expected {width}",
                record.id,
                row.len()
            )));
        }
        width = Some(row.len());
        ids.push(record.id.clone());
        features.push(row);
        let mut target_row = Vec::with_capacity(targets.len());
        for target in targets {
            let Some(value) = record.properties.get(target) else {
                return Err(Error::InvalidInput(format!(
                    "record `{}` has no property `{target}`",
                    record.id
                )));
            };
            let number = value.as_f64().ok_or_else(|| {
                Error::InvalidInput(format!(
                    "property `{target}` of `{}` is not numeric",
                    record.id
                ))
            })?;
            target_row.push(number);
        }
        target_rows.push(target_row);
    }
    Ok(FeatureTable {
        ids,
        features,
        targets: target_rows,
        feature_names: Vec::new(),
        target_names: targets.to_vec(),
    })
}

/// Writes a dataset as Parquet bytes.
///
/// # Errors
///
/// Returns [`Error::InvalidInput`] for Arrow/Parquet errors.
#[cfg(feature = "parquet")]
pub fn write_parquet_bytes(dataset: &Dataset) -> Result<Vec<u8>, Error> {
    use std::sync::Arc;

    use arrow::array::StringArray;
    use arrow::datatypes::{DataType, Field, Schema};
    use arrow::record_batch::RecordBatch;
    use parquet::arrow::ArrowWriter;

    let mut ids = Vec::with_capacity(dataset.len());
    let mut structures = Vec::with_capacity(dataset.len());
    let mut properties = Vec::with_capacity(dataset.len());
    for record in dataset.iter() {
        ids.push(record.id.clone());
        structures.push(StructureDocument::new(record.structure.clone()).to_json()?);
        properties.push(
            serde_json::to_string(&record.properties)
                .map_err(|error| Error::InvalidInput(format!("properties JSON: {error}")))?,
        );
    }
    let schema = Arc::new(Schema::new(vec![
        Field::new("id", DataType::Utf8, false),
        Field::new("structure", DataType::Utf8, false),
        Field::new("properties", DataType::Utf8, false),
    ]));
    let batch = RecordBatch::try_new(
        schema.clone(),
        vec![
            Arc::new(StringArray::from(ids)),
            Arc::new(StringArray::from(structures)),
            Arc::new(StringArray::from(properties)),
        ],
    )
    .map_err(|error| Error::InvalidInput(format!("Parquet batch: {error}")))?;
    let mut buffer = Vec::new();
    let mut writer = ArrowWriter::try_new(&mut buffer, schema, None)
        .map_err(|error| Error::InvalidInput(format!("Parquet writer: {error}")))?;
    writer
        .write(&batch)
        .map_err(|error| Error::InvalidInput(format!("Parquet write: {error}")))?;
    writer
        .close()
        .map_err(|error| Error::InvalidInput(format!("Parquet close: {error}")))?;
    Ok(buffer)
}

/// Reads Parquet bytes written by [`write_parquet_bytes`].
///
/// # Errors
///
/// Returns [`Error::InvalidInput`] for Parquet or JSON errors.
#[cfg(feature = "parquet")]
pub fn read_parquet_bytes(bytes: &[u8]) -> Result<Dataset, Error> {
    use parquet::arrow::arrow_reader::ParquetRecordBatchReaderBuilder;

    let data = bytes::Bytes::copy_from_slice(bytes);
    let builder = ParquetRecordBatchReaderBuilder::try_new(data)
        .map_err(|error| Error::InvalidInput(format!("Parquet reader: {error}")))?;
    let reader = builder
        .build()
        .map_err(|error| Error::InvalidInput(format!("Parquet reader: {error}")))?;
    let mut dataset = Dataset::new("parquet");
    for batch in reader {
        let batch =
            batch.map_err(|error| Error::InvalidInput(format!("Parquet batch: {error}")))?;
        let ids = batch
            .column(0)
            .as_any()
            .downcast_ref::<arrow::array::StringArray>()
            .ok_or_else(|| Error::InvalidInput("Parquet id column is not UTF-8".to_owned()))?;
        let structures = batch
            .column(1)
            .as_any()
            .downcast_ref::<arrow::array::StringArray>()
            .ok_or_else(|| {
                Error::InvalidInput("Parquet structure column is not UTF-8".to_owned())
            })?;
        let properties = batch
            .column(2)
            .as_any()
            .downcast_ref::<arrow::array::StringArray>()
            .ok_or_else(|| {
                Error::InvalidInput("Parquet properties column is not UTF-8".to_owned())
            })?;
        for row in 0..batch.num_rows() {
            let document = StructureDocument::from_json(structures.value(row))?;
            let properties: BTreeMap<String, serde_json::Value> =
                serde_json::from_str(properties.value(row))
                    .map_err(|error| Error::InvalidInput(format!("properties JSON: {error}")))?;
            let mut record = xtalforge_analysis::dataset::DatasetRecord::new(
                ids.value(row),
                document.structure,
            )?;
            for (key, value) in properties {
                record = record.with_property(key, crate::import::property_from_json(&value)?)?;
            }
            dataset.push(record)?;
        }
    }
    Ok(dataset)
}

/// Writes a dataset as a Parquet file.
///
/// The schema is three UTF-8 columns: `id`, `structure` (the versioned
/// structure document JSON) and `properties` (a JSON object). This keeps the
/// conversion lossless without prescribing an Arrow schema for every future
/// structure field.
///
/// # Errors
///
/// Returns [`Error::InvalidInput`] for I/O or Arrow/Parquet errors.
#[cfg(feature = "parquet")]
pub fn write_parquet(dataset: &Dataset, path: impl AsRef<std::path::Path>) -> Result<(), Error> {
    use std::sync::Arc;

    use arrow::array::StringArray;
    use arrow::datatypes::{DataType, Field, Schema};
    use arrow::record_batch::RecordBatch;
    use parquet::arrow::ArrowWriter;

    let mut ids = Vec::with_capacity(dataset.len());
    let mut structures = Vec::with_capacity(dataset.len());
    let mut properties = Vec::with_capacity(dataset.len());
    for record in dataset.iter() {
        ids.push(record.id.clone());
        structures.push(StructureDocument::new(record.structure.clone()).to_json()?);
        properties.push(
            serde_json::to_string(&record.properties)
                .map_err(|error| Error::InvalidInput(format!("properties JSON: {error}")))?,
        );
    }
    let schema = Arc::new(Schema::new(vec![
        Field::new("id", DataType::Utf8, false),
        Field::new("structure", DataType::Utf8, false),
        Field::new("properties", DataType::Utf8, false),
    ]));
    let batch = RecordBatch::try_new(
        schema.clone(),
        vec![
            Arc::new(StringArray::from(ids)),
            Arc::new(StringArray::from(structures)),
            Arc::new(StringArray::from(properties)),
        ],
    )
    .map_err(|error| Error::InvalidInput(format!("Parquet batch: {error}")))?;
    let file = std::fs::File::create(path.as_ref()).map_err(|error| {
        Error::InvalidInput(format!(
            "could not create {}: {error}",
            path.as_ref().display()
        ))
    })?;
    let mut writer = ArrowWriter::try_new(file, schema, None)
        .map_err(|error| Error::InvalidInput(format!("Parquet writer: {error}")))?;
    writer
        .write(&batch)
        .map_err(|error| Error::InvalidInput(format!("Parquet write: {error}")))?;
    writer
        .close()
        .map_err(|error| Error::InvalidInput(format!("Parquet close: {error}")))?;
    Ok(())
}

/// Reads a Parquet file written by [`write_parquet`].
///
/// # Errors
///
/// Returns [`Error::InvalidInput`] for I/O, schema or JSON errors.
#[cfg(feature = "parquet")]
pub fn read_parquet(path: impl AsRef<std::path::Path>) -> Result<Dataset, Error> {
    use parquet::arrow::arrow_reader::ParquetRecordBatchReaderBuilder;

    let file = std::fs::File::open(path.as_ref()).map_err(|error| {
        Error::InvalidInput(format!(
            "could not open {}: {error}",
            path.as_ref().display()
        ))
    })?;
    let builder = ParquetRecordBatchReaderBuilder::try_new(file)
        .map_err(|error| Error::InvalidInput(format!("Parquet reader: {error}")))?;
    let reader = builder
        .build()
        .map_err(|error| Error::InvalidInput(format!("Parquet reader: {error}")))?;
    let mut dataset = Dataset::new(
        path.as_ref()
            .file_stem()
            .and_then(std::ffi::OsStr::to_str)
            .unwrap_or("parquet"),
    );
    for batch in reader {
        let batch =
            batch.map_err(|error| Error::InvalidInput(format!("Parquet batch: {error}")))?;
        let ids = batch
            .column(0)
            .as_any()
            .downcast_ref::<arrow::array::StringArray>()
            .ok_or_else(|| Error::InvalidInput("Parquet id column is not UTF-8".to_owned()))?;
        let structures = batch
            .column(1)
            .as_any()
            .downcast_ref::<arrow::array::StringArray>()
            .ok_or_else(|| {
                Error::InvalidInput("Parquet structure column is not UTF-8".to_owned())
            })?;
        let properties = batch
            .column(2)
            .as_any()
            .downcast_ref::<arrow::array::StringArray>()
            .ok_or_else(|| {
                Error::InvalidInput("Parquet properties column is not UTF-8".to_owned())
            })?;
        for row in 0..batch.num_rows() {
            let document = StructureDocument::from_json(structures.value(row))?;
            let properties: BTreeMap<String, serde_json::Value> =
                serde_json::from_str(properties.value(row))
                    .map_err(|error| Error::InvalidInput(format!("properties JSON: {error}")))?;
            let mut record = xtalforge_analysis::dataset::DatasetRecord::new(
                ids.value(row),
                document.structure,
            )?;
            for (key, value) in properties {
                record = record.with_property(key, crate::import::property_from_json(&value)?)?;
            }
            dataset.push(record)?;
        }
    }
    Ok(dataset)
}
