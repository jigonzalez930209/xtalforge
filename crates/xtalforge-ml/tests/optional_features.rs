//! Parquet and ONNX feature tests (Tasks 21.2.1 and 21.3.1).

#![cfg(any(feature = "parquet", feature = "onnx"))]
#![allow(
    clippy::unwrap_used,
    clippy::expect_used,
    clippy::panic,
    clippy::indexing_slicing
)]

use xtalforge_analysis::dataset::{Dataset, DatasetRecord, PropertyValue};
use xtalforge_core::{Lattice, Site, StructureBuilder};

fn structure() -> xtalforge_core::Structure {
    let lattice = Lattice::from_parameters(5.0, 5.0, 5.0, 90.0, 90.0, 90.0).expect("cell");
    StructureBuilder::new("C", lattice)
        .site(Site::new("C1", "C", [0.0, 0.0, 0.0]).expect("C"))
        .build()
        .expect("structure")
}

#[cfg(feature = "parquet")]
#[test]
fn parquet_round_trip_is_lossless() {
    use xtalforge_ml::{read_parquet, write_parquet};

    let mut dataset = Dataset::new("parquet");
    let mut first = DatasetRecord::new("a", structure()).expect("record");
    first = first
        .with_property("energy", PropertyValue::Number(-1.25))
        .expect("property");
    dataset.push(first).expect("push");
    dataset
        .push(DatasetRecord::new("b", structure()).expect("record"))
        .expect("push");
    let path = std::env::temp_dir().join(format!("xtalforge-ml-{}.parquet", std::process::id()));
    write_parquet(&dataset, &path).expect("write");
    let read = read_parquet(&path).expect("read");
    assert_eq!(read.len(), dataset.len());
    let first = read.get(0).expect("first");
    assert_eq!(first.id, "a");
    assert_eq!(
        first
            .properties
            .get("energy")
            .and_then(PropertyValue::as_f64),
        Some(-1.25)
    );
    let _ = std::fs::remove_file(path);
}

#[cfg(feature = "parquet")]
#[test]
fn parquet_missing_file_errors() {
    use xtalforge_ml::read_parquet;
    assert!(read_parquet("/definitely/missing/file.parquet").is_err());
}

#[cfg(feature = "onnx")]
#[test]
fn invalid_onnx_model_is_rejected() {
    use xtalforge_ml::{OnnxModel, OnnxOptions};
    assert!(
        OnnxModel::from_bytes(
            b"not an onnx model",
            &OnnxOptions {
                n_features: 4,
                input: None
            }
        )
        .is_err()
    );
    assert!(
        OnnxModel::from_bytes(
            b"",
            &OnnxOptions {
                n_features: 0,
                input: None
            }
        )
        .is_err()
    );
}
