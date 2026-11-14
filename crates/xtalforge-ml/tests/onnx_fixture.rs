//! ONNX reference test against the exported PyTorch linear model (Task 21.3.1).

#![cfg(feature = "onnx")]
#![allow(
    clippy::unwrap_used,
    clippy::expect_used,
    clippy::panic,
    clippy::indexing_slicing
)]

use xtalforge_core::{Lattice, Site, StructureBuilder};
use xtalforge_ff::Calculator;
use xtalforge_ml::{OnnxCalculator, OnnxDescriptor, OnnxModel, OnnxOptions};

const MODEL: &[u8] = include_bytes!("../data/onnx_linear.onnx");
const REFERENCE: &str = include_str!("../data/onnx_linear.json");

fn structure() -> xtalforge_core::Structure {
    let lattice =
        Lattice::new([[20.0, 0.0, 0.0], [0.0, 20.0, 0.0], [0.0, 0.0, 20.0]]).expect("cell");
    StructureBuilder::new("onnx", lattice)
        .site(Site::new("C1", "C", [0.0, 0.0, 0.0]).expect("C"))
        .site(Site::new("C2", "C", [0.1, 0.0, 0.0]).expect("C"))
        .build()
        .expect("structure")
}

#[test]
fn onnx_energy_and_forces_match_pytorch() {
    let reference: serde_json::Value = serde_json::from_str(REFERENCE).expect("reference");
    let model = OnnxModel::from_bytes(
        MODEL,
        &OnnxOptions {
            n_features: 6,
            input: Some("features".to_owned()),
        },
    )
    .expect("model");
    let calculator =
        OnnxCalculator::new(model, OnnxDescriptor::Cartesian, 1e-3, 1e-3).expect("calculator");
    let structure = structure();
    let energy = Calculator::energy(&calculator, &structure).expect("energy");
    let expected_energy = reference["energy"].as_f64().expect("energy");
    assert!(
        (energy - expected_energy).abs() < 1e-5,
        "energy {energy} differs from PyTorch {expected_energy}"
    );
    let (_, forces) = Calculator::energy_forces(&calculator, &structure).expect("forces");
    let expected = reference["forces"].as_array().expect("forces");
    for (atom, expected) in forces.iter().zip(expected) {
        let row = expected.as_array().expect("row");
        for (axis, value) in atom.iter().enumerate() {
            let reference = row[axis].as_f64().expect("force");
            assert!(
                (value - reference).abs() < 2e-4,
                "force {value} differs from PyTorch {reference}"
            );
        }
    }
}

#[test]
fn onnx_rejects_wrong_width() {
    let model = OnnxModel::from_bytes(
        MODEL,
        &OnnxOptions {
            n_features: 6,
            input: None,
        },
    )
    .expect("model");
    assert!(model.predict(&[0.0; 5]).is_err());
}
