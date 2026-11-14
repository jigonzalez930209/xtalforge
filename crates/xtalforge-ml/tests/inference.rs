//! Property, ensemble, selection and active-learning tests (Tasks 21.3-21.4).

#![allow(
    clippy::unwrap_used,
    clippy::expect_used,
    clippy::panic,
    clippy::indexing_slicing
)]

use xtalforge_analysis::dataset::{Dataset, DatasetRecord, PropertyValue};
use xtalforge_core::{Lattice, Site, StructureBuilder};
use xtalforge_ml::{
    ActiveLearningOptions, Ensemble, Predictor, PropertyModel, SelectionOptions, SelectionStrategy,
    active_learning_loop, active_learning_workflow, append_predictions, predict_dataset,
    select_batch, uncertainty,
};

fn structure(index: usize) -> xtalforge_core::Structure {
    let lattice = Lattice::from_parameters(5.0, 5.0, 5.0, 90.0, 90.0, 90.0).expect("cell");
    StructureBuilder::new(format!("s{index}"), lattice)
        .site(Site::new("C1", "C", [0.1 * index as f64, 0.0, 0.0]).expect("C"))
        .build()
        .expect("structure")
}

fn dataset() -> Dataset {
    let mut dataset = Dataset::new("inference");
    for index in 0..4 {
        let mut record = DatasetRecord::new(format!("r{index}"), structure(index)).expect("record");
        record = record
            .with_property("target", PropertyValue::Number(index as f64))
            .expect("property");
        dataset.push(record).expect("push");
    }
    dataset
}

#[test]
fn linear_model_predicts_and_reports_errors() {
    let model = PropertyModel::new("linear", vec![2.0, -1.0], 0.5).expect("model");
    assert_eq!(model.predict(&[1.0, 2.0]).expect("predict"), 0.5);
    assert!(model.predict(&[1.0]).is_err());
    let normalized = model
        .clone()
        .with_normalization(vec![1.0, 1.0], vec![2.0, 2.0])
        .expect("normalized");
    assert!(normalized.predict(&[1.0, 1.0]).is_ok());
    let text = model.to_json().expect("json");
    assert_eq!(PropertyModel::from_json(&text).expect("parse"), model);
}

#[test]
fn dataset_prediction_and_append() {
    let model = PropertyModel::new("m", vec![1.0], 0.0).expect("model");
    let prediction = predict_dataset(&dataset(), &model, |structure| {
        Ok(vec![structure.lattice().volume() / 1000.0])
    })
    .expect("predict");
    assert_eq!(prediction.predictions.len(), 4);
    let appended = append_predictions(&dataset(), &prediction, "predicted").expect("append");
    assert!(
        appended
            .iter()
            .all(|record| record.properties.contains_key("predicted"))
    );
}

#[test]
fn ensemble_reports_disagreement() {
    let first = PropertyModel::new("a", vec![1.0], 0.0).expect("a");
    let second = PropertyModel::new("b", vec![1.0], 1.0).expect("b");
    let ensemble = Ensemble::new(vec![first, second]).expect("ensemble");
    let prediction = ensemble.predict(&[2.0]).expect("predict");
    assert!((prediction.mean - 2.5).abs() < 1e-12);
    assert!((prediction.std - 0.5).abs() < 1e-12);
    assert_eq!(uncertainty(&[1.0, 1.0]).std, 0.0);
    assert!(Ensemble::new(Vec::new()).is_err());
}

#[test]
fn selection_strategies_return_unique_indices() {
    let features = vec![
        vec![0.0, 0.0],
        vec![1.0, 0.0],
        vec![5.0, 5.0],
        vec![5.1, 5.0],
    ];
    let uncertainties = vec![0.1, 0.9, 0.2, 0.8];
    for strategy in [
        SelectionStrategy::Uncertainty,
        SelectionStrategy::FarthestPoint,
        SelectionStrategy::Diversity,
    ] {
        let selected = select_batch(
            &features,
            &uncertainties,
            &SelectionOptions {
                batch: 2,
                strategy,
                diversity_radius: 0.5,
            },
        )
        .expect("select");
        assert_eq!(selected.len(), 2);
        assert_ne!(selected[0], selected[1]);
    }
    assert!(
        select_batch(
            &features,
            &uncertainties,
            &SelectionOptions {
                batch: 0,
                ..SelectionOptions::default()
            }
        )
        .is_err()
    );
}

#[test]
fn active_learning_runs_rounds_and_appends_labels() {
    let model = Box::new(PropertyModel::new("m", vec![1.0], 0.0).expect("model"));
    let report = active_learning_loop(
        &dataset(),
        model,
        |structure| Ok(vec![structure.sites()[0].frac()[0]]),
        |model, features| Ok(model.predict(features).unwrap().abs()),
        |structure| Ok(structure.sites()[0].frac()[0] * 10.0),
        |_dataset| Ok(None),
        &ActiveLearningOptions {
            rounds: 2,
            selection: SelectionOptions {
                batch: 2,
                strategy: SelectionStrategy::Uncertainty,
                diversity_radius: 0.0,
            },
            target: "computed".to_owned(),
            stop_when_empty: true,
        },
    )
    .expect("active learning");
    assert_eq!(report.rounds.len(), 2);
    assert_eq!(report.dataset.len(), 8);
    assert!(
        report
            .dataset
            .iter()
            .any(|record| record.properties.contains_key("computed"))
    );
}

#[test]
fn workflow_document_has_all_steps() {
    let document = active_learning_workflow("al", "data.json", "model.json", 4, 2, None);
    assert_eq!(document["workflow_version"], 1);
    let steps = document["steps"].as_array().expect("steps");
    assert_eq!(steps.len(), 6);
    assert_eq!(steps[1]["verb"], "search");
    assert_eq!(steps[2]["verb"], "ml_predict");
    assert_eq!(steps[3]["verb"], "ml_select");
}
