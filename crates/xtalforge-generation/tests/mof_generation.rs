//! MOF dataset, orientation and functionalisation tests (Task 20.4).

#![allow(
    clippy::unwrap_used,
    clippy::expect_used,
    clippy::panic,
    clippy::indexing_slicing
)]

use xtalforge_generation::{
    MofDatasetOptions, PrototypeOptions, build_prototype, functionalization_sweep,
    generate_mof_dataset, global_orientation_fit,
};
use xtalforge_mof::FunctionalGroup;

#[test]
fn orientation_fit_recovers_a_known_rotation() {
    let connections = [
        [1.0, 0.0, 0.0],
        [-1.0, 0.0, 0.0],
        [0.0, 1.0, 0.0],
        [0.0, -1.0, 0.0],
    ];
    let targets = [
        [0.0, 1.0, 0.0],
        [0.0, -1.0, 0.0],
        [-1.0, 0.0, 0.0],
        [1.0, 0.0, 0.0],
    ];
    let fit = global_orientation_fit(&connections, &targets).expect("fit");
    assert!(fit.rmsd < 1e-9, "rmsd {}", fit.rmsd);
    assert_eq!(fit.assignment.len(), 4);
}

#[test]
fn orientation_fit_rejects_bad_input() {
    assert!(global_orientation_fit(&[], &[]).is_err());
    assert!(global_orientation_fit(&[[1.0, 0.0, 0.0]], &[]).is_err());
    assert!(global_orientation_fit(&[[0.0, 0.0, 0.0]], &[[1.0, 0.0, 0.0]]).is_err());
    let many = vec![[1.0, 0.0, 0.0]; 7];
    assert!(global_orientation_fit(&many, &many).is_err());
}

#[test]
fn mof_dataset_builds_curated_combinations() {
    let options = MofDatasetOptions {
        limit: 2,
        descriptors: false,
        relax: false,
        seed: 1,
        ..MofDatasetOptions::default()
    };
    let report = generate_mof_dataset(&options).expect("dataset");
    assert_eq!(report.records.len(), 2);
    assert_eq!(report.combinations, 2);
    assert_eq!(report.built + report.failed, 2);
}

#[test]
fn mof_dataset_can_select_names() {
    let options = MofDatasetOptions {
        topologies: vec!["pcu".to_owned()],
        limit: 1,
        descriptors: false,
        relax: false,
        seed: 0,
        ..MofDatasetOptions::default()
    };
    let report = generate_mof_dataset(&options).expect("dataset");
    assert_eq!(report.records.len(), 1);
    assert_eq!(report.records[0].topology, "pcu");
}

#[test]
fn mof_dataset_validates_options() {
    assert!(
        generate_mof_dataset(&MofDatasetOptions {
            limit: 0,
            ..MofDatasetOptions::default()
        })
        .is_err()
    );
    assert!(
        generate_mof_dataset(&MofDatasetOptions {
            pore_spacing: 0.0,
            ..MofDatasetOptions::default()
        })
        .is_err()
    );
}

#[test]
fn functionalization_sweep_covers_group_and_fraction() {
    let structure = build_prototype("A1", &PrototypeOptions::default()).expect("A1");
    let sweep = functionalization_sweep(
        &structure,
        &[FunctionalGroup::Amino, FunctionalGroup::Fluorine],
        &[0.0, 1.0],
        7,
        false,
    )
    .expect("sweep");
    assert_eq!(sweep.records.len(), 4);
    assert_eq!(sweep.built + sweep.failed, 4);
}

#[test]
fn functionalization_sweep_validates_options() {
    let structure = build_prototype("A1", &PrototypeOptions::default()).expect("A1");
    assert!(functionalization_sweep(&structure, &[], &[1.0], 1, false).is_err());
    assert!(functionalization_sweep(&structure, &[FunctionalGroup::Amino], &[], 1, false).is_err());
    assert!(
        functionalization_sweep(&structure, &[FunctionalGroup::Amino], &[2.0], 1, false).is_err()
    );
}
