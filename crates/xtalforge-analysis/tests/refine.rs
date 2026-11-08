//! Integration tests for whole-range and structural refinement
//! (Tasks 7.3.2, 7.8 and 7.9.1).

#![allow(
    clippy::unwrap_used,
    clippy::expect_used,
    clippy::panic,
    clippy::indexing_slicing
)]

use xtalforge_analysis::{
    PxrdOptions, RefineOptions, RefinedPhase, RietveldOptions, calculate_pxrd, combined_objective,
    refine_phase,
};
use xtalforge_core::{Lattice, Site, Structure, StructureBuilder};
use xtalforge_io::xy::PatternPoint;

fn silicon() -> Structure {
    StructureBuilder::new("Si", Lattice::cubic(5.430941).unwrap())
        .site(Site::new("Si1", "Si", [0.0, 0.0, 0.0]).unwrap())
        .build()
        .unwrap()
}

fn two_atom(position: [f64; 3]) -> Structure {
    StructureBuilder::new("Si2", Lattice::cubic(5.430941).unwrap())
        .site(Site::new("Si1", "Si", [0.0, 0.0, 0.0]).unwrap())
        .site(Site::new("Si2", "Si", position).unwrap())
        .build()
        .unwrap()
}

fn observed_from(structure: &Structure, scale: f64, background: f64) -> Vec<PatternPoint> {
    calculate_pxrd(structure, &PxrdOptions::default())
        .unwrap()
        .points
        .iter()
        .map(|point| PatternPoint {
            two_theta: point.two_theta,
            intensity: scale * point.intensity + background,
            uncertainty: None,
        })
        .collect()
}

fn parameter(result: &RefinedPhase, name: &str) -> f64 {
    result
        .parameters
        .iter()
        .find(|(key, _)| key == name)
        .map(|(_, value)| *value)
        .unwrap()
}

#[test]
fn synthetic_pattern_recovers_scale_and_background() {
    let structure = silicon();
    let observed = observed_from(&structure, 2.5, 7.0);
    let result = refine_phase(
        &structure,
        &observed,
        &RietveldOptions::default(),
        &RefineOptions::default(),
    )
    .unwrap();
    assert!(result.converged);
    assert_eq!(result.iterations, 0);
    assert!(result.rwp < 1e-6, "rwp = {}", result.rwp);
    assert!((parameter(&result, "scale") - 2.5).abs() < 1e-6);
    assert!((parameter(&result, "background_0") - 7.0).abs() < 1e-6);
    assert_eq!(result.structure, structure);
}

#[test]
fn refining_a_displaced_coordinate_lowers_rwp_or_converges() {
    let truth = two_atom([0.30, 0.20, 0.10]);
    let start = two_atom([0.33, 0.20, 0.10]);
    let observed = observed_from(&truth, 1.0, 0.0);
    let options = RietveldOptions::default();
    let initial = refine_phase(
        &start,
        &observed,
        &options,
        &RefineOptions {
            max_iterations: 0,
            refine_positions: true,
            ..RefineOptions::default()
        },
    )
    .unwrap();
    let refined = refine_phase(
        &start,
        &observed,
        &options,
        &RefineOptions {
            max_iterations: 10,
            tolerance: 1e-6,
            refine_positions: true,
            ..RefineOptions::default()
        },
    )
    .unwrap();
    assert!(
        refined.converged || refined.rwp < initial.rwp,
        "initial rwp = {}, refined rwp = {}",
        initial.rwp,
        refined.rwp
    );
    assert!(refined.rwp.is_finite());
    assert!(refined.parameters.iter().any(|(name, _)| name == "Si2.x"));
}

#[test]
fn too_many_free_parameters_are_rejected() {
    let structure = two_atom([0.30, 0.20, 0.10]);
    let observed: Vec<PatternPoint> = (0..5)
        .map(|index| PatternPoint {
            two_theta: 10.0 + index as f64 * 0.02,
            intensity: if index == 2 { 1.0 } else { 0.0 },
            uncertainty: None,
        })
        .collect();
    let error = refine_phase(
        &structure,
        &observed,
        &RietveldOptions::default(),
        &RefineOptions {
            refine_positions: true,
            refine_adps: true,
            refine_occupancies: true,
            ..RefineOptions::default()
        },
    )
    .unwrap_err();
    assert!(error.to_string().contains("structural parameters"));
    assert!(
        refine_phase(
            &structure,
            &[],
            &RietveldOptions::default(),
            &RefineOptions::default()
        )
        .is_err()
    );
    assert!(
        refine_phase(
            &structure,
            &observed,
            &RietveldOptions::default(),
            &RefineOptions {
                tolerance: 0.0,
                ..RefineOptions::default()
            }
        )
        .is_err()
    );
}

#[test]
fn combined_objective_is_monotonic_in_both_terms() {
    let base = combined_objective(1.0, 2.0, 0.5);
    assert!((base - 1.5).abs() < 1e-12);
    assert!(combined_objective(2.0, 2.0, 0.5) > base);
    assert!(combined_objective(1.0, 3.0, 0.5) > base);
    assert!((combined_objective(5.0, 2.0, 0.0) - 2.0).abs() < 1e-12);
    assert!((combined_objective(5.0, 2.0, 1.0) - 5.0).abs() < 1e-12);
    assert!(combined_objective(f64::NAN, 2.0, 0.5).is_infinite());

    let structure = silicon();
    let observed = observed_from(&structure, 2.0, 1.0);
    let result = refine_phase(
        &structure,
        &observed,
        &RietveldOptions::default(),
        &RefineOptions {
            combined_energy: Some(0.5),
            ..RefineOptions::default()
        },
    )
    .unwrap();
    let merit = parameter(&result, "combined_objective");
    assert!(merit.is_finite());
    assert!(merit >= 0.0);
}
