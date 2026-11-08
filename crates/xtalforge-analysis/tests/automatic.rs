//! Integration tests for the automatic peaks-index-Le Bail-Rietveld run
//! (Task 7.10).

#![allow(
    clippy::unwrap_used,
    clippy::expect_used,
    clippy::panic,
    clippy::indexing_slicing
)]

use xtalforge_analysis::{AutomaticOptions, PxrdOptions, Radiation, calculate_pxrd, run_automatic};
use xtalforge_core::{Lattice, Site, Structure, StructureBuilder};
use xtalforge_io::xy::PatternPoint;

fn silicon() -> Structure {
    StructureBuilder::new("Si", Lattice::cubic(5.430941).unwrap())
        .site(Site::new("Si1", "Si", [0.0, 0.0, 0.0]).unwrap())
        .build()
        .unwrap()
}

fn pattern() -> Vec<PatternPoint> {
    calculate_pxrd(
        &silicon(),
        &PxrdOptions {
            radiation: Radiation::CuKa,
            two_theta_range: [10.0, 80.0],
            step: 0.01,
            fwhm: 0.05,
        },
    )
    .unwrap()
    .points
    .iter()
    .map(|point| PatternPoint {
        two_theta: point.two_theta,
        intensity: point.intensity,
        uncertainty: None,
    })
    .collect()
}

#[test]
fn automatic_run_records_steps_in_order() {
    let observed = pattern();
    let result = run_automatic(&silicon(), &observed, &AutomaticOptions::default()).unwrap();
    assert_eq!(result.steps, vec!["peaks", "index", "le_bail", "rietveld"]);
    assert!(!result.peaks.is_empty());
    assert!(result.lebail.is_some());
    assert!(result.phase.is_some());
    let phase = result.phase.unwrap();
    assert!(phase.rwp.is_finite());
    assert!(phase.rwp >= 0.0);
    if let Some(index) = result.index {
        assert!(index.a.is_finite());
        assert!(index.completeness >= AutomaticOptions::default().min_completeness);
    }
}

#[test]
fn empty_and_peakless_inputs_error() {
    let structure = silicon();
    assert!(run_automatic(&structure, &[], &AutomaticOptions::default()).is_err());
    let flat: Vec<PatternPoint> = (0..100)
        .map(|index| PatternPoint {
            two_theta: 10.0 + index as f64 * 0.02,
            intensity: 0.0,
            uncertainty: None,
        })
        .collect();
    assert!(run_automatic(&structure, &flat, &AutomaticOptions::default()).is_err());
    assert!(
        run_automatic(
            &structure,
            &flat,
            &AutomaticOptions {
                min_completeness: 2.0,
                ..AutomaticOptions::default()
            }
        )
        .is_err()
    );
}
