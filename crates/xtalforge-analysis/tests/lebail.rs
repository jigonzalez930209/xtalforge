//! Integration tests for Le Bail extraction (Task 7.6.1).

#![allow(
    clippy::unwrap_used,
    clippy::expect_used,
    clippy::panic,
    clippy::indexing_slicing
)]

use xtalforge_analysis::{LebailOptions, le_bail};
use xtalforge_core::{Lattice, Site, StructureBuilder};
use xtalforge_io::xy::PatternPoint;

fn silicon() -> xtalforge_core::Structure {
    StructureBuilder::new("Si", Lattice::cubic(5.430941).unwrap())
        .site(Site::new("Si1", "Si", [0.0, 0.0, 0.0]).unwrap())
        .build()
        .unwrap()
}

#[test]
fn le_bail_extracts_positive_intensities_and_fits() {
    let structure = silicon();
    let points: Vec<PatternPoint> = (0..2000)
        .map(|index| {
            let two_theta = 5.0 + index as f64 * 0.02;
            let intensity = (-((two_theta - 28.44).powi(2)) / (2.0 * 0.08 * 0.08)).exp()
                + 0.5 * (-((two_theta - 47.3).powi(2)) / (2.0 * 0.08 * 0.08)).exp();
            PatternPoint {
                two_theta,
                intensity,
                uncertainty: None,
            }
        })
        .collect();
    let result = le_bail(&structure, &points, &LebailOptions::default()).unwrap();
    assert!(!result.intensities.is_empty());
    assert!(
        result
            .intensities
            .iter()
            .all(|value| value.is_finite() && *value >= 0.0)
    );
    assert!(result.r_factor >= 0.0);
    assert_eq!(result.pattern.len(), points.len());
}

#[test]
fn invalid_inputs_are_rejected() {
    let structure = silicon();
    let points = vec![PatternPoint {
        two_theta: 10.0,
        intensity: 1.0,
        uncertainty: None,
    }];
    assert!(le_bail(&structure, &points, &LebailOptions::default()).is_err());
    assert!(
        le_bail(
            &structure,
            &points,
            &LebailOptions {
                fwhm: 0.0,
                ..Default::default()
            }
        )
        .is_err()
    );
}

#[test]
fn pawley_recovers_positive_intensities() {
    use xtalforge_analysis::pawley;
    let structure = silicon();
    let points: Vec<PatternPoint> = (0..2000)
        .map(|index| {
            let two_theta = 5.0 + index as f64 * 0.02;
            let intensity = (-((two_theta - 28.44).powi(2)) / (2.0 * 0.08 * 0.08)).exp();
            PatternPoint {
                two_theta,
                intensity,
                uncertainty: None,
            }
        })
        .collect();
    let result = pawley(&structure, &points, &LebailOptions::default()).unwrap();
    assert!(result.intensities.iter().all(|value| value.is_finite()));
    assert!(result.intensities.iter().any(|value| *value > 0.0));
    assert!(result.r_factor >= 0.0);
}
