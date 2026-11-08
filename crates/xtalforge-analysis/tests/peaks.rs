//! Integration tests for measured-pattern peak finding (Tasks 7.2/7.3).

#![allow(
    clippy::unwrap_used,
    clippy::expect_used,
    clippy::panic,
    clippy::indexing_slicing
)]

use xtalforge_analysis::{PeakOptions, find_peaks};
use xtalforge_io::xy::PatternPoint;

fn gaussian(two_theta: f64, center: f64, height: f64, width: f64) -> f64 {
    height * (-((two_theta - center).powi(2)) / (2.0 * width * width)).exp()
}

#[test]
fn finds_two_gaussian_peaks() {
    let points: Vec<PatternPoint> = (0..600)
        .map(|index| {
            let two_theta = 5.0 + index as f64 * 0.02;
            PatternPoint {
                two_theta,
                intensity: gaussian(two_theta, 10.0, 10.0, 0.2)
                    + gaussian(two_theta, 13.0, 5.0, 0.3),
                uncertainty: None,
            }
        })
        .collect();
    let peaks = find_peaks(
        &points,
        &PeakOptions {
            min_height: 1.0,
            min_distance: 0.5,
            smoothing: 1,
        },
    )
    .unwrap();
    assert_eq!(peaks.len(), 2);
    assert!((peaks[0].two_theta - 10.0).abs() < 0.05);
    assert!((peaks[1].two_theta - 13.0).abs() < 0.05);
    assert!(peaks[0].intensity > peaks[1].intensity);
    assert!(peaks[0].fwhm > 0.0);
}

#[test]
fn bad_inputs_are_rejected() {
    let points = vec![
        PatternPoint {
            two_theta: 1.0,
            intensity: 1.0,
            uncertainty: None,
        },
        PatternPoint {
            two_theta: 2.0,
            intensity: 2.0,
            uncertainty: None,
        },
        PatternPoint {
            two_theta: 3.0,
            intensity: 1.0,
            uncertainty: None,
        },
    ];
    assert!(find_peaks(&points, &PeakOptions::default()).is_ok());
    assert!(find_peaks(&points[..2], &PeakOptions::default()).is_err());
    assert!(
        find_peaks(
            &points,
            &PeakOptions {
                min_distance: 0.0,
                ..Default::default()
            }
        )
        .is_err()
    );
}
