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
            ..PeakOptions::default()
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

#[test]
fn background_envelope_follows_the_baseline_not_the_peaks() {
    let points: Vec<PatternPoint> = (0..600)
        .map(|index| {
            let two_theta = 10.0 + index as f64 * 0.02;
            let baseline = 50.0 + 1.0 * two_theta;
            PatternPoint {
                two_theta,
                intensity: baseline + gaussian(two_theta, 16.0, 30.0, 0.08),
                uncertainty: None,
            }
        })
        .collect();
    let envelope = xtalforge_analysis::background_envelope(&points, 40);
    assert_eq!(envelope.len(), points.len());
    let peak = 300;
    assert!(envelope[peak] < 50.0 + points[peak].two_theta + 2.0);
    for index in [50, 150, 450, 550] {
        let baseline = 50.0 + points[index].two_theta;
        assert!(
            (envelope[index] - baseline).abs() < 0.05 * baseline,
            "{index}: {} != {baseline}",
            envelope[index]
        );
    }
}

#[test]
fn fit_peaks_recovers_center_height_and_width() {
    let points: Vec<PatternPoint> = (0..400)
        .map(|index| {
            let two_theta = 18.0 + index as f64 * 0.01;
            PatternPoint {
                two_theta,
                intensity: 2.0 + 0.5 * (two_theta - 18.0) + gaussian(two_theta, 20.0, 10.0, 0.08),
                uncertainty: None,
            }
        })
        .collect();
    let options = PeakOptions {
        fit_window: 0.5,
        refine: true,
        ..PeakOptions::default()
    };
    let peaks = xtalforge_analysis::fit_peaks(&points, &[20.0], &options).unwrap();
    assert_eq!(peaks.len(), 1);
    assert!((peaks[0].two_theta - 20.0).abs() < 2e-3);
    assert!((peaks[0].intensity - 10.0).abs() < 0.2);
    assert!((peaks[0].fwhm - 0.188).abs() < 0.02);
}

#[test]
fn background_and_refine_find_a_weak_peak() {
    let points: Vec<PatternPoint> = (0..900)
        .map(|index| {
            let two_theta = 10.0 + index as f64 * 0.02;
            let baseline = 40.0 + 15.0 * ((two_theta - 10.0) / 6.0).sin();
            PatternPoint {
                two_theta,
                intensity: baseline + gaussian(two_theta, 22.0, 1.5, 0.12),
                uncertainty: None,
            }
        })
        .collect();
    let peaks = find_peaks(
        &points,
        &PeakOptions {
            min_height: 0.5,
            min_distance: 0.3,
            background: true,
            background_window: 60,
            refine: true,
            ..PeakOptions::default()
        },
    )
    .unwrap();
    let found = peaks
        .iter()
        .find(|peak| (peak.two_theta - 22.0).abs() < 0.2)
        .expect("weak peak");
    assert!((found.two_theta - 22.0).abs() < 0.02);
}

#[test]
fn fit_peaks_validates_its_input() {
    let points: Vec<PatternPoint> = (0..10)
        .map(|index| PatternPoint {
            two_theta: index as f64 * 0.1,
            intensity: 1.0,
            uncertainty: None,
        })
        .collect();
    assert!(xtalforge_analysis::fit_peaks(&points, &[], &PeakOptions::default()).is_err());
    assert!(xtalforge_analysis::fit_peaks(&points[..4], &[0.5], &PeakOptions::default()).is_err());
    assert!(
        xtalforge_analysis::fit_peaks(
            &points,
            &[0.5],
            &PeakOptions {
                fit_window: 0.0,
                ..PeakOptions::default()
            }
        )
        .is_err()
    );
}
