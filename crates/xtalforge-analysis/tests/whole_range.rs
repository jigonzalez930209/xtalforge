//! Integration tests for whole-range peak refinement, persistent parameter
//! sets, Pareto sweeps and cancellation (Task 14.1.5).

#![allow(
    clippy::unwrap_used,
    clippy::expect_used,
    clippy::panic,
    clippy::indexing_slicing
)]

use std::sync::Arc;
use std::sync::atomic::{AtomicUsize, Ordering};

use xtalforge_analysis::{
    Cancellation, InstrumentProfile, LineShape, ParameterSet, PeakOptions, RunControl,
    WholeRangeOptions, knee_index, pareto_sweep, refine_whole_range,
    refine_whole_range_with_control,
};
use xtalforge_io::xy::PatternPoint;

/// Synthetic lines: `(position, height)`.
const LINES: [(f64, f64); 3] = [(20.0, 10.0), (24.0, 8.0), (28.0, 6.0)];
/// Shared FWHM of the synthetic lines.
const FWHM: f64 = 0.25;

fn synthetic() -> Vec<PatternPoint> {
    let mut points = Vec::new();
    let mut x = 15.0;
    while x <= 33.0 {
        let mut intensity = 5.0 + 0.01 * (x - 15.0);
        for (position, height) in LINES {
            let delta = (x - position) / FWHM;
            intensity += height * (-4.0 * std::f64::consts::LN_2 * delta * delta).exp();
        }
        points.push(PatternPoint {
            two_theta: x,
            intensity,
            uncertainty: None,
        });
        x += 0.01;
    }
    points
}

fn options() -> WholeRangeOptions {
    WholeRangeOptions {
        detection: PeakOptions {
            min_height: 6.0,
            min_distance: 1.0,
            smoothing: 1,
            background: false,
            refine: false,
            ..PeakOptions::default()
        },
        shape: LineShape::Gaussian,
        eta: 0.0,
        refine_eta: false,
        max_peaks: 10,
        max_iterations: 100,
        tolerance: 1e-10,
    }
}

#[test]
fn whole_range_fit_recovers_three_overlapping_lines() {
    let points = synthetic();
    let fit = refine_whole_range(&points, &options()).unwrap();
    assert_eq!(fit.peaks.len(), LINES.len(), "peaks {:?}", fit.peaks);
    for (index, (position, height)) in LINES.iter().enumerate() {
        let peak = &fit.peaks[index];
        assert!(
            (peak.two_theta - position).abs() < 0.02,
            "line {index} at {} != {position}",
            peak.two_theta
        );
        assert!(
            (peak.intensity - height).abs() < 0.5,
            "line {index} height {} != {height}",
            peak.intensity
        );
        assert!((peak.fwhm - FWHM).abs() < 0.03);
    }
    assert!((fit.fwhm - FWHM).abs() < 0.02, "fwhm {}", fit.fwhm);
    assert!(
        (fit.background[0] - 5.0).abs() < 0.2,
        "b0 {:?}",
        fit.background
    );
    assert!(fit.background[1].abs() < 0.2, "b1 {:?}", fit.background);
    assert!(fit.rwp < 0.01, "Rwp {}", fit.rwp);
    assert!(fit.chi_squared.is_finite());
}

#[test]
fn cancellation_and_progress_are_cooperative() {
    let points = synthetic();
    let token = Cancellation::new();
    token.cancel();
    let control = RunControl::new().with_cancellation(token.clone());
    let error = refine_whole_range_with_control(&points, &options(), &control).unwrap_err();
    assert_eq!(error.code(), "invalid_input");
    assert!(error.message().contains("cancelled"), "{}", error.message());
    assert!(token.is_cancelled());
    token.reset();
    assert!(!token.is_cancelled());

    let events = Arc::new(AtomicUsize::new(0));
    let counter = Arc::clone(&events);
    let control = RunControl::new().with_progress(move |progress| {
        assert_eq!(progress.stage, "whole_range");
        assert!(progress.value.is_some());
        counter.fetch_add(1, Ordering::SeqCst);
    });
    let fit = refine_whole_range_with_control(&points, &options(), &control).unwrap();
    assert!(fit.rwp < 0.01);
    assert!(events.load(Ordering::SeqCst) > 0);
}

#[test]
fn whole_range_option_errors_are_typed() {
    let points = synthetic();
    assert!(refine_whole_range(&points[..4], &options()).is_err());
    assert!(
        refine_whole_range(
            &points,
            &WholeRangeOptions {
                shape: LineShape::Tch,
                ..options()
            }
        )
        .is_err()
    );
    assert!(
        refine_whole_range(
            &points,
            &WholeRangeOptions {
                max_peaks: 0,
                ..options()
            }
        )
        .is_err()
    );
    let quiet: Vec<PatternPoint> = points
        .iter()
        .map(|point| PatternPoint {
            intensity: 1.0,
            ..*point
        })
        .collect();
    assert!(refine_whole_range(&quiet, &options()).is_err());
}

#[test]
fn parameter_sets_round_trip_json_and_files() {
    let mut set = ParameterSet::new("rutile");
    set.set("eta", 0.62).unwrap();
    set.set("profile.fwhm", 0.11).unwrap();
    set.set("caglioti.u", 0.004).unwrap();
    set.set("caglioti.v", -0.002).unwrap();
    set.set("caglioti.w", 0.012).unwrap();
    set.set("asymmetry.fcj1", 0.02).unwrap();
    set.set("asymmetry.fcj2", 0.03).unwrap();
    set.set("monochromator.polarization", 0.8).unwrap();
    assert!(set.set("bad", f64::NAN).is_err());
    assert!(set.set("", 1.0).is_err());
    assert_eq!(set.get("eta"), Some(0.62));

    let text = set.to_json().unwrap();
    let parsed = ParameterSet::from_json(&text).unwrap();
    assert_eq!(parsed, set);
    assert_eq!(parsed.version, 1);
    assert!(ParameterSet::from_json("not json").is_err());

    let mut other = set.clone();
    other.version = 2;
    assert!(set.merge(&other).is_err());
    let mut merged = ParameterSet::new("merge");
    merged.merge(&set).unwrap();
    assert_eq!(merged.get("eta"), Some(0.62));

    let profile = set.instrument_profile(InstrumentProfile::default());
    assert_eq!(profile.caglioti.unwrap().u, 0.004);
    assert_eq!(profile.asymmetry.unwrap().fcj1, 0.02);
    assert_eq!(profile.monochromator.polarization, 0.8);
    let collected = ParameterSet::from_instrument_profile("instrument", &profile);
    assert_eq!(collected.get("caglioti.w"), Some(0.012));

    let path = std::env::temp_dir().join(format!(
        "xtalforge-parameters-{}-{}.json",
        std::process::id(),
        "whole-range"
    ));
    set.save(&path).unwrap();
    let loaded = ParameterSet::load(&path).unwrap();
    assert_eq!(loaded, set);
    std::fs::remove_file(&path).unwrap();
}

#[test]
fn pareto_sweep_finds_the_front_and_the_knee() {
    let weights = [0.0, 0.25, 0.5, 0.75, 1.0];
    let sweep = pareto_sweep(&weights, |weight| Ok((weight, (1.0 - weight).powi(2)))).unwrap();
    assert_eq!(sweep.samples.len(), weights.len());
    assert_eq!(sweep.front.len(), weights.len());
    assert_eq!(sweep.knee, Some(2));
    assert!((sweep.samples[2].combined - (0.5 * 0.5 + 0.5 * 0.25)).abs() < 1e-12);

    let straight: Vec<(f64, f64)> = (0..5)
        .map(|index| (index as f64, 4.0 - index as f64))
        .collect();
    assert_eq!(knee_index(&straight), None);
    assert_eq!(knee_index(&[]), None);
    assert_eq!(knee_index(&[(0.0, 1.0), (1.0, 0.0)]), None);

    let knee_points = [(0.0, 1.0), (0.5, 0.25), (1.0, 0.0)];
    assert_eq!(knee_index(&knee_points), Some(1));

    assert!(pareto_sweep(&[], |weight| Ok((weight, weight))).is_err());
    assert!(pareto_sweep(&[0.5], |_| Ok((f64::NAN, 0.0))).is_err());
    assert!(pareto_sweep(&[1.5], |weight| Ok((weight, weight))).is_err());
}
