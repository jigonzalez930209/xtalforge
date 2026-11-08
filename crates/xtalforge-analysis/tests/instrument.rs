//! Integration tests for instrument corrections and the K alpha doublet
//! (Task 7.2.2).

#![allow(
    clippy::unwrap_used,
    clippy::expect_used,
    clippy::panic,
    clippy::indexing_slicing
)]

use xtalforge_analysis::{Geometry, Instrument, doublet, doublet_checked};

#[test]
fn zero_shift_moves_the_reported_peak() {
    let instrument = Instrument {
        zero_shift: 0.1,
        ..Instrument::default()
    };
    let corrected = instrument.correct(30.0);
    assert!((corrected - 30.1).abs() < 1e-12);
    assert!((Instrument::default().correct(30.0) - 30.0).abs() < 1e-12);
}

#[test]
fn displacement_has_the_documented_sign() {
    let radius = 240.0;
    let displacement = 1.0;
    let instrument = Instrument {
        displacement,
        sample_radius: radius,
        ..Instrument::default()
    };
    let two_theta: f64 = 30.0;
    let theta = (0.5 * two_theta).to_radians();
    let expected = two_theta - (2.0 * displacement * theta.cos() / radius).to_degrees();
    let corrected = instrument.correct(two_theta);
    assert!(corrected < two_theta);
    assert!((corrected - expected).abs() < 1e-12);

    let capillary = Instrument {
        geometry: Geometry::DebyeScherrer,
        ..instrument
    };
    assert!((capillary.correct(two_theta) - two_theta).abs() < 1e-12);
}

#[test]
fn doublet_weights_sum_to_one_and_ka2_follows_bragg() {
    let ka1 = 1.540_598;
    let ka2 = 1.544_39;
    let pair = doublet(28.44, ka1, ka2, 0.5);
    assert_eq!(pair.len(), 2);
    assert!((pair[0].1 + pair[1].1 - 1.0).abs() < 1e-12);
    assert!((pair[0].1 - 2.0 / 3.0).abs() < 1e-12);
    assert!((pair[1].1 - 1.0 / 3.0).abs() < 1e-12);
    assert!(pair[1].0 > pair[0].0);
    let expected = 2.0
        * ((0.5 * 28.44_f64.to_radians()).sin() * ka2 / ka1)
            .asin()
            .to_degrees();
    assert!((pair[1].0 - expected).abs() < 1e-12);

    let only_ka1 = doublet(28.44, ka1, ka2, 0.0);
    assert!((only_ka1[0].1 - 1.0).abs() < 1e-12);
    assert!(only_ka1[1].1.abs() < 1e-12);
}

#[test]
fn invalid_values_error() {
    let instrument = Instrument::default();
    assert!(instrument.validate().is_ok());
    assert!(instrument.doublet(28.44, 1.540_598, 1.544_39).is_ok());
    assert!(doublet_checked(0.0, 1.540_598, 1.544_39, 0.5).is_err());
    assert!(doublet_checked(28.44, -1.0, 1.544_39, 0.5).is_err());
    assert!(doublet_checked(28.44, 1.540_598, 1.544_39, 2.0).is_err());
    assert!(
        Instrument {
            sample_radius: 0.0,
            ..Instrument::default()
        }
        .validate()
        .is_err()
    );
    assert!(
        Instrument {
            ka2_ratio: 1.5,
            ..Instrument::default()
        }
        .validate()
        .is_err()
    );
    assert!(
        Instrument {
            polarization: -0.1,
            ..Instrument::default()
        }
        .validate()
        .is_err()
    );
    assert!(
        Instrument {
            zero_shift: f64::NAN,
            ..Instrument::default()
        }
        .validate()
        .is_err()
    );
}
