//! Integration tests for profiles and backgrounds (Task 7.4.2/7.4.3).

#![allow(
    clippy::unwrap_used,
    clippy::expect_used,
    clippy::panic,
    clippy::indexing_slicing
)]

use xtalforge_analysis::{Background, Profile, ProfileKind};

#[test]
fn profiles_are_normalised_to_their_height() {
    let gaussian = Profile {
        kind: ProfileKind::Gaussian,
        center: 20.0,
        fwhm: 0.2,
        eta: 0.0,
        height: 10.0,
    };
    assert!((gaussian.value(20.0) - 10.0).abs() < 1e-12);
    assert!((gaussian.value(20.0 + 0.1) - 5.0).abs() < 1e-9);
    let lorentzian = Profile {
        kind: ProfileKind::Lorentzian,
        ..gaussian
    };
    assert!((lorentzian.value(20.0) - 10.0).abs() < 1e-12);
    assert!((lorentzian.value(20.1) - 5.0).abs() < 1e-9);
    let voigt = Profile {
        kind: ProfileKind::PseudoVoigt,
        eta: 0.5,
        ..gaussian
    };
    assert!((voigt.value(20.1) - 5.0).abs() < 1e-9);
}

#[test]
fn sampling_and_backgrounds() {
    let profile = Profile {
        kind: ProfileKind::Gaussian,
        center: 1.0,
        fwhm: 0.1,
        eta: 0.0,
        height: 1.0,
    };
    let points = profile.sample([0.0, 2.0], 0.5).unwrap();
    assert_eq!(points.len(), 5);
    assert!(profile.sample([2.0, 0.0], 0.5).is_err());
    assert!(profile.sample([0.0, 2.0], 0.0).is_err());
    assert_eq!(Background::Constant(2.0).value(9.0), 2.0);
    assert_eq!(Background::Linear(1.0, 2.0).value(3.0), 7.0);
    assert_eq!(Background::Quadratic(1.0, 2.0, 3.0).value(2.0), 17.0);
}
