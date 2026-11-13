//! Integration tests for the instrument profile, radiation and presets
//! (Task 14.1.1).

#![allow(
    clippy::unwrap_used,
    clippy::expect_used,
    clippy::panic,
    clippy::indexing_slicing
)]

use xtalforge_analysis::{
    Caglioti, FcjAsymmetry, INSTRUMENT_PRESET_NAMES, InstrumentProfile, LineShape, Monochromator,
    PxrdOptions, PxrdProfileOptions, Radiation, RadiationKind, RadiationSource, Tch,
    calculate_pxrd, calculate_pxrd_profile, electron_scattering_factor, electron_wavelength,
    instrument_presets, neutron_scattering_length, preset_by_name, tch_preset_by_name,
};
use xtalforge_core::{Lattice, Site, StructureBuilder};

fn silicon() -> xtalforge_core::Structure {
    let positions = [
        [0.0, 0.0, 0.0],
        [0.0, 0.5, 0.5],
        [0.5, 0.0, 0.5],
        [0.5, 0.5, 0.0],
        [0.25, 0.25, 0.25],
        [0.25, 0.75, 0.75],
        [0.75, 0.25, 0.75],
        [0.75, 0.75, 0.25],
    ];
    let mut builder = StructureBuilder::new("Si", Lattice::cubic(5.431).unwrap());
    for (index, frac) in positions.iter().enumerate() {
        builder = builder.site(Site::new(format!("Si{}", index + 1), "Si", *frac).unwrap());
    }
    builder.build().unwrap()
}

#[test]
fn caglioti_and_tch_widths_match_the_documented_formulas() {
    let caglioti = Caglioti {
        u: 0.004,
        v: -0.002,
        w: 0.012,
    };
    let two_theta = 30.0_f64;
    let tangent = (0.5 * two_theta.to_radians()).tan();
    let expected = (0.004 * tangent * tangent - 0.002 * tangent + 0.012).sqrt();
    assert!((caglioti.fwhm(two_theta) - expected).abs() < 1e-12);
    assert!(caglioti.fwhm(120.0) > caglioti.fwhm(20.0));

    let tch = Tch {
        gaussian: caglioti,
        x: 0.002,
        y: 0.004,
    };
    let gaussian = caglioti.fwhm(two_theta);
    let lorentzian = 0.002 * tangent + 0.004 / (0.5 * two_theta.to_radians()).cos();
    let ratio = lorentzian / gaussian;
    let polynomial = 1.0
        + 2.692_69 * ratio
        + 2.428_43 * ratio.powi(2)
        + 4.471_63 * ratio.powi(3)
        + 0.078_42 * ratio.powi(4)
        + ratio.powi(5);
    let expected_fwhm = gaussian * polynomial.powf(0.2);
    let expected_eta = 1.366_03 * ratio - 0.477_19 * ratio.powi(2) + 0.111_16 * ratio.powi(3);
    let (fwhm, eta) = tch.width(two_theta);
    assert!((fwhm - expected_fwhm).abs() < 1e-12);
    assert!((eta - expected_eta).abs() < 1e-12);
    assert!((0.0..=1.0).contains(&eta));
}

#[test]
fn profiles_are_normalised_and_obey_the_half_width_definition() {
    let gaussian = InstrumentProfile {
        shape: LineShape::Gaussian,
        fwhm: 0.2,
        ..InstrumentProfile::default()
    };
    assert!((gaussian.value(30.0, 30.0) - 1.0).abs() < 1e-12);
    assert!((gaussian.value(30.0, 30.1) - 0.5).abs() < 1e-12);
    assert!((gaussian.value(30.0, 29.9) - 0.5).abs() < 1e-12);

    let lorentzian = InstrumentProfile {
        shape: LineShape::Lorentzian,
        fwhm: 0.2,
        ..InstrumentProfile::default()
    };
    assert!((lorentzian.value(30.0, 30.1) - 0.5).abs() < 1e-12);

    let mut area = 0.0;
    let step = 0.0005;
    let mut x = 30.0 - 8.0 * 0.2;
    while x <= 30.0 + 8.0 * 0.2 {
        area += gaussian.area_normalised(30.0, x) * step;
        x += step;
    }
    assert!((area - 1.0).abs() < 1e-4, "area {area}");
}

#[test]
fn asymmetry_enhances_only_the_low_angle_tail() {
    let profile = InstrumentProfile {
        shape: LineShape::Lorentzian,
        fwhm: 0.4,
        asymmetry: Some(FcjAsymmetry {
            fcj1: 0.3,
            fcj2: 0.2,
        }),
        ..InstrumentProfile::default()
    };
    let center = 40.0;
    let low = profile.value(center, center - 0.2);
    let high = profile.value(center, center + 0.2);
    assert!(low > high, "low {low} high {high}");
    assert!((profile.value(center, center) - 1.0).abs() < 1e-12);
    let plain = InstrumentProfile {
        asymmetry: None,
        ..profile
    };
    assert!((plain.value(center, center - 0.2) - plain.value(center, center + 0.2)).abs() < 1e-12);
    assert!(
        FcjAsymmetry {
            fcj1: -0.1,
            fcj2: 0.0
        }
        .validate()
        .is_err()
    );
}

#[test]
fn polarisation_and_lorentz_factor_follow_the_hand_formula() {
    let unmonochromated = InstrumentProfile::default();
    let two_theta = 30.0_f64;
    let theta = 0.5 * two_theta.to_radians();
    let denom = theta.sin().powi(2) * theta.cos();
    let c2 = two_theta.to_radians().cos().powi(2);
    let expected = (1.0 + c2) / denom;
    assert!((unmonochromated.lorentz_polarization(two_theta) - expected).abs() < 1e-12);

    let polarised = InstrumentProfile {
        monochromator: Monochromator { polarization: 1.0 },
        ..unmonochromated
    };
    let expected = 2.0 * c2 / denom;
    assert!((polarised.lorentz_polarization(two_theta) - expected).abs() < 1e-12);
    assert!(
        polarised.lorentz_polarization(two_theta) < unmonochromated.lorentz_polarization(two_theta)
    );

    let neutron = InstrumentProfile {
        radiation: RadiationSource::neutron(1.5).unwrap(),
        ..unmonochromated
    };
    assert!((neutron.lorentz_polarization(two_theta) - 1.0 / denom).abs() < 1e-12);
}

#[test]
fn profile_pattern_reproduces_the_historical_gaussian() {
    let structure = silicon();
    let options = PxrdOptions::default()
        .with_radiation(Radiation::CuKa1)
        .with_fwhm(0.15);
    let legacy = calculate_pxrd(&structure, &options).unwrap();
    let profile_options = PxrdProfileOptions::from_pxrd(&options).unwrap();
    let modern = calculate_pxrd_profile(&structure, &profile_options).unwrap();
    assert_eq!(legacy.reflections, modern.reflections);
    assert_eq!(legacy.points.len(), modern.points.len());
    for (old, new) in legacy.points.iter().zip(modern.points.iter()) {
        let scale = old.intensity.abs().max(1.0);
        assert!((old.intensity - new.intensity).abs() < 1e-9 * scale);
    }
}

#[test]
fn caglioti_patterns_broaden_with_angle_and_stay_finite() {
    let structure = silicon();
    let profile = InstrumentProfile {
        shape: LineShape::PseudoVoigt,
        fwhm: 0.1,
        caglioti: Some(Caglioti {
            u: 0.004,
            v: -0.002,
            w: 0.012,
        }),
        eta: 0.5,
        ..InstrumentProfile::default()
    };
    let options = PxrdProfileOptions {
        profile,
        two_theta_range: [10.0, 100.0],
        step: 0.01,
    };
    let pattern = calculate_pxrd_profile(&structure, &options).unwrap();
    assert!(!pattern.points.is_empty());
    assert!(
        pattern
            .points
            .iter()
            .all(|point| point.intensity.is_finite())
    );
    assert!(pattern.max_intensity() > 0.0);
    assert!(profile.fwhm_at(110.0) > profile.fwhm_at(20.0));
}

#[test]
fn neutron_lengths_and_electron_wavelengths_match_reference_values() {
    let silicon = neutron_scattering_length("Si").unwrap();
    assert!((silicon - 4.1491e-5).abs() < 1e-12);
    let hydrogen = neutron_scattering_length("H").unwrap();
    assert!(hydrogen < 0.0);
    assert!((hydrogen + 3.739e-5).abs() < 1e-12);
    assert!(neutron_scattering_length("Xx").is_err());

    let wavelength = electron_wavelength(200.0).unwrap();
    assert!(
        (wavelength - 0.025_08).abs() < 2e-4,
        "200 kV wavelength {wavelength}"
    );
    assert!(electron_wavelength(0.0).is_err());
    assert!(electron_wavelength(f64::NAN).is_err());
    let source = RadiationSource::electron(200.0).unwrap();
    assert_eq!(source.kind, RadiationKind::Electron);

    let low_angle = electron_scattering_factor("Si", 0.2).unwrap();
    let high_angle = electron_scattering_factor("Si", 0.6).unwrap();
    assert!(low_angle > high_angle);
    assert!(high_angle > 0.0);
    assert!(electron_scattering_factor("Si", -0.1).is_err());
}

#[test]
fn instrument_profile_and_tch_patterns_are_wired() {
    let instrument = xtalforge_analysis::Instrument {
        polarization: 0.8,
        ..xtalforge_analysis::Instrument::default()
    };
    let profile = instrument
        .profile(RadiationSource::xray(1.0).unwrap())
        .unwrap();
    assert_eq!(profile.monochromator.polarization, 0.8);
    assert!(
        instrument
            .profile(RadiationSource::neutron(1.5).unwrap())
            .is_ok()
    );
    assert!(
        instrument
            .profile(RadiationSource {
                kind: RadiationKind::XRay,
                wavelength: f64::NAN,
            })
            .is_err()
    );

    let tch = tch_preset_by_name("lab-cu-ka").unwrap();
    let pattern = calculate_pxrd_profile(
        &silicon(),
        &PxrdProfileOptions {
            profile: tch,
            two_theta_range: [10.0, 60.0],
            step: 0.02,
        },
    )
    .unwrap();
    assert!(!pattern.points.is_empty());
    assert!(
        pattern
            .points
            .iter()
            .all(|point| point.intensity.is_finite())
    );
    assert!(pattern.max_intensity() > 0.0);
}

#[test]
fn profile_validation_rejects_bad_values() {
    let valid = InstrumentProfile::default();
    assert!(valid.validate().is_ok());
    let bad = [
        InstrumentProfile { fwhm: 0.0, ..valid },
        InstrumentProfile { eta: 1.5, ..valid },
        InstrumentProfile {
            shape: LineShape::Tch,
            ..valid
        },
        InstrumentProfile {
            radiation: RadiationSource {
                kind: RadiationKind::XRay,
                wavelength: 0.0,
            },
            ..valid
        },
        InstrumentProfile {
            monochromator: Monochromator { polarization: 2.0 },
            ..valid
        },
        InstrumentProfile {
            caglioti: Some(Caglioti {
                u: f64::NAN,
                v: 0.0,
                w: 0.0,
            }),
            ..valid
        },
        InstrumentProfile {
            tch: Some(Tch {
                gaussian: Caglioti::default(),
                x: f64::NAN,
                y: 0.0,
            }),
            ..valid
        },
        InstrumentProfile {
            asymmetry: Some(FcjAsymmetry {
                fcj1: f64::NAN,
                fcj2: 0.0,
            }),
            ..valid
        },
    ];
    for profile in bad {
        assert!(profile.validate().is_err(), "{profile:?}");
    }
    assert!(RadiationSource::xray(0.0).is_err());
    assert!(RadiationSource::neutron(-1.0).is_err());
    assert!(RadiationSource::xray(f64::NAN).is_err());
    assert!(Monochromator { polarization: -0.1 }.validate().is_err());
    assert!(
        Caglioti {
            u: 0.0,
            v: 0.0,
            w: f64::NAN
        }
        .validate()
        .is_err()
    );
    assert!(
        Tch {
            gaussian: Caglioti::default(),
            x: f64::NAN,
            y: 0.0
        }
        .validate()
        .is_err()
    );
    assert!(tch_preset_by_name("unknown").is_err());

    let bad_options = [
        PxrdProfileOptions::default().with_step(0.0),
        PxrdProfileOptions::default().with_range(80.0, 5.0),
        PxrdProfileOptions::default().with_range(-1.0, 20.0),
        PxrdProfileOptions::default().with_step(-0.02),
    ];
    for options in bad_options {
        assert!(options.validate().is_err());
        assert!(
            calculate_pxrd_profile(&silicon(), &options)
                .unwrap_err()
                .code()
                == "invalid_input"
        );
    }
    let broken_profile = PxrdProfileOptions {
        profile: InstrumentProfile {
            fwhm: -1.0,
            ..InstrumentProfile::default()
        },
        ..PxrdProfileOptions::default()
    };
    assert!(calculate_pxrd_profile(&silicon(), &broken_profile).is_err());
}

#[test]
fn instrument_presets_are_complete_and_named() {
    let presets = instrument_presets();
    assert_eq!(presets.len(), INSTRUMENT_PRESET_NAMES.len());
    for (preset, name) in presets.iter().zip(INSTRUMENT_PRESET_NAMES.iter()) {
        assert_eq!(preset.name, *name);
        assert!(preset.validate().is_ok(), "{}", preset.name);
        assert!(!preset.label.is_empty());
    }
    let copper = preset_by_name("LAB-CU-KA").unwrap();
    assert!((copper.profile.radiation.wavelength - 1.541_84).abs() < 1e-12);
    assert_eq!(copper.profile.radiation.kind, RadiationKind::XRay);
    let tch = tch_preset_by_name("lab-cu-ka").unwrap();
    assert_eq!(tch.shape, LineShape::Tch);
    assert!(tch.tch.is_some());
    assert!(tch.fwhm_at(30.0) > 0.0);
    assert!(preset_by_name("unknown").is_err());
}
