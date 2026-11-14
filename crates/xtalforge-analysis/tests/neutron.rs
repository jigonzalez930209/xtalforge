//! Integration tests for nuclear neutron powder diffraction (Task 19.3.1).
//!
//! # Provenance of the hand calculations
//!
//! The cells are written by hand from the textbook rock-salt and CsCl
//! descriptions. The scattering lengths are the values of the bundled NIST
//! NCNR table (`crates/xtalforge-analysis/data/neutron.json`, public-domain
//! factual data): `b_Na = 3.630 fm`, `b_Cl = 9.577 fm`, `b_H = -3.739 fm`,
//! `b_Ti = -3.438 fm`; `1 fm = 1e-5 Angstrom`. Expected structure factors are
//! evaluated in the test from these numbers, not read from a reference file.

#![allow(
    clippy::unwrap_used,
    clippy::expect_used,
    clippy::panic,
    clippy::indexing_slicing
)]

use xtalforge_analysis::{
    NeutronAxis, NeutronGeometry, NeutronOptions, NeutronPattern, NeutronReflection,
    neutron_powder, neutron_scattering_length,
};
use xtalforge_core::{Lattice, Site, Structure, StructureBuilder};

const NACL_A: f64 = 5.6402;
const B_NA: f64 = 3.630e-5;
const B_CL: f64 = 9.577e-5;

fn site(label: &str, element: &str, frac: [f64; 3]) -> Site {
    Site::new(label, element, frac).unwrap()
}

fn neutron_profile(wavelength: f64) -> xtalforge_analysis::InstrumentProfile {
    xtalforge_analysis::InstrumentProfile {
        radiation: xtalforge_analysis::RadiationSource::neutron(wavelength).unwrap(),
        shape: xtalforge_analysis::LineShape::Gaussian,
        fwhm: 0.2,
        ..xtalforge_analysis::InstrumentProfile::default()
    }
}

/// Rock salt NaCl in P1: four Na and four Cl atoms of the conventional cell.
fn halite() -> Structure {
    let sodium = [
        [0.0, 0.0, 0.0],
        [0.0, 0.5, 0.5],
        [0.5, 0.0, 0.5],
        [0.5, 0.5, 0.0],
    ];
    let chlorine = [
        [0.5, 0.0, 0.0],
        [0.0, 0.5, 0.0],
        [0.0, 0.0, 0.5],
        [0.5, 0.5, 0.5],
    ];
    let mut builder = StructureBuilder::new("NaCl", Lattice::cubic(NACL_A).unwrap());
    for (index, frac) in sodium.iter().enumerate() {
        builder = builder.site(site(&format!("Na{}", index + 1), "Na", *frac));
    }
    for (index, frac) in chlorine.iter().enumerate() {
        builder = builder.site(site(&format!("Cl{}", index + 1), "Cl", *frac));
    }
    builder.build().unwrap()
}

/// CsCl-like H-Ti cell with H at the origin and Ti at the body centre.
fn hydrogen_titanium() -> Structure {
    StructureBuilder::new("HTi", Lattice::cubic(4.0).unwrap())
        .site(site("H1", "H", [0.0, 0.0, 0.0]))
        .site(site("Ti1", "Ti", [0.5, 0.5, 0.5]))
        .build()
        .unwrap()
}

fn find(pattern: &NeutronPattern, hkl: [i32; 3]) -> &NeutronReflection {
    pattern
        .reflections
        .iter()
        .find(|reflection| reflection.hkl == hkl)
        .unwrap_or_else(|| panic!("missing {hkl:?}"))
}

fn two_theta(wavelength: f64, d: f64) -> f64 {
    2.0 * (0.5 * wavelength / d).asin().to_degrees()
}

fn lorentz(two_theta: f64) -> f64 {
    let theta = 0.5 * two_theta.to_radians();
    1.0 / (theta.sin().powi(2) * theta.cos())
}

#[test]
fn rock_salt_positions_and_intensities_match_the_hand_calculation() {
    let wavelength = 1.5;
    let options = NeutronOptions {
        geometry: NeutronGeometry::ConstantWavelength { wavelength },
        profile: neutron_profile(wavelength),
        range: [10.0, 60.0],
        step: 0.02,
    };
    let pattern = neutron_powder(&halite(), &options).unwrap();

    let d111 = NACL_A / 3.0_f64.sqrt();
    let d200 = NACL_A / 2.0;
    let f111 = 4.0 * (B_NA - B_CL);
    let f200 = 4.0 * (B_NA + B_CL);
    let reflection_111 = find(&pattern, [1, 1, 1]);
    let reflection_200 = find(&pattern, [2, 0, 0]);

    assert!((reflection_111.d - d111).abs() < 1e-12);
    assert!((reflection_200.d - d200).abs() < 1e-12);
    assert!((reflection_111.two_theta - two_theta(wavelength, d111)).abs() < 1e-12);
    assert!((reflection_200.two_theta - two_theta(wavelength, d200)).abs() < 1e-12);
    assert!((reflection_111.nuclear_amplitude_squared - f111 * f111).abs() < 1e-12 * f111 * f111);
    assert!((reflection_200.nuclear_amplitude_squared - f200 * f200).abs() < 1e-12 * f200 * f200);
    assert!(
        (reflection_111.intensity - f111 * f111 * lorentz(reflection_111.two_theta)).abs()
            < 1e-12 * reflection_111.intensity
    );
    assert!(
        (reflection_200.intensity - f200 * f200 * lorentz(reflection_200.two_theta)).abs()
            < 1e-12 * reflection_200.intensity
    );

    // A P1 description lists each member of the orbit; their sum carries the
    // multiplicity of the conventional rock-salt reflection.
    let orbit_111: Vec<&NeutronReflection> = pattern
        .reflections
        .iter()
        .filter(|reflection| (reflection.d - d111).abs() < 1e-9)
        .collect();
    let orbit_200: Vec<&NeutronReflection> = pattern
        .reflections
        .iter()
        .filter(|reflection| (reflection.d - d200).abs() < 1e-9)
        .collect();
    assert_eq!(orbit_111.len(), 8);
    assert_eq!(orbit_200.len(), 6);
    assert_eq!(reflection_111.label(), "(1 1 1)");
    assert!(
        (pattern
            .normalised()
            .iter()
            .map(|point| point.intensity)
            .fold(0.0, f64::max)
            - 100.0)
            .abs()
            < 1e-9
    );
    let sum_111: f64 = orbit_111.iter().map(|value| value.intensity).sum();
    let sum_200: f64 = orbit_200.iter().map(|value| value.intensity).sum();
    let expected_111 = f111 * f111 * 8.0 * lorentz(two_theta(wavelength, d111));
    let expected_200 = f200 * f200 * 6.0 * lorentz(two_theta(wavelength, d200));
    assert!((sum_111 - expected_111).abs() < 1e-12 * expected_111);
    assert!((sum_200 - expected_200).abs() < 1e-12 * expected_200);

    assert_eq!(pattern.axis, NeutronAxis::TwoTheta { wavelength });
    assert_eq!(pattern.axis.unit(), "degrees 2theta");
    assert!(pattern.max_intensity() > 0.0);
}

#[test]
fn negative_scattering_lengths_change_the_phase() {
    let wavelength = 1.5;
    let pattern = neutron_powder(
        &hydrogen_titanium(),
        &NeutronOptions {
            geometry: NeutronGeometry::ConstantWavelength { wavelength },
            profile: neutron_profile(wavelength),
            range: [10.0, 70.0],
            step: 0.02,
        },
    )
    .unwrap();
    let b_h: f64 = -3.739e-5;
    let b_ti: f64 = -3.438e-5;
    assert!(neutron_scattering_length("H").unwrap() < 0.0);
    assert!(neutron_scattering_length("Ti").unwrap() < 0.0);

    // (100): the Ti phase is exp(i pi) = -1, so the two negative lengths
    // subtract. (110): both phases are +1 and the lengths add.
    let difference = b_h - b_ti;
    let sum = b_h + b_ti;
    let reflection_100 = find(&pattern, [1, 0, 0]);
    let reflection_110 = find(&pattern, [1, 1, 0]);
    assert!(
        (reflection_100.nuclear_amplitude_squared - difference * difference).abs()
            < 1e-14 * difference * difference
    );
    assert!((reflection_110.nuclear_amplitude_squared - sum * sum).abs() < 1e-12 * sum * sum);
    // If H scattered positively the (100) line would be about 568 times
    // stronger: this is the phase effect the test protects.
    let positive_hypothesis = (b_h.abs() + b_ti.abs()).powi(2);
    assert!(reflection_100.nuclear_amplitude_squared < 0.01 * positive_hypothesis);
}

#[test]
fn time_of_flight_axis_is_monotonic_and_calibrated() {
    let path: f64 = 20.0;
    let bank: f64 = 90.0;
    let relative = 0.002;
    let range = [1000.0, 30_000.0];
    let step = 20.0;
    let pattern = neutron_powder(
        &halite(),
        &NeutronOptions {
            geometry: NeutronGeometry::TimeOfFlight {
                total_path_m: path,
                two_theta_deg: bank,
                relative_d_fwhm: relative,
            },
            profile: neutron_profile(1.0),
            range,
            step,
        },
    )
    .unwrap();
    assert_eq!(pattern.axis.unit(), "microseconds");
    assert_eq!(pattern.points.len(), 1451);
    for pair in pattern.points.windows(2) {
        assert!(pair[1].x > pair[0].x);
        assert!((pair[1].x - pair[0].x - step).abs() < 1e-9);
    }
    let scale = 2.0 * 252.777 * path * (0.5 * bank.to_radians()).sin();
    for reflection in &pattern.reflections {
        let tof = reflection.tof_us.expect("time-of-flight reflection");
        assert!((tof - scale * reflection.d).abs() < 1e-9);
        assert!((reflection.two_theta - bank).abs() < 1e-12);
        assert!(tof >= range[0] - 1e-6 && tof <= range[1] + 1e-6);
    }
    let mut by_spacing: Vec<&NeutronReflection> = pattern.reflections.iter().collect();
    by_spacing.sort_by(|left, right| left.d.total_cmp(&right.d));
    for pair in by_spacing.windows(2) {
        assert!(pair[1].tof_us.unwrap() >= pair[0].tof_us.unwrap());
    }
    assert!(!pattern.reflections.is_empty());
    assert!(pattern.max_intensity() > 0.0);
}

#[test]
fn invalid_neutron_requests_return_typed_errors() {
    let structure = halite();
    let base = NeutronOptions {
        geometry: NeutronGeometry::ConstantWavelength { wavelength: 1.5 },
        profile: neutron_profile(1.5),
        range: [10.0, 60.0],
        step: 0.02,
    };
    let mut bad = base;
    bad.geometry = NeutronGeometry::ConstantWavelength { wavelength: 0.0 };
    assert_eq!(
        neutron_powder(&structure, &bad).unwrap_err().code(),
        "invalid_input"
    );
    bad = base;
    bad.range = [0.0, 60.0];
    assert!(neutron_powder(&structure, &bad).is_err());
    bad = base;
    bad.range = [60.0, 10.0];
    assert!(neutron_powder(&structure, &bad).is_err());
    bad = base;
    bad.step = 0.0;
    assert!(neutron_powder(&structure, &bad).is_err());
    bad = base;
    bad.profile = xtalforge_analysis::InstrumentProfile::default();
    assert!(neutron_powder(&structure, &bad).is_err());
    bad = base;
    bad.range = [5.0, 179.0];
    bad.step = 1e-6;
    assert!(neutron_powder(&structure, &bad).is_err());

    let bad_tof = NeutronOptions {
        geometry: NeutronGeometry::TimeOfFlight {
            total_path_m: 0.0,
            two_theta_deg: 90.0,
            relative_d_fwhm: 0.002,
        },
        profile: neutron_profile(1.0),
        range: [1000.0, 10_000.0],
        step: 10.0,
    };
    assert!(neutron_powder(&structure, &bad_tof).is_err());
    let bad_angle = NeutronOptions {
        geometry: NeutronGeometry::TimeOfFlight {
            total_path_m: 20.0,
            two_theta_deg: 180.0,
            relative_d_fwhm: 0.002,
        },
        ..bad_tof
    };
    assert!(neutron_powder(&structure, &bad_angle).is_err());
    let bad_resolution = NeutronOptions {
        geometry: NeutronGeometry::TimeOfFlight {
            total_path_m: 20.0,
            two_theta_deg: 90.0,
            relative_d_fwhm: 0.0,
        },
        ..bad_tof
    };
    assert!(neutron_powder(&structure, &bad_resolution).is_err());

    let dummy = StructureBuilder::new("X", Lattice::cubic(4.0).unwrap())
        .site(site("X1", "X", [0.0, 0.0, 0.0]))
        .build()
        .unwrap();
    assert_eq!(
        neutron_powder(&dummy, &base).unwrap_err().code(),
        "invalid_input"
    );
}
