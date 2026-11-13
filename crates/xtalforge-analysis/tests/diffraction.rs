//! Integration tests for single-crystal structure factors, HKL export and
//! electron diffraction (Task 14.1.4).

#![allow(
    clippy::unwrap_used,
    clippy::expect_used,
    clippy::panic,
    clippy::indexing_slicing
)]

use xtalforge_analysis::{
    ElectronDiffractionOptions, PxrdOptions, Radiation, RadiationSource, StructureFactor,
    StructureFactorOptions, calculate_pxrd, electron_diffraction, hkl_text,
    neutron_scattering_length, structure_factors,
};
use xtalforge_core::{Lattice, Site, Structure, StructureBuilder};

fn silicon() -> Structure {
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

fn xray_factors(structure: &Structure) -> Vec<StructureFactor> {
    structure_factors(
        structure,
        &StructureFactorOptions {
            radiation: RadiationSource::xray(1.540_598).unwrap(),
            two_theta_range: [5.0, 80.0],
            min_amplitude: 0.0,
        },
    )
    .unwrap()
}

fn find(factors: &[StructureFactor], hkl: [i32; 3]) -> &StructureFactor {
    factors
        .iter()
        .find(|reflection| reflection.hkl == hkl)
        .unwrap_or_else(|| panic!("missing {hkl:?}"))
}

#[test]
fn single_crystal_amplitudes_match_the_powder_engine() {
    let structure = silicon();
    let factors = xray_factors(&structure);
    let pattern = calculate_pxrd(
        &structure,
        &PxrdOptions {
            radiation: Radiation::CuKa1,
            two_theta_range: [5.0, 80.0],
            step: 0.02,
            fwhm: 0.1,
        },
    )
    .unwrap();
    let reflection = pattern
        .reflections
        .iter()
        .find(|reflection| reflection.hkl == [1, 1, 1])
        .unwrap();
    let factor = find(&factors, [1, 1, 1]);
    let theta = 0.5 * reflection.two_theta.to_radians();
    let denominator = theta.sin().powi(2) * theta.cos();
    let c2 = reflection.two_theta.to_radians().cos().powi(2);
    let lp = (1.0 + c2) / denominator;
    assert!((reflection.two_theta - factor.two_theta).abs() < 1e-12);
    assert!((reflection.d - factor.d).abs() < 1e-12);
    let expected = reflection.intensity / lp;
    assert!(
        (expected - factor.amplitude_squared).abs() < 1e-9 * expected,
        "{} vs {}",
        expected,
        factor.amplitude_squared
    );
    assert!(factor.real > 0.0);
}

#[test]
fn friedel_pairs_are_conjugate_and_absences_are_weak() {
    let factors = xray_factors(&silicon());
    let first = find(&factors, [1, 1, 1]);
    let second = find(&factors, [-1, -1, -1]);
    assert!((first.real - second.real).abs() < 1e-12);
    assert!((first.imaginary + second.imaginary).abs() < 1e-12);
    assert!((first.amplitude_squared - second.amplitude_squared).abs() < 1e-15);

    let absent = factors
        .iter()
        .find(|reflection| reflection.hkl == [2, 0, 0]);
    assert!(
        absent.is_none_or(|reflection| reflection.amplitude_squared < 1e-6),
        "diamond 200 is absent"
    );
}

#[test]
fn hkl_text_has_one_hklf4_row_per_reflection() {
    let factors = xray_factors(&silicon());
    let text = hkl_text(&factors);
    assert_eq!(text.lines().count(), factors.len());
    let first_line = text.lines().next().unwrap();
    let columns: Vec<&str> = first_line.split_whitespace().collect();
    assert_eq!(columns.len(), 5);
    let hkl: Vec<i32> = columns[..3]
        .iter()
        .map(|value| value.parse().unwrap())
        .collect();
    assert_eq!(hkl, factors[0].hkl);
    let amplitude: f64 = columns[3].parse().unwrap();
    assert!((amplitude - factors[0].amplitude_squared).abs() < 1e-6);
    let sigma: f64 = columns[4].parse().unwrap();
    assert!((sigma - factors[0].amplitude()).abs() < 1e-6);
    for line in text.lines() {
        assert_eq!(line.split_whitespace().count(), 5);
    }
    let path = std::env::temp_dir().join(format!(
        "xtalforge-hkl-{}-{}.hkl",
        std::process::id(),
        "silicon"
    ));
    xtalforge_analysis::write_hkl(&path, &factors).unwrap();
    assert_eq!(std::fs::read_to_string(&path).unwrap(), text);
    std::fs::remove_file(&path).unwrap();
}

#[test]
fn form_factors_and_thresholds_validate() {
    let low = xtalforge_analysis::xray_form_factor("Si", 0.2).unwrap();
    let high = xtalforge_analysis::xray_form_factor("Si", 0.6).unwrap();
    assert!(low > high);
    assert!(high > 0.0);
    assert!(xtalforge_analysis::xray_form_factor("Si", -0.1).is_err());
    assert!(xtalforge_analysis::xray_form_factor("Qq", 0.1).is_err());

    let structure = silicon();
    let bad_threshold = StructureFactorOptions {
        min_amplitude: -1.0,
        ..StructureFactorOptions::default()
    };
    assert_eq!(
        structure_factors(&structure, &bad_threshold)
            .unwrap_err()
            .code(),
        "invalid_input"
    );
}

#[test]
fn neutron_structure_factors_use_the_scattering_length() {
    let structure = StructureBuilder::new("Si", Lattice::cubic(5.431).unwrap())
        .site(Site::new("Si1", "Si", [0.0, 0.0, 0.0]).unwrap())
        .build()
        .unwrap();
    let factors = structure_factors(
        &structure,
        &StructureFactorOptions {
            radiation: RadiationSource::neutron(1.0).unwrap(),
            two_theta_range: [10.0, 40.0],
            min_amplitude: 0.0,
        },
    )
    .unwrap();
    let expected = neutron_scattering_length("Si").unwrap().powi(2);
    let factor = find(&factors, [1, 0, 0]);
    assert!(
        (factor.amplitude_squared - expected).abs() < 1e-12 * expected,
        "{} vs {expected}",
        factor.amplitude_squared
    );
    assert!(factor.imaginary.abs() < 1e-20);
}

#[test]
fn electron_patterns_report_wavelength_and_ring_radius() {
    let pattern = electron_diffraction(
        &silicon(),
        &ElectronDiffractionOptions {
            acceleration_kv: 200.0,
            two_theta_max: 3.0,
            camera_length_mm: Some(100.0),
            min_intensity: 0.0,
        },
    )
    .unwrap();
    assert!((pattern.wavelength - 0.025_08).abs() < 2e-4);
    assert!(!pattern.reflections.is_empty());
    for reflection in &pattern.reflections {
        assert!(reflection.d > 0.0);
        assert!((reflection.g - 1.0 / reflection.d).abs() < 1e-12);
        let radius = reflection
            .ring_radius_mm
            .expect("a camera length was requested");
        let expected = 100.0 * reflection.two_theta.to_radians().tan();
        assert!(
            (radius - expected).abs() < 1e-9,
            "radius {radius} vs {expected}"
        );
        assert!(reflection.intensity > 0.0);
    }
    assert_eq!(pattern.camera_length_mm, Some(100.0));
    assert_eq!(pattern.acceleration_kv, 200.0);

    let no_camera = electron_diffraction(
        &silicon(),
        &ElectronDiffractionOptions {
            two_theta_max: 1.0,
            ..ElectronDiffractionOptions::default()
        },
    )
    .unwrap();
    assert!(
        no_camera
            .reflections
            .iter()
            .all(|reflection| reflection.ring_radius_mm.is_none())
    );
}

#[test]
fn invalid_requests_return_typed_errors() {
    let structure = silicon();
    let bad_range = StructureFactorOptions {
        radiation: RadiationSource::default(),
        two_theta_range: [80.0, 5.0],
        min_amplitude: 0.0,
    };
    assert_eq!(
        structure_factors(&structure, &bad_range)
            .unwrap_err()
            .code(),
        "invalid_input"
    );
    let bad_angle = ElectronDiffractionOptions {
        two_theta_max: 0.0,
        ..ElectronDiffractionOptions::default()
    };
    assert!(electron_diffraction(&structure, &bad_angle).is_err());
    let bad_voltage = ElectronDiffractionOptions {
        acceleration_kv: -200.0,
        ..ElectronDiffractionOptions::default()
    };
    assert!(electron_diffraction(&structure, &bad_voltage).is_err());
    let missing_element = StructureBuilder::new("X", Lattice::cubic(5.0).unwrap())
        .site(Site::new("X1", "X", [0.0, 0.0, 0.0]).unwrap())
        .build()
        .unwrap();
    let neutron = StructureFactorOptions {
        radiation: RadiationSource::neutron(1.0).unwrap(),
        two_theta_range: [10.0, 40.0],
        min_amplitude: 0.0,
    };
    assert!(structure_factors(&missing_element, &neutron).is_err());
}
