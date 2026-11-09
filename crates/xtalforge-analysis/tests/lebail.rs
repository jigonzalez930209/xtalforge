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
    assert!(result.rwp.is_finite() && result.rwp >= 0.0);
    assert!(result.chi_squared.is_finite() && result.chi_squared >= 0.0);
    assert_eq!(result.intensity_esds.len(), result.intensities.len());
    assert!(
        result
            .intensity_esds
            .iter()
            .all(|value| value.is_finite() && *value >= 0.0)
    );
    assert!(result.intensity_esds.iter().any(|value| *value > 0.0));
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
    assert_eq!(result.intensity_esds.len(), result.intensities.len());
    assert!(result.intensity_esds.iter().any(|value| *value > 0.0));
}

#[test]
fn le_bail_intensity_esds_cover_a_noisy_pattern() {
    use xtalforge_analysis::{
        LebailOptions, Profile, ProfileKind, PxrdOptions, Radiation, calculate_pxrd,
    };
    let structure = silicon();
    let pattern = calculate_pxrd(
        &structure,
        &PxrdOptions {
            radiation: Radiation::CuKa,
            two_theta_range: [10.0, 80.0],
            step: 0.02,
            fwhm: 0.1,
        },
    )
    .unwrap();
    let mut seed = 0x1234_5678_u64;
    let mut noise = |scale: f64| {
        let mut sum = 0.0;
        for _ in 0..12 {
            seed = seed
                .wrapping_mul(6_364_136_223_846_793_005)
                .wrapping_add(1_442_695_040_888_963_407);
            let value = ((seed >> 11) as f64) / ((1_u64 << 53) as f64);
            sum += value;
        }
        (sum - 6.0) * scale
    };
    let points: Vec<PatternPoint> = (0..3501)
        .map(|index| {
            let two_theta = 10.0 + index as f64 * 0.02;
            let mut intensity = 5.0;
            for reflection in &pattern.reflections {
                intensity += Profile {
                    kind: ProfileKind::Gaussian,
                    center: reflection.two_theta,
                    fwhm: 0.1,
                    eta: 0.5,
                    height: reflection.intensity,
                }
                .value(two_theta);
            }
            let noisy = intensity + noise(intensity.sqrt());
            PatternPoint {
                two_theta,
                intensity: noisy.max(0.0),
                uncertainty: None,
            }
        })
        .collect();
    let result = le_bail(
        &structure,
        &points,
        &LebailOptions {
            kind: ProfileKind::Gaussian,
            fwhm: 0.1,
            ..LebailOptions::default()
        },
    )
    .unwrap();
    assert_eq!(result.intensities.len(), pattern.reflections.len());
    let mut covered = 0;
    for (index, reflection) in pattern.reflections.iter().enumerate() {
        let found = result.intensities[index];
        let esd = result.intensity_esds[index];
        assert!(esd > 0.0, "reflection {index} has no esd");
        if (found - reflection.intensity).abs() <= 4.0 * esd {
            covered += 1;
        }
    }
    assert!(
        covered * 10 >= result.intensities.len() * 9,
        "only {covered} of {} intensities are within four esds",
        result.intensities.len()
    );
}
