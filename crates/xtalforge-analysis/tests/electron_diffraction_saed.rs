//! Integration tests for kinematic zone-axis electron diffraction
//! (Task 19.3.3).
//!
//! # Provenance of the hand calculation
//!
//! Face-centred copper is written by hand with `a = 3.615 Angstrom`. The
//! 200 kV wavelength is checked against the tabulated `0.02508 Angstrom`. The
//! expected spot positions are `L lambda g` with `L = 500 mm`, and the
//! expected intensities are `|4 f_e|^2` from the public
//! `electron_scattering_factor`; no reference diffraction image is used.

#![allow(
    clippy::unwrap_used,
    clippy::expect_used,
    clippy::panic,
    clippy::indexing_slicing
)]

use xtalforge_analysis::{
    SaedOptions, SaedPattern, SaedSpot, electron_scattering_factor, zone_axis_pattern,
};
use xtalforge_core::{Lattice, Site, Structure, StructureBuilder};

const CU_A: f64 = 3.615;

fn copper() -> Structure {
    let positions = [
        [0.0, 0.0, 0.0],
        [0.0, 0.5, 0.5],
        [0.5, 0.0, 0.5],
        [0.5, 0.5, 0.0],
    ];
    let mut builder = StructureBuilder::new("Cu", Lattice::cubic(CU_A).unwrap());
    for (index, frac) in positions.iter().enumerate() {
        builder = builder.site(Site::new(format!("Cu{}", index + 1), "Cu", *frac).unwrap());
    }
    builder.build().unwrap()
}

fn options() -> SaedOptions {
    SaedOptions {
        acceleration_kv: 200.0,
        zone_axis: [0, 0, 1],
        camera_length_mm: 500.0,
        field_of_view_mm: 20.0,
        image_size: 128,
        spot_sigma_pixels: 2.0,
        min_intensity: 0.0,
        central_intensity: 4.0,
        max_laue_index: 0,
    }
}

fn find(pattern: &SaedPattern, hkl: [i32; 3]) -> &SaedSpot {
    pattern
        .spots
        .iter()
        .find(|spot| spot.hkl == hkl)
        .unwrap_or_else(|| panic!("missing {hkl:?}"))
}

#[test]
fn copper_001_has_eight_centred_spots_with_hand_calculated_positions() {
    let pattern = zone_axis_pattern(&copper(), &options()).unwrap();
    assert!((pattern.wavelength - 0.025_08).abs() < 2e-4);
    for (actual, expected) in pattern.zone_normal.iter().zip([0.0, 0.0, 1.0]) {
        assert!((actual - expected).abs() < 1e-12);
    }
    for (actual, expected) in pattern.detector_x.iter().zip([1.0, 0.0, 0.0]) {
        assert!((actual - expected).abs() < 1e-12);
    }
    for (actual, expected) in pattern.detector_y.iter().zip([0.0, 1.0, 0.0]) {
        assert!((actual - expected).abs() < 1e-12);
    }
    assert_eq!(pattern.spots.len(), 8);
    assert!(pattern.spots.iter().all(|spot| spot.hkl[2] == 0));

    let expected: Vec<[i32; 3]> = vec![
        [2, 0, 0],
        [-2, 0, 0],
        [0, 2, 0],
        [0, -2, 0],
        [2, 2, 0],
        [2, -2, 0],
        [-2, 2, 0],
        [-2, -2, 0],
    ];
    let mut actual: Vec<[i32; 3]> = pattern.spots.iter().map(|spot| spot.hkl).collect();
    actual.sort_unstable();
    let mut expected = expected;
    expected.sort_unstable();
    assert_eq!(actual, expected);

    let scale = 500.0 * pattern.wavelength;
    for spot in &pattern.spots {
        let d = CU_A / f64::from(spot.hkl[0] * spot.hkl[0] + spot.hkl[1] * spot.hkl[1]).sqrt();
        assert!(
            (spot.d - d).abs() < 1e-12,
            "{:?}: {} vs {}",
            spot.hkl,
            spot.d,
            d
        );
        assert!((spot.g - 1.0 / d).abs() < 1e-12);
        assert!((spot.x_mm - scale * f64::from(spot.hkl[0]) / CU_A).abs() < 1e-12);
        assert!((spot.y_mm - scale * f64::from(spot.hkl[1]) / CU_A).abs() < 1e-12);
        let s = 0.5 / d;
        let form = electron_scattering_factor("Cu", s).unwrap();
        let expected_intensity = (4.0 * form).powi(2);
        assert!(
            (spot.intensity - expected_intensity).abs() < 1e-12 * expected_intensity,
            "{} vs {expected_intensity}",
            spot.intensity
        );
        assert!(spot.excitation_error > 0.0);
    }
    let strongest = find(&pattern, [2, 0, 0]);
    let weaker = find(&pattern, [2, 2, 0]);
    assert_eq!(strongest.label(), "(2 0 0)");
    assert!((strongest.radius_mm() - strongest.x_mm.abs()).abs() < 1e-12);
    assert!(strongest.intensity > weaker.intensity);
    assert!((pattern.max_intensity() - strongest.intensity).abs() < 1e-12);

    let centroid: [f64; 2] = pattern.spots.iter().fold([0.0, 0.0], |sum, spot| {
        [sum[0] + spot.x_mm, sum[1] + spot.y_mm]
    });
    assert!(centroid[0].abs() < 1e-12);
    assert!(centroid[1].abs() < 1e-12);
    for spot in &pattern.spots {
        let partner = find(&pattern, spot.hkl.map(|value| -value));
        assert!((spot.intensity - partner.intensity).abs() < 1e-12 * spot.intensity);
        assert!((spot.x_mm + partner.x_mm).abs() < 1e-12);
        assert!((spot.y_mm + partner.y_mm).abs() < 1e-12);
    }
}

#[test]
fn image_buffer_is_symmetric_deterministic_and_eight_bit() {
    let pattern = zone_axis_pattern(&copper(), &options()).unwrap();
    let image = pattern.image().unwrap();
    assert_eq!(image.width, 128);
    assert_eq!(image.height, 128);
    assert_eq!(image.pixels.len(), 128 * 128);
    assert!((image.max_intensity() - 1.0).abs() < 1e-6);
    assert!(image.pixel(64, 64).unwrap() > 0.9);
    assert!(image.pixel(108, 64).unwrap() > 0.2);

    for y in 0..image.height {
        for x in 0..image.width / 2 {
            let left = image.pixel(x, y).unwrap();
            let right = image.pixel(image.width - 1 - x, y).unwrap();
            assert!((left - right).abs() < 1e-5);
        }
    }
    for x in 0..image.width {
        for y in 0..image.height / 2 {
            let top = image.pixel(x, y).unwrap();
            let bottom = image.pixel(x, image.height - 1 - y).unwrap();
            assert!((top - bottom).abs() < 1e-5);
        }
    }
    let again = pattern.image().unwrap();
    assert_eq!(image, again);
    let bytes = image.to_u8();
    assert_eq!(bytes.len(), image.pixels.len());
    assert_eq!(bytes.iter().copied().max(), Some(255));
    assert!(image.pixel(128, 0).is_none());
}

#[test]
fn general_zone_axis_keeps_only_zero_laue_zone_spots() {
    let structure = StructureBuilder::new("Fe", Lattice::cubic(2.8665).unwrap())
        .site(Site::new("Fe1", "Fe", [0.0, 0.0, 0.0]).unwrap())
        .site(Site::new("Fe2", "Fe", [0.5, 0.5, 0.5]).unwrap())
        .build()
        .unwrap();
    let pattern = zone_axis_pattern(
        &structure,
        &SaedOptions {
            zone_axis: [1, 1, 1],
            ..options()
        },
    )
    .unwrap();
    assert!(!pattern.spots.is_empty());
    for spot in &pattern.spots {
        assert_eq!(spot.hkl[0] + spot.hkl[1] + spot.hkl[2], 0);
    }
    let norm = 1.0 / 3.0_f64.sqrt();
    for (actual, expected) in pattern.zone_normal.iter().zip([norm, norm, norm]) {
        assert!((actual - expected).abs() < 1e-12);
    }
    let centroid: [f64; 2] = pattern.spots.iter().fold([0.0, 0.0], |sum, spot| {
        [sum[0] + spot.x_mm, sum[1] + spot.y_mm]
    });
    assert!(centroid[0].abs() < 1e-12);
    assert!(centroid[1].abs() < 1e-12);
}

#[test]
fn invalid_electron_requests_return_typed_errors() {
    let structure = copper();
    let mut bad = options();
    bad.zone_axis = [0, 0, 0];
    assert_eq!(
        zone_axis_pattern(&structure, &bad).unwrap_err().code(),
        "invalid_input"
    );
    bad = options();
    bad.image_size = 8;
    assert!(zone_axis_pattern(&structure, &bad).is_err());
    bad = options();
    bad.spot_sigma_pixels = 0.0;
    assert!(zone_axis_pattern(&structure, &bad).is_err());
    bad = options();
    bad.field_of_view_mm = 0.0;
    assert!(zone_axis_pattern(&structure, &bad).is_err());
    bad = options();
    bad.camera_length_mm = 0.0;
    assert!(zone_axis_pattern(&structure, &bad).is_err());
    bad = options();
    bad.acceleration_kv = 0.0;
    assert!(zone_axis_pattern(&structure, &bad).is_err());
    bad = options();
    bad.max_laue_index = 101;
    assert!(zone_axis_pattern(&structure, &bad).is_err());
    bad = options();
    bad.central_intensity = -1.0;
    assert!(zone_axis_pattern(&structure, &bad).is_err());

    let huge = StructureBuilder::new("huge", Lattice::cubic(10_000.0).unwrap())
        .site(Site::new("Cu1", "Cu", [0.0, 0.0, 0.0]).unwrap())
        .build()
        .unwrap();
    assert!(zone_axis_pattern(&huge, &options()).is_err());
}
