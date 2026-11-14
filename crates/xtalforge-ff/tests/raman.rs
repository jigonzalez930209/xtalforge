//! Integration tests for powder-averaged Raman spectra (Task 19.4.2).
//!
//! The activity checks use the invariant formula
//! `(45 alpha^2 + 7 beta^2) / 45` with tensors whose answer is known by hand:
//! an isotropic tensor gives `a^2`, and `diag(a, 0, 0)` gives `4 a^2 / 15`.

#![allow(
    clippy::unwrap_used,
    clippy::expect_used,
    clippy::panic,
    clippy::indexing_slicing
)]

use xtalforge_ff::{
    Broadening, RamanOptions, RamanTensorSet, VibrationalModes, powder_average, raman_spectrum,
};

fn modes(frequencies: Vec<f64>, vectors: Vec<Vec<[f64; 3]>>) -> VibrationalModes {
    let zero_modes = frequencies
        .iter()
        .enumerate()
        .filter(|(_, frequency)| frequency.abs() <= 5.0)
        .map(|(index, _)| index)
        .collect();
    let imaginary = frequencies
        .iter()
        .enumerate()
        .filter(|(_, frequency)| **frequency < -5.0)
        .map(|(index, _)| index)
        .collect();
    VibrationalModes {
        frequencies,
        modes: vectors,
        zero_modes,
        imaginary,
    }
}

fn one_atom_modes(frequencies: Vec<f64>) -> VibrationalModes {
    let vectors = frequencies
        .iter()
        .map(|_| vec![[1.0, 0.0, 0.0]])
        .collect::<Vec<_>>();
    modes(frequencies, vectors)
}

fn diagonal(values: [f64; 3]) -> [[f64; 3]; 3] {
    [
        [values[0], 0.0, 0.0],
        [0.0, values[1], 0.0],
        [0.0, 0.0, values[2]],
    ]
}

#[test]
fn isotropic_tensor_powder_average_is_the_squared_amplitude() {
    for amplitude in [0.0, 1.0, 2.5, -3.0] {
        let activity = powder_average(&diagonal([amplitude, amplitude, amplitude]));
        assert!(
            (activity - amplitude * amplitude).abs() < 1e-12,
            "{activity} != {}",
            amplitude * amplitude
        );
    }
}

#[test]
fn uniaxial_tensor_powder_average_matches_the_closed_form() {
    let activity = powder_average(&diagonal([1.0, 0.0, 0.0]));
    assert!((activity - 4.0 / 15.0).abs() < 1e-12, "{activity}");
    let scaled = powder_average(&diagonal([2.0, 0.0, 0.0]));
    assert!((scaled - 16.0 / 15.0).abs() < 1e-12, "{scaled}");
    let traceless = powder_average(&diagonal([1.0, -1.0, 0.0]));
    assert!((traceless - 7.0 / 15.0).abs() < 1e-12, "{traceless}");
}

#[test]
fn raman_spectrum_uses_the_real_mode_and_broadens_it() {
    let modes = one_atom_modes(vec![0.0, 0.0, 1000.0]);
    let tensors = RamanTensorSet::from_tensors(vec![diagonal([1.0, 1.0, 1.0])]);
    let spectrum = raman_spectrum(&modes, &tensors, &RamanOptions::default()).unwrap();

    assert_eq!(spectrum.peaks.len(), 1);
    assert_eq!(spectrum.peaks[0].mode, 2);
    assert!((spectrum.peaks[0].activity - 1.0).abs() < 1e-12);
    assert!(!spectrum.points.is_empty());
    assert!(spectrum.max_activity() >= 0.999_999);
    // The grid is 950, 955, ..., 1050 and the stick lands exactly on 1000.
    let peak = spectrum
        .points
        .iter()
        .min_by(|first, second| {
            (first[0] - 1000.0)
                .abs()
                .partial_cmp(&(second[0] - 1000.0).abs())
                .unwrap()
        })
        .unwrap();
    assert!((peak[0] - 1000.0).abs() < 1e-12);
    assert!((peak[1] - 1.0).abs() < 1e-12);
}

#[test]
fn per_mode_layout_keeps_the_mode_indices_and_skips_imaginary_modes() {
    let modes = one_atom_modes(vec![0.0, -20.0, 500.0, 800.0]);
    let tensors = RamanTensorSet::from_tensors(vec![
        diagonal([0.0, 0.0, 0.0]),
        diagonal([9.0, 9.0, 9.0]),
        diagonal([1.0, 1.0, 1.0]),
        diagonal([2.0, 2.0, 2.0]),
    ]);
    let spectrum = raman_spectrum(&modes, &tensors, &RamanOptions::default()).unwrap();
    assert_eq!(spectrum.peaks.len(), 2);
    assert_eq!(spectrum.peaks[0].mode, 2);
    assert_eq!(spectrum.peaks[1].mode, 3);
    assert!((spectrum.peaks[0].activity - 1.0).abs() < 1e-12);
    assert!((spectrum.peaks[1].activity - 4.0).abs() < 1e-12);
}

#[test]
fn gaussian_option_and_a_larger_width_change_the_curve() {
    let modes = one_atom_modes(vec![0.0, 0.0, 1000.0]);
    let tensors = RamanTensorSet::from_tensors(vec![diagonal([1.0, 1.0, 1.0])]);
    let lorentzian = raman_spectrum(&modes, &tensors, &RamanOptions::default()).unwrap();
    let gaussian = raman_spectrum(
        &modes,
        &tensors,
        &RamanOptions {
            broadening: Broadening::Gaussian,
            ..RamanOptions::default()
        },
    )
    .unwrap();
    assert_eq!(lorentzian.points.len(), gaussian.points.len());
    assert!(
        (lorentzian.points.iter().map(|point| point[1]).sum::<f64>()
            - gaussian.points.iter().map(|point| point[1]).sum::<f64>())
        .abs()
            > 1e-6
    );
    let wide = raman_spectrum(
        &modes,
        &tensors,
        &RamanOptions {
            fwhm: 40.0,
            ..RamanOptions::default()
        },
    )
    .unwrap();
    assert!(wide.points.len() > lorentzian.points.len());
}

#[test]
fn text_and_json_tensor_files_parse_identically() {
    let text = "\
        # two modes\n\
        1  1 0 0  0 1 0  0 0 1\n\
        2  2 0 0  0 2 0  0 0 2\n";
    let from_text = RamanTensorSet::from_text(text).unwrap();
    let json = r#"{"natoms": 2, "tensors": [
        [[1, 0, 0], [0, 1, 0], [0, 0, 1]],
        [[2, 0, 0], [0, 2, 0], [0, 0, 2]]
    ]}"#;
    let from_json = RamanTensorSet::from_text(json).unwrap();
    assert_eq!(from_text, from_json);
    assert_eq!(from_text.len(), 2);
    assert!(!from_text.is_empty());
}

#[test]
fn empty_and_mismatched_inputs_are_rejected() {
    let silent = one_atom_modes(vec![0.0, 0.0, 0.0]);
    let empty = RamanTensorSet::from_tensors(Vec::new());
    let spectrum = raman_spectrum(&silent, &empty, &RamanOptions::default()).unwrap();
    assert!(spectrum.points.is_empty());
    assert!(spectrum.peaks.is_empty());

    let active = one_atom_modes(vec![0.0, 0.0, 500.0]);
    let error = raman_spectrum(&active, &empty, &RamanOptions::default()).unwrap_err();
    assert_eq!(error.code(), "invalid_input");

    let wrong_count =
        RamanTensorSet::from_tensors(vec![diagonal([1.0, 1.0, 1.0]), diagonal([1.0, 1.0, 1.0])]);
    let error = raman_spectrum(&active, &wrong_count, &RamanOptions::default()).unwrap_err();
    assert_eq!(error.code(), "invalid_input");

    let mut not_finite = diagonal([1.0, 1.0, 1.0]);
    not_finite[0][0] = f64::INFINITY;
    let error = raman_spectrum(
        &active,
        &RamanTensorSet::from_tensors(vec![not_finite]),
        &RamanOptions::default(),
    )
    .unwrap_err();
    assert_eq!(error.code(), "invalid_input");

    let error = raman_spectrum(
        &active,
        &RamanTensorSet::from_tensors(vec![diagonal([1.0, 1.0, 1.0])]),
        &RamanOptions {
            fwhm: -1.0,
            ..RamanOptions::default()
        },
    )
    .unwrap_err();
    assert_eq!(error.code(), "invalid_input");

    let error = RamanTensorSet::from_text("1 0 0 0 0 0 0 0").unwrap_err();
    assert_eq!(error.code(), "invalid_input");

    let error = RamanTensorSet::from_text("{\"natoms\": 3, \"tensors\": []}").unwrap_err();
    assert_eq!(error.code(), "invalid_input");

    let error = RamanTensorSet::from_text("   ").unwrap_err();
    assert_eq!(error.code(), "invalid_input");
}
