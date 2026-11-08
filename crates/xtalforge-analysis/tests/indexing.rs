//! Integration tests for cubic indexing (Task 7.5).

#![allow(
    clippy::unwrap_used,
    clippy::expect_used,
    clippy::panic,
    clippy::indexing_slicing
)]

use xtalforge_analysis::{Peak, index_cubic};

#[test]
fn indexes_cubic_silicon_peaks() {
    let a = 5.430941;
    let wavelength = 1.540598;
    let hkls = [[1, 1, 1], [2, 2, 0], [3, 1, 1], [4, 0, 0], [3, 3, 1]];
    let peaks: Vec<Peak> = hkls
        .iter()
        .map(|hkl| {
            let sum = (hkl[0] * hkl[0] + hkl[1] * hkl[1] + hkl[2] * hkl[2]) as f64;
            let d = a / sum.sqrt();
            let two_theta = 2.0 * (wavelength / (2.0 * d)).asin().to_degrees();
            Peak {
                two_theta,
                intensity: 1.0,
                fwhm: 0.1,
            }
        })
        .collect();
    let index = index_cubic(&peaks, wavelength, 5, 0.05).unwrap().unwrap();
    assert!((index.a - a).abs() < 0.05);
    assert!(index.completeness >= 0.8);
    assert!(index.mean_error < 0.05);
}
