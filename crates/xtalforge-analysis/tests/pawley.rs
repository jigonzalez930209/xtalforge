//! Integration tests for the Pawley fit (Task 7.6.2).

#![allow(
    clippy::unwrap_used,
    clippy::expect_used,
    clippy::panic,
    clippy::indexing_slicing
)]

use xtalforge_analysis::{
    PawleyOptions, ProfileKind, PxrdOptions, Radiation, calculate_pxrd, pawley_refine,
};
use xtalforge_core::{Lattice, Site, Structure, StructureBuilder};
use xtalforge_io::xy::PatternPoint;

const FWHM: f64 = 0.1;

fn caesium_chloride() -> Structure {
    StructureBuilder::new("CsCl", Lattice::cubic(4.0).unwrap())
        .site(Site::new("Cs1", "Cs", [0.0, 0.0, 0.0]).unwrap())
        .site(Site::new("Cl1", "Cl", [0.5, 0.5, 0.5]).unwrap())
        .build()
        .unwrap()
}

/// Synthetic pattern from the PXRD reflections plus a zero shift and a
/// quadratic background, sampled on a 0.02 degree grid.
fn synthetic(structure: &Structure, zero: f64, range: [f64; 2]) -> Vec<PatternPoint> {
    let pattern = calculate_pxrd(
        structure,
        &PxrdOptions {
            radiation: Radiation::CuKa,
            two_theta_range: range,
            step: 0.02,
            fwhm: FWHM,
        },
    )
    .unwrap();
    let normalisation = 2.0 / FWHM * (std::f64::consts::LN_2 / std::f64::consts::PI).sqrt();
    let count = ((range[1] - range[0]) / 0.02).floor() as usize + 1;
    let mut points = Vec::with_capacity(count);
    for index in 0..count {
        let two_theta = range[0] + index as f64 * 0.02;
        let mut intensity = 4.0 + 0.05 * (two_theta - range[0]);
        for reflection in &pattern.reflections {
            let offset = two_theta - (reflection.two_theta + zero);
            let scaled = offset / FWHM;
            intensity += reflection.intensity
                * normalisation
                * (-4.0 * std::f64::consts::LN_2 * scaled * scaled).exp();
        }
        points.push(PatternPoint {
            two_theta,
            intensity,
            uncertainty: None,
        });
    }
    points
}

#[test]
fn pawley_recovers_intensities_cell_and_zero_shift() {
    let structure = caesium_chloride();
    let observed = synthetic(&structure, 0.03, [10.0, 70.0]);
    let options = PawleyOptions {
        kind: ProfileKind::Gaussian,
        fwhm: FWHM,
        background_order: 2,
        refine_cell: true,
        refine_zero: true,
        ..PawleyOptions::default()
    };
    let result = pawley_refine(&structure, &observed, &options).unwrap();
    assert!(result.converged, "the fit did not converge");
    assert!(result.rwp < 0.05, "Rwp {}", result.rwp);
    assert!(
        (result.zero_shift - 0.03).abs() < 5e-3,
        "zero shift {}",
        result.zero_shift
    );
    assert!(
        (result.cell[0] - 4.0).abs() < 5e-3,
        "a = {}",
        result.cell[0]
    );
    assert!(result.cell_esds[0] > 0.0, "no cell esd");
    assert!(result.zero_shift_esd > 0.0, "no zero-shift esd");
    assert!(
        result
            .reflections
            .iter()
            .all(|reflection| reflection.esd > 0.0)
    );
    assert!(result.chi_squared.is_finite());
    assert!(result.rexp.is_finite() && result.rexp > 0.0);

    let reference = calculate_pxrd(
        &structure,
        &PxrdOptions {
            radiation: Radiation::CuKa,
            two_theta_range: [10.0, 70.0],
            step: 0.02,
            fwhm: FWHM,
        },
    )
    .unwrap();
    for (found, expected) in result.reflections.iter().zip(reference.reflections.iter()) {
        assert_eq!(found.hkl, expected.hkl);
        assert!(
            (found.intensity - expected.intensity).abs()
                <= 0.02 * expected.intensity.abs().max(1.0),
            "{:?}: {} != {}",
            found.hkl,
            found.intensity,
            expected.intensity
        );
    }
}

#[test]
fn pawley_profiles_size_and_strain() {
    let structure = caesium_chloride();
    // A width model with size 0.08 deg and strain 0.03 deg in tan(theta).
    let range = [10.0, 70.0];
    let pattern = calculate_pxrd(
        &structure,
        &PxrdOptions {
            radiation: Radiation::CuKa,
            two_theta_range: range,
            step: 0.02,
            fwhm: 0.1,
        },
    )
    .unwrap();
    let count = ((range[1] - range[0]) / 0.02).floor() as usize + 1;
    let mut observed = Vec::with_capacity(count);
    for index in 0..count {
        let two_theta = range[0] + index as f64 * 0.02;
        let mut intensity = 2.0;
        for reflection in &pattern.reflections {
            let theta = reflection.two_theta.to_radians() / 2.0;
            let width = (0.08_f64.powi(2) + (0.03 * theta.tan()).powi(2)).sqrt();
            let normalisation =
                2.0 / width * (std::f64::consts::LN_2 / std::f64::consts::PI).sqrt();
            let offset = two_theta - reflection.two_theta;
            let scaled = offset / width;
            intensity += reflection.intensity
                * normalisation
                * (-4.0 * std::f64::consts::LN_2 * scaled * scaled).exp();
        }
        observed.push(PatternPoint {
            two_theta,
            intensity,
            uncertainty: None,
        });
    }
    let options = PawleyOptions {
        kind: ProfileKind::Gaussian,
        fwhm: 0.1,
        refine_cell: true,
        refine_zero: false,
        refine_size: true,
        refine_strain: true,
        max_iterations: 30,
        ..PawleyOptions::default()
    };
    let result = pawley_refine(&structure, &observed, &options).unwrap();
    assert!(result.rwp < 0.05, "Rwp {}", result.rwp);
    assert!((result.size - 0.08).abs() < 0.02, "size {}", result.size);
    assert!(
        (result.strain - 0.03).abs() < 0.02,
        "strain {}",
        result.strain
    );
}

#[test]
fn pawley_validates_its_input() {
    let structure = caesium_chloride();
    let observed = synthetic(&structure, 0.0, [10.0, 30.0]);
    assert!(pawley_refine(&structure, &observed[..4], &PawleyOptions::default()).is_err());
    assert!(
        pawley_refine(
            &structure,
            &observed,
            &PawleyOptions {
                fwhm: 0.0,
                ..PawleyOptions::default()
            }
        )
        .is_err()
    );
    let mut reversed = observed.clone();
    reversed.reverse();
    assert!(pawley_refine(&structure, &reversed, &PawleyOptions::default()).is_err());
}
