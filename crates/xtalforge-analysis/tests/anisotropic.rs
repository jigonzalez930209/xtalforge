//! Integration tests for anisotropic displacement parameters in the
//! structure factors and Rietveld model (Task 14.1.2).

#![allow(
    clippy::unwrap_used,
    clippy::expect_used,
    clippy::panic,
    clippy::indexing_slicing
)]

use xtalforge_analysis::{
    InstrumentProfile, LineShape, PxrdOptions, PxrdProfileOptions, Radiation, RadiationSource,
    RefineOptions, RietveldOptions, RietveldProfileOptions, StructureFactor,
    StructureFactorOptions, calculate_pxrd, calculate_pxrd_profile, fit_phase_profile,
    refine_phase, structure_factors,
};
use xtalforge_core::{Lattice, Site, Structure, StructureBuilder};
use xtalforge_io::xy::PatternPoint;

const A: f64 = 4.0;
const C: f64 = 6.0;

fn tetragonal(u_iso: Option<f64>, u_aniso: Option<[f64; 6]>) -> Structure {
    let mut site = Site::new("Si1", "Si", [0.0, 0.0, 0.0]).unwrap();
    if let Some(value) = u_iso {
        site = site.with_u_iso(value).unwrap();
    }
    if let Some(value) = u_aniso {
        site = site.with_u_aniso(value).unwrap();
    }
    StructureBuilder::new(
        "tetragonal",
        Lattice::from_parameters(A, A, C, 90.0, 90.0, 90.0).unwrap(),
    )
    .site(site)
    .build()
    .unwrap()
}

fn factors(structure: &Structure) -> Vec<StructureFactor> {
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

fn amplitude(factors: &[StructureFactor], hkl: [i32; 3]) -> f64 {
    factors
        .iter()
        .find(|reflection| reflection.hkl == hkl)
        .unwrap_or_else(|| panic!("missing {hkl:?}"))
        .amplitude_squared
}

#[test]
fn isotropic_anisotropic_parameters_reduce_to_the_same_factor() {
    let isotropic = factors(&tetragonal(Some(0.012), None));
    let anisotropic = factors(&tetragonal(
        None,
        Some([0.012, 0.012, 0.012, 0.0, 0.0, 0.0]),
    ));
    assert_eq!(isotropic.len(), anisotropic.len());
    for (first, second) in isotropic.iter().zip(anisotropic.iter()) {
        assert_eq!(first.hkl, second.hkl);
        let scale = first.amplitude_squared.abs().max(1.0);
        assert!(
            (first.amplitude_squared - second.amplitude_squared).abs() < 1e-10 * scale,
            "{} {} vs {}",
            first.label(),
            first.amplitude_squared,
            second.amplitude_squared
        );
    }
}

#[test]
fn anisotropic_tensor_changes_only_the_affected_reflections() {
    let reference = factors(&tetragonal(Some(0.02), None));
    let anisotropic = factors(&tetragonal(None, Some([0.02, 0.02, 0.005, 0.0, 0.0, 0.0])));
    let c_star = 1.0 / C;
    let expected_ratio =
        (4.0 * std::f64::consts::PI.powi(2) * c_star * c_star * (0.02 - 0.005)).exp();
    let ratio = amplitude(&anisotropic, [0, 0, 1]) / amplitude(&reference, [0, 0, 1]);
    assert!(
        (ratio - expected_ratio).abs() < 1e-10,
        "001 ratio {ratio} != {expected_ratio}"
    );
    let in_plane = amplitude(&anisotropic, [1, 0, 0]) / amplitude(&reference, [1, 0, 0]);
    assert!((in_plane - 1.0).abs() < 1e-12, "100 ratio {in_plane}");
}

fn synthetic(structure: &Structure, fwhm: f64) -> Vec<PatternPoint> {
    let profile = InstrumentProfile {
        radiation: RadiationSource::from_radiation(Radiation::CuKa).unwrap(),
        shape: LineShape::Gaussian,
        fwhm,
        ..InstrumentProfile::default()
    };
    let pattern = calculate_pxrd_profile(
        structure,
        &PxrdProfileOptions {
            profile,
            two_theta_range: [10.0, 80.0],
            step: 0.02,
        },
    )
    .unwrap();
    pattern
        .points
        .into_iter()
        .map(|point| PatternPoint {
            two_theta: point.two_theta,
            intensity: point.intensity + 2.0,
            uncertainty: None,
        })
        .collect()
}

#[test]
fn rietveld_profile_model_uses_the_anisotropic_tensor() {
    let truth = tetragonal(None, Some([0.02, 0.02, 0.005, 0.0, 0.0, 0.0]));
    let observed = synthetic(&truth, 0.12);
    let options = RietveldProfileOptions {
        base: RietveldOptions {
            radiation: Radiation::CuKa,
            two_theta_range: [10.0, 80.0],
            step: 0.02,
            fwhm: 0.12,
            ..RietveldOptions::default()
        },
        profile: InstrumentProfile {
            radiation: RadiationSource::from_radiation(Radiation::CuKa).unwrap(),
            shape: LineShape::Gaussian,
            fwhm: 0.12,
            ..InstrumentProfile::default()
        },
        refine_eta: false,
        refine_width: false,
        ..RietveldProfileOptions::default()
    };
    let correct = fit_phase_profile(&truth, &observed, &options).unwrap();
    assert!(
        correct.phase.rwp < 1e-3,
        "correct model Rwp {}",
        correct.phase.rwp
    );
    let wrong = tetragonal(Some(0.02), None);
    let wrong_fit = fit_phase_profile(&wrong, &observed, &options).unwrap();
    assert!(
        wrong_fit.phase.rwp > correct.phase.rwp * 10.0,
        "wrong model Rwp {} correct {}",
        wrong_fit.phase.rwp,
        correct.phase.rwp
    );
}

#[test]
fn the_adp_multiplier_scales_every_anisotropic_component() {
    let truth = tetragonal(None, Some([0.02, 0.02, 0.02, 0.001, 0.0, 0.0]));
    let start = tetragonal(None, Some([0.005, 0.005, 0.005, 0.000_25, 0.0, 0.0]));
    let observed = synthetic(&truth, 0.12);
    let result = refine_phase(
        &start,
        &observed,
        &RietveldOptions {
            radiation: Radiation::CuKa,
            two_theta_range: [10.0, 80.0],
            step: 0.02,
            fwhm: 0.12,
            ..RietveldOptions::default()
        },
        &RefineOptions {
            refine_adps: true,
            max_iterations: 30,
            ..RefineOptions::default()
        },
    )
    .unwrap();
    let scale = result
        .parameters
        .iter()
        .find(|(key, _)| key == "adp_scale")
        .map(|(_, value)| *value)
        .unwrap();
    assert!(
        (scale - 4.0).abs() < 0.05,
        "adp_scale {scale} with Rwp {}",
        result.rwp
    );
    let site = &result.structure.sites()[0];
    let aniso = site.u_aniso().unwrap();
    assert!((aniso[0] - 0.02).abs() < 1e-3);
    assert!((aniso[3] - 0.001).abs() < 1e-4);
    assert!(result.rwp < 1e-3, "Rwp {}", result.rwp);
}

#[test]
fn legacy_pxrd_options_still_apply_isotropic_parameters() {
    let structure = tetragonal(Some(0.02), None);
    let pattern = calculate_pxrd(
        &structure,
        &PxrdOptions {
            radiation: Radiation::CuKa,
            two_theta_range: [10.0, 80.0],
            step: 0.02,
            fwhm: 0.12,
        },
    )
    .unwrap();
    let zero = calculate_pxrd(
        &tetragonal(None, None),
        &PxrdOptions {
            radiation: Radiation::CuKa,
            two_theta_range: [10.0, 80.0],
            step: 0.02,
            fwhm: 0.12,
        },
    )
    .unwrap();
    let with_adp = pattern
        .reflections
        .iter()
        .find(|reflection| reflection.hkl == [0, 0, 1])
        .unwrap();
    let without = zero
        .reflections
        .iter()
        .find(|reflection| reflection.hkl == [0, 0, 1])
        .unwrap();
    assert!(with_adp.intensity < without.intensity);
}
