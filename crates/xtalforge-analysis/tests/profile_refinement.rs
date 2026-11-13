//! Integration tests for refinable mixing and profile parameters in the
//! Pawley, Le Bail and Rietveld engines (Task 14.1.3).

#![allow(
    clippy::unwrap_used,
    clippy::expect_used,
    clippy::panic,
    clippy::indexing_slicing
)]

use xtalforge_analysis::{
    InstrumentProfile, LebailOptions, LebailProfileOptions, LineShape, PawleyOptions,
    PawleyProfileOptions, ProfileKind, PxrdOptions, Radiation, RadiationSource, RietveldOptions,
    RietveldProfileOptions, calculate_pxrd, fit_phase_profile, le_bail_profile,
    pawley_refine_profile,
};
use xtalforge_core::{Lattice, Site, Structure, StructureBuilder};
use xtalforge_io::xy::PatternPoint;

/// Pseudo-Voigt mixing of the synthetic patterns.
const TRUE_ETA: f64 = 0.7;
/// FWHM of the synthetic patterns in degrees.
const FWHM: f64 = 0.12;

fn cscl() -> Structure {
    StructureBuilder::new("CsCl", Lattice::cubic(4.0).unwrap())
        .site(Site::new("Cs1", "Cs", [0.0, 0.0, 0.0]).unwrap())
        .site(Site::new("Cl1", "Cl", [0.5, 0.5, 0.5]).unwrap())
        .build()
        .unwrap()
}

fn synthetic(eta: f64) -> Vec<PatternPoint> {
    synthetic_with(eta, false)
}

/// Synthetic pattern from the Pawley active window (six FWHM), which is the
/// model that `pawley_refine` evaluates.
fn synthetic_pawley(eta: f64) -> Vec<PatternPoint> {
    synthetic_with(eta, true)
}

fn synthetic_with(eta: f64, active_window: bool) -> Vec<PatternPoint> {
    let profile = InstrumentProfile {
        radiation: RadiationSource::from_radiation(Radiation::CuKa).unwrap(),
        shape: LineShape::PseudoVoigt,
        fwhm: FWHM,
        eta,
        ..InstrumentProfile::default()
    };
    let reflections = calculate_pxrd(
        &cscl(),
        &PxrdOptions {
            radiation: Radiation::CuKa,
            two_theta_range: [15.0, 70.0],
            step: 0.02,
            fwhm: FWHM,
        },
    )
    .unwrap()
    .reflections;
    let mut points = Vec::new();
    let mut x = 15.0;
    while x <= 70.0 {
        let mut intensity = 0.0;
        for reflection in &reflections {
            if active_window && (x - reflection.two_theta).abs() > 6.0 * FWHM {
                continue;
            }
            intensity += reflection.intensity * profile.area_normalised(reflection.two_theta, x);
        }
        points.push(PatternPoint {
            two_theta: x,
            intensity: intensity + 5.0,
            uncertainty: None,
        });
        x += 0.02;
    }
    points
}

fn fixed_pawley(eta: f64) -> PawleyOptions {
    PawleyOptions {
        kind: ProfileKind::PseudoVoigt,
        eta,
        fwhm: FWHM,
        background_order: 0,
        refine_cell: false,
        refine_zero: false,
        refine_size: false,
        refine_strain: false,
        max_iterations: 10,
        max_reflections: 200,
    }
}

#[test]
fn pawley_recovers_a_pseudo_voigt_mixing() {
    let observed = synthetic_pawley(TRUE_ETA);
    let options = PawleyProfileOptions {
        base: fixed_pawley(0.2),
        refine_eta: true,
        refine_width: false,
        eta_range: [0.0, 1.0],
        ..PawleyProfileOptions::default()
    };
    let result = pawley_refine_profile(&cscl(), &observed, &options).unwrap();
    assert!(
        (result.eta - TRUE_ETA).abs() < 0.02,
        "eta {} != {TRUE_ETA} with Rwp {}",
        result.eta,
        result.fit.rwp
    );
    assert!(result.fit.rwp < 1e-3, "Rwp {}", result.fit.rwp);
    assert!((result.fwhm - FWHM).abs() < 1e-12);

    let fixed = pawley_refine_profile(
        &cscl(),
        &observed,
        &PawleyProfileOptions {
            base: fixed_pawley(0.2),
            refine_eta: false,
            ..PawleyProfileOptions::default()
        },
    )
    .unwrap();
    assert_eq!(fixed.eta, 0.2);
    assert!(
        fixed.fit.rwp > result.fit.rwp * 100.0,
        "fixed Rwp {} refined {}",
        fixed.fit.rwp,
        result.fit.rwp
    );
}

#[test]
fn pawley_refines_the_width_when_asked() {
    let observed = synthetic_pawley(TRUE_ETA);
    let options = PawleyProfileOptions {
        base: fixed_pawley(TRUE_ETA),
        refine_eta: false,
        refine_width: true,
        ..PawleyProfileOptions::default()
    };
    let result = pawley_refine_profile(&cscl(), &observed, &options).unwrap();
    assert!(
        (result.fwhm - FWHM).abs() < 5e-3,
        "fwhm {} != {FWHM}",
        result.fwhm
    );
    assert!(result.fit.rwp < 1e-3, "Rwp {}", result.fit.rwp);
}

#[test]
fn le_bail_recovers_a_pseudo_voigt_mixing() {
    let observed = synthetic(TRUE_ETA);
    let options = LebailProfileOptions {
        base: LebailOptions {
            fwhm: FWHM,
            window: 3.0,
            iterations: 24,
            kind: ProfileKind::PseudoVoigt,
            eta: 0.2,
        },
        refine_eta: true,
        eta_range: [0.0, 1.0],
        ..LebailProfileOptions::default()
    };
    let result = le_bail_profile(&cscl(), &observed, &options).unwrap();
    assert!(
        (result.eta - TRUE_ETA).abs() < 0.05,
        "eta {} != {TRUE_ETA} with Rwp {}",
        result.eta,
        result.fit.rwp
    );
    assert!(
        result.fit.rwp < 0.12,
        "Rwp {} with Rp {}",
        result.fit.rwp,
        result.fit.r_factor
    );
    assert!(result.fit.r_factor < 0.03, "Rp {}", result.fit.r_factor);
    assert_eq!(result.fit.intensities.len(), result.fit.positions.len());
}

#[test]
fn rietveld_profile_recovers_a_pseudo_voigt_mixing() {
    let observed = synthetic(TRUE_ETA);
    let options = RietveldProfileOptions {
        base: RietveldOptions {
            radiation: Radiation::CuKa,
            two_theta_range: [15.0, 70.0],
            step: 0.02,
            fwhm: FWHM,
            ..RietveldOptions::default()
        },
        profile: InstrumentProfile {
            radiation: RadiationSource::from_radiation(Radiation::CuKa).unwrap(),
            shape: LineShape::PseudoVoigt,
            fwhm: FWHM,
            eta: 0.2,
            ..InstrumentProfile::default()
        },
        refine_eta: true,
        refine_width: false,
        eta_range: [0.0, 1.0],
        ..RietveldProfileOptions::default()
    };
    let result = fit_phase_profile(&cscl(), &observed, &options).unwrap();
    assert!(
        (result.eta - TRUE_ETA).abs() < 0.03,
        "eta {} != {TRUE_ETA} with Rwp {}",
        result.eta,
        result.phase.rwp
    );
    assert!(result.phase.rwp < 1e-3, "Rwp {}", result.phase.rwp);
    assert!(result.phase.scale.is_finite() && result.phase.scale > 0.0);
    assert!((result.width_scale - 1.0).abs() < 1e-12);
}

#[test]
fn profile_option_errors_are_typed() {
    let observed = synthetic(0.5);
    let bad_lebail = LebailProfileOptions {
        eta_range: [0.8, 0.2],
        ..LebailProfileOptions::default()
    };
    assert!(le_bail_profile(&cscl(), &observed, &bad_lebail).is_err());
    let bad_rietveld = RietveldProfileOptions {
        eta_range: [0.8, 0.2],
        ..RietveldProfileOptions::default()
    };
    assert!(fit_phase_profile(&cscl(), &observed, &bad_rietveld).is_err());
    assert!(fit_phase_profile(&cscl(), &[], &RietveldProfileOptions::default()).is_err());
}

#[test]
fn old_options_and_defaults_are_unchanged() {
    let defaults = PawleyOptions::default();
    assert_eq!(defaults.eta, 0.5);
    assert_eq!(defaults.fwhm, 0.1);
    assert!(defaults.refine_cell && defaults.refine_zero);
    assert!(!defaults.refine_size && !defaults.refine_strain);
    let profile_defaults = PawleyProfileOptions::default();
    assert_eq!(profile_defaults.base, defaults);
    assert!(profile_defaults.refine_eta);
    assert_eq!(profile_defaults.eta_range, [0.0, 1.0]);

    // An invalid eta window is refused before any fit runs.
    let invalid = PawleyProfileOptions {
        base: defaults,
        eta_range: [0.8, 0.2],
        ..PawleyProfileOptions::default()
    };
    assert!(pawley_refine_profile(&cscl(), &synthetic(0.5), &invalid).is_err());
}
