//! Integration tests for the McCusker refinement stages (Task 7.8.1).
//!
//! Every synthetic pattern is generated with the same forward model the
//! refinement uses (a Gaussian profile, a constant background and an optional
//! March-Dollase factor), so the correct answer is known exactly.

#![allow(
    clippy::unwrap_used,
    clippy::expect_used,
    clippy::panic,
    clippy::indexing_slicing
)]

use xtalforge_analysis::{
    Background, ProfileKind, PxrdOptions, Radiation, RefineOptions, RefinedPhase, RietveldOptions,
    calculate_pxrd, free_directions, refine_phase,
};
use xtalforge_core::{Lattice, Site, Structure, StructureBuilder};
use xtalforge_io::xy::PatternPoint;

const BASE_FWHM: f64 = 0.1;

fn cubic(a: f64) -> Lattice {
    Lattice::from_parameters(a, a, a, 90.0, 90.0, 90.0).unwrap()
}

fn cscl(a: f64, u: Option<f64>, cesium: f64, chlorine: f64) -> Structure {
    let site = |label: &str, element: &str, frac: [f64; 3], occupancy: f64| {
        let mut site = Site::new(label, element, frac)
            .unwrap()
            .with_occupancy(occupancy)
            .unwrap();
        if let Some(u) = u {
            site = site.with_u_iso(u).unwrap();
        }
        site
    };
    StructureBuilder::new("CsCl", cubic(a))
        .site(site("Cs1", "Cs", [0.0, 0.0, 0.0], cesium))
        .site(site("Cl1", "Cl", [0.5, 0.5, 0.5], chlorine))
        .build()
        .unwrap()
}

fn rietveld(fwhm: f64) -> RietveldOptions {
    RietveldOptions {
        radiation: Radiation::CuKa1,
        two_theta_range: [10.0, 80.0],
        step: 0.02,
        fwhm,
        background: Background::Constant(0.0),
        kind: ProfileKind::Gaussian,
    }
}

/// March-Dollase factor for a cubic cell and the `c` axis.
fn march_dollase(hkl: [i32; 3], ratio: f64) -> f64 {
    let [h, k, l] = hkl.map(f64::from);
    let norm = (h * h + k * k + l * l).sqrt();
    let cosine = l / norm;
    let squared = cosine * cosine;
    let sine = 1.0 - squared;
    (ratio * ratio * squared + sine / ratio).powf(-1.5)
}

/// Samples a synthetic pattern at the observed grid.
fn synthetic(
    structure: &Structure,
    zero: f64,
    fwhm: f64,
    texture: Option<f64>,
) -> Vec<PatternPoint> {
    let pattern = calculate_pxrd(
        structure,
        &PxrdOptions {
            radiation: Radiation::CuKa1,
            two_theta_range: [10.0, 80.0],
            step: 0.02,
            fwhm,
        },
    )
    .unwrap();
    let normalisation = 2.0 / fwhm * (std::f64::consts::LN_2 / std::f64::consts::PI).sqrt();
    let count = ((80.0_f64 - 10.0) / 0.02).floor() as usize + 1;
    let mut points = Vec::with_capacity(count);
    for index in 0..count {
        let two_theta = 10.0 + index as f64 * 0.02;
        let mut intensity = 5.0;
        for reflection in &pattern.reflections {
            let factor = texture.map_or(1.0, |ratio| march_dollase(reflection.hkl, ratio));
            let delta = (two_theta - (reflection.two_theta + zero)) / fwhm;
            if delta.abs() <= 8.0 {
                intensity += reflection.intensity
                    * factor
                    * normalisation
                    * (-4.0 * std::f64::consts::LN_2 * delta * delta).exp();
            }
        }
        points.push(PatternPoint {
            two_theta,
            intensity,
            uncertainty: None,
        });
    }
    points
}

fn parameter(result: &RefinedPhase, name: &str) -> f64 {
    result
        .parameters
        .iter()
        .find(|(key, _)| key == name)
        .map(|(_, value)| *value)
        .unwrap_or_else(|| panic!("missing parameter {name}"))
}

fn esd(result: &RefinedPhase, name: &str) -> f64 {
    result
        .parameter_esds
        .iter()
        .find(|(key, _)| key == name)
        .map(|(_, value)| *value)
        .unwrap_or_else(|| panic!("missing esd {name}"))
}

#[test]
fn cell_zero_and_width_recover_a_synthetic_pattern() {
    let truth = cscl(4.0, None, 1.0, 1.0);
    let observed = synthetic(&truth, 0.03, 0.08, None);
    let options = RefineOptions {
        refine_zero: true,
        refine_cell: true,
        refine_width: true,
        max_iterations: 50,
        ..RefineOptions::default()
    };
    let result = refine_phase(&truth, &observed, &rietveld(BASE_FWHM), &options).unwrap();
    assert_eq!(
        result.stages,
        vec!["background_scale", "zero_cell", "width"]
    );
    assert!(result.rwp < 1e-4, "Rwp {}", result.rwp);
    assert!((parameter(&result, "zero") - 0.03).abs() < 2e-3);
    assert!((parameter(&result, "cell.0") - 4.0).abs() < 2e-3);
    assert!(
        (parameter(&result, "width_scale") - 0.8).abs() < 5e-3,
        "zero {} cell {} width {} rwp {}",
        parameter(&result, "zero"),
        parameter(&result, "cell.0"),
        parameter(&result, "width_scale"),
        result.rwp
    );
    assert!(result.rwp < 1e-4, "Rwp {}", result.rwp);
    for name in ["zero", "cell.0", "width_scale"] {
        assert!(esd(&result, name) > 0.0, "{name} has no esd");
    }
}

#[test]
fn texture_recovers_the_march_dollase_ratio() {
    let truth = cscl(4.0, None, 1.0, 1.0);
    let observed = synthetic(&truth, 0.0, BASE_FWHM, Some(1.5));
    let options = RefineOptions {
        refine_texture: true,
        texture_axis: [0.0, 0.0, 1.0],
        max_iterations: 50,
        ..RefineOptions::default()
    };
    let result = refine_phase(&truth, &observed, &rietveld(BASE_FWHM), &options).unwrap();
    assert_eq!(result.stages, vec!["background_scale", "structure"]);
    assert!(
        (parameter(&result, "texture") - 1.5).abs() < 2e-3,
        "texture {} rwp {}",
        parameter(&result, "texture"),
        result.rwp
    );
    assert!(result.rwp < 1e-5, "Rwp {}", result.rwp);
    assert!(esd(&result, "texture") > 0.0);
}

#[test]
fn adp_scale_recovers_a_known_thermal_multiplier() {
    let truth = cscl(4.0, Some(0.02), 1.0, 1.0);
    let start = cscl(4.0, Some(0.005), 1.0, 1.0);
    let observed = synthetic(&truth, 0.0, BASE_FWHM, None);
    let options = RefineOptions {
        refine_adps: true,
        max_iterations: 20,
        ..RefineOptions::default()
    };
    let result = refine_phase(&start, &observed, &rietveld(BASE_FWHM), &options).unwrap();
    assert!(
        (parameter(&result, "adp_scale") - 4.0).abs() < 0.02,
        "adp_scale {}",
        parameter(&result, "adp_scale")
    );
    assert!(result.rwp < 1e-5, "Rwp {}", result.rwp);
    for site in result.structure.sites() {
        assert!((site.u_iso().unwrap_or(0.0) - 0.02).abs() < 2e-4);
    }
    assert!(esd(&result, "adp_scale") > 0.0);
}

#[test]
fn occupancy_ratios_are_identifiable_under_the_scale() {
    // The scale multiplies every reflection, so only the occupancy ratios are
    // observable. The test checks the ratio and the bounds, not the absolute
    // values that the scale absorbs.
    let truth = cscl(4.0, Some(0.02), 0.9, 0.5);
    let start = cscl(4.0, Some(0.005), 0.7, 0.5);
    let observed = synthetic(&truth, 0.0, BASE_FWHM, None);
    let options = RefineOptions {
        refine_adps: true,
        refine_occupancies: true,
        max_iterations: 30,
        ..RefineOptions::default()
    };
    let result = refine_phase(&start, &observed, &rietveld(BASE_FWHM), &options).unwrap();
    let cesium = result.structure.sites()[0].occupancy();
    let chlorine = result.structure.sites()[1].occupancy();
    assert!(
        (cesium / chlorine - 1.8).abs() < 0.05,
        "occupancy ratio {}",
        cesium / chlorine
    );
    for occupancy in [cesium, chlorine] {
        assert!((0.0..=1.0).contains(&occupancy));
        assert!(occupancy > 0.0);
    }
    assert!(result.rwp < 1e-3, "Rwp {}", result.rwp);
    assert!(esd(&result, "Cs1.occupancy") > 0.0);
}

#[test]
fn every_stage_runs_in_mccusker_order() {
    let rutile = include_str!("../../../tests/corpus/simple/TiO2.cif");
    let base = xtalforge_io::cif::read_cif(rutile).unwrap();
    let mut builder = StructureBuilder::new(
        "rutile",
        Lattice::from_parameters(4.5937, 4.5937, 2.9587, 90.0, 90.0, 90.0).unwrap(),
    )
    .hall(base.hall());
    builder = builder.site(Site::new("Ti", "Ti", [0.0, 0.0, 0.0]).unwrap());
    builder = builder.site(Site::new("O", "O", [0.30478, 0.30478, 0.0]).unwrap());
    for operation in base.operations() {
        builder = builder.operation(*operation);
    }
    let structure = builder.build().unwrap();
    let observed = synthetic(&structure, 0.03, BASE_FWHM, Some(1.5));
    let options = RefineOptions {
        refine_positions: true,
        refine_adps: true,
        refine_occupancies: true,
        refine_zero: true,
        refine_cell: true,
        refine_width: true,
        refine_texture: true,
        max_iterations: 10,
        ..RefineOptions::default()
    };
    let result = refine_phase(&structure, &observed, &rietveld(BASE_FWHM), &options).unwrap();
    assert_eq!(
        result.stages,
        vec!["background_scale", "zero_cell", "width", "structure"]
    );
    for name in [
        "zero",
        "cell.0",
        "cell.2",
        "width_scale",
        "adp_scale",
        "texture",
        "Ti.occupancy",
        "O.occupancy",
        "O.d0",
    ] {
        assert!(
            result.parameters.iter().any(|(key, _)| key == name),
            "missing {name}"
        );
    }
}

#[test]
fn wyckoff_directions_follow_the_site_symmetry() {
    let rutile = include_str!("../../../tests/corpus/simple/TiO2.cif");
    let base = xtalforge_io::cif::read_cif(rutile).unwrap();
    assert!(free_directions(base.operations(), [0.0, 0.0, 0.0]).is_empty());
    let oxygen = free_directions(base.operations(), [0.30478, 0.30478, 0.0]);
    assert_eq!(oxygen.len(), 1);
    assert!((oxygen[0][0] - oxygen[0][1]).abs() < 1e-12);
    assert!(oxygen[0][2].abs() < 1e-12);
    let general = free_directions(base.operations(), [0.1, 0.2, 0.3]);
    assert_eq!(general.len(), 3);
    let p1 = StructureBuilder::new("p1", cubic(4.0))
        .site(Site::new("C1", "C", [0.0, 0.0, 0.0]).unwrap())
        .build()
        .unwrap();
    assert_eq!(free_directions(p1.operations(), [0.1, 0.2, 0.3]).len(), 3);
}
