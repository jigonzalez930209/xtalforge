//! Integration tests for calculated powder diffraction (Task 5.1).

#![allow(
    clippy::unwrap_used,
    clippy::expect_used,
    clippy::panic,
    clippy::indexing_slicing
)]

use xtalforge_analysis::pxrd::{Pattern, PxrdOptions, Radiation, Reflection, calculate_pxrd};
use xtalforge_core::{Lattice, Site, Structure, StructureBuilder, SymOp};

const SI_A: f64 = 5.431;
const NACL_A: f64 = 5.6402;

fn site(label: &str, element: &str, frac: [f64; 3]) -> Site {
    Site::new(label, element, frac).unwrap()
}

/// Diamond cubic silicon in P1: the eight atoms of the conventional cell.
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
    let mut builder = StructureBuilder::new("Si", Lattice::cubic(SI_A).unwrap());
    for (index, frac) in positions.iter().enumerate() {
        builder = builder.site(site(&format!("Si{}", index + 1), "Si", *frac));
    }
    builder.build().unwrap()
}

/// Rock salt NaCl in P1: four Na and four Cl atoms.
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

fn reflection(pattern: &Pattern, hkl: [i32; 3]) -> Option<&Reflection> {
    pattern
        .reflections
        .iter()
        .find(|reflection| reflection.hkl == hkl)
}

fn bragg_two_theta(d: f64, wavelength: f64) -> f64 {
    2.0 * (wavelength / (2.0 * d)).asin().to_degrees()
}

#[test]
fn radiation_wavelengths_match_the_reference_values() {
    let expected = [
        (Radiation::CuKa1, 1.540_598),
        (Radiation::CuKa, 1.541_84),
        (Radiation::MoKa, 0.710_73),
        (Radiation::CoKa, 1.788_97),
        (Radiation::FeKa, 1.936_04),
        (Radiation::CrKa, 2.289_70),
        (Radiation::AgKa, 0.559_41),
    ];
    for (radiation, wavelength) in expected {
        assert_eq!(radiation.wavelength().unwrap(), wavelength);
        assert!(radiation.label().contains('A'));
    }
    assert_eq!(Radiation::Custom(0.9).wavelength().unwrap(), 0.9);
    for invalid in [0.0, -1.0, f64::NAN, f64::INFINITY] {
        assert!(Radiation::Custom(invalid).wavelength().is_err());
    }
}

#[test]
fn silicon_reflections_obey_bragg_and_are_sorted() {
    let options = PxrdOptions::default().with_radiation(Radiation::CuKa1);
    let pattern = calculate_pxrd(&silicon(), &options).unwrap();
    assert_eq!(pattern.reflections.len(), 74);
    let repeated = calculate_pxrd(&silicon(), &options).unwrap();
    assert_eq!(repeated.reflections.len(), pattern.reflections.len());
    let wavelength = Radiation::CuKa1.wavelength().unwrap();
    for reflection in &pattern.reflections {
        let expected = bragg_two_theta(reflection.d, wavelength);
        assert!((reflection.two_theta - expected).abs() < 1e-9);
        assert!(reflection.two_theta >= 5.0 && reflection.two_theta <= 80.0);
        assert!(reflection.intensity > 0.0);
        assert_eq!(reflection.multiplicity, 1);
    }
    for pair in pattern.reflections.windows(2) {
        assert!(pair[0].two_theta <= pair[1].two_theta + 1e-12);
    }
}

#[test]
fn silicon_111_matches_the_hand_calculation() {
    let options = PxrdOptions::default().with_radiation(Radiation::CuKa1);
    let pattern = calculate_pxrd(&silicon(), &options).unwrap();
    let reflection = reflection(&pattern, [1, 1, 1]).unwrap();
    let d_hand = SI_A / 3.0_f64.sqrt();
    let wavelength = 1.540_598;
    let expected = bragg_two_theta(d_hand, wavelength);
    assert!((reflection.d - d_hand).abs() < 1e-12);
    assert!(
        (reflection.two_theta - expected).abs() < 1e-3,
        "Si 111: {} vs {}",
        reflection.two_theta,
        expected
    );
    // Published Cu Ka1 position of the Si 111 peak.
    assert!((reflection.two_theta - 28.4429).abs() < 0.01);
}

#[test]
fn diamond_systematic_absences_are_absent() {
    let pattern = calculate_pxrd(&silicon(), &PxrdOptions::default()).unwrap();
    for absent in [[2, 0, 0], [0, 2, 0], [0, 0, 2], [2, 2, 2]] {
        assert!(
            reflection(&pattern, absent).is_none(),
            "diamond reflection {absent:?} should be systematically absent"
        );
    }
    assert!(reflection(&pattern, [1, 1, 1]).is_some());
    assert!(reflection(&pattern, [2, 2, 0]).is_some());
}

#[test]
fn halite_200_is_stronger_than_111_and_obeys_bragg() {
    let options = PxrdOptions::default().with_radiation(Radiation::CuKa);
    let pattern = calculate_pxrd(&halite(), &options).unwrap();
    assert_eq!(pattern.reflections.len(), 112);
    let repeated = calculate_pxrd(&halite(), &options).unwrap();
    assert_eq!(repeated.reflections.len(), pattern.reflections.len());
    let strong = reflection(&pattern, [2, 0, 0]).unwrap();
    let weak = reflection(&pattern, [1, 1, 1]).unwrap();
    assert!(strong.intensity > weak.intensity);
    let d_hand = NACL_A / 2.0;
    assert!((strong.d - d_hand).abs() < 1e-12);
    let expected = bragg_two_theta(d_hand, 1.541_84);
    assert!((strong.two_theta - expected).abs() < 1e-3);
    // (222) is the next allowed even reflection in the rock-salt lattice.
    let next = reflection(&pattern, [2, 2, 0]).unwrap();
    assert!(next.d < strong.d);
    assert!(next.two_theta > strong.two_theta);
}

#[test]
fn pattern_is_deterministic_and_sampled_on_the_step() {
    let options = PxrdOptions::default();
    let first = calculate_pxrd(&silicon(), &options).unwrap();
    let second = calculate_pxrd(&silicon(), &options).unwrap();
    assert_eq!(first, second);
    assert_eq!(first.points.len(), 3751);
    assert!((first.points[0].two_theta - 5.0).abs() < 1e-12);
    assert!((first.points[3750].two_theta - 80.0).abs() < 1e-9);
    for pair in first.points.windows(2) {
        assert!((pair[1].two_theta - pair[0].two_theta - 0.02).abs() < 1e-12);
        assert!(pair[0].intensity >= 0.0);
    }
    assert!(first.max_intensity() > 0.0);
    let normalised = first.normalised();
    assert!((normalised.iter().map(|p| p.intensity).fold(0.0, f64::max) - 100.0).abs() < 1e-9);
}

#[test]
fn symmetry_operations_give_multiplicities_and_merge_friedel_pairs() {
    let inversion = SymOp::new([[-1, 0, 0], [0, -1, 0], [0, 0, -1]], [0.0; 3]).unwrap();
    let merged_structure = StructureBuilder::new("P-1", Lattice::cubic(4.0).unwrap())
        .hall("P -1")
        .site(site("Si1", "Si", [0.1, 0.2, 0.3]))
        .operation(SymOp::identity())
        .operation(inversion)
        .build()
        .unwrap();
    let p1_structure = StructureBuilder::new("P1", Lattice::cubic(4.0).unwrap())
        .site(site("Si1", "Si", [0.1, 0.2, 0.3]))
        .site(site("Si2", "Si", [0.9, 0.8, 0.7]))
        .build()
        .unwrap();
    let options = PxrdOptions::default().with_radiation(Radiation::CuKa1);
    let merged = calculate_pxrd(&merged_structure, &options).unwrap();
    let plain = calculate_pxrd(&p1_structure, &options).unwrap();
    assert!(!merged.reflections.is_empty());
    for reflection in &merged.reflections {
        assert_eq!(reflection.multiplicity, 2);
    }
    assert_eq!(merged.reflections.len() * 2, plain.reflections.len());
}

#[test]
fn invalid_options_and_empty_structures_are_rejected() {
    let structure = silicon();
    let valid = PxrdOptions::default();
    let bad_options = [
        valid.with_step(0.0),
        valid.with_step(-0.02),
        valid.with_step(1e-300),
        valid.with_fwhm(0.0),
        valid.with_range(80.0, 5.0),
        valid.with_range(-1.0, 20.0),
        valid.with_range(5.0, 181.0),
        valid.with_radiation(Radiation::Custom(f64::NAN)),
    ];
    for options in bad_options {
        let error = calculate_pxrd(&structure, &options).unwrap_err();
        assert_eq!(error.code(), "invalid_input");
    }
    let empty = StructureBuilder::new("empty", Lattice::cubic(5.0).unwrap())
        .build()
        .unwrap();
    let error = calculate_pxrd(&empty, &valid).unwrap_err();
    assert_eq!(error.code(), "invalid_input");
}

#[test]
fn triclinic_cell_is_consistent_with_the_reciprocal_cell() {
    // A distorted cell exercises the general reciprocal length path.
    let lattice = Lattice::from_parameters(6.0, 7.0, 8.0, 80.0, 100.0, 110.0).unwrap();
    let structure = StructureBuilder::new("triclinic", lattice)
        .site(site("Si1", "Si", [0.13, 0.27, 0.41]))
        .build()
        .unwrap();
    let pattern = calculate_pxrd(
        &structure,
        &PxrdOptions::default().with_radiation(Radiation::MoKa),
    )
    .unwrap();
    assert!(!pattern.reflections.is_empty());
    for reflection in &pattern.reflections {
        let d = lattice.d_spacing(reflection.hkl.map(f64::from)).unwrap();
        assert!((reflection.d - d).abs() < 1e-12);
        let expected = bragg_two_theta(d, Radiation::MoKa.wavelength().unwrap());
        assert!((reflection.two_theta - expected).abs() < 1e-9);
    }
}
