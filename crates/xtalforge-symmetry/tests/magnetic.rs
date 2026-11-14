//! Magnetic symmetry detection and moment symmetrisation (Task 19.1.2).
//!
//! Golden data from spglib's magnetic tests is not redistributed with this
//! repository, so the cases are hand-built and deterministic. The deviation is
//! documented in `crates/xtalforge-symmetry/src/magnetic.rs`.

#![allow(
    clippy::unwrap_used,
    clippy::expect_used,
    clippy::panic,
    clippy::indexing_slicing
)]

use serde_json::json;
use xtalforge_core::{Error, Lattice, MagneticAxis, Site, Structure, StructureBuilder};
use xtalforge_symmetry::{
    MagneticSymmetryOptions, detect_magnetic_symmetry, detect_magnetic_symmetry_with_options,
    format_magmom, magmom_values, symmetrize_moments, symmetrize_structure_moments,
};

fn site(label: &str, frac: [f64; 3], magmom: f64) -> Site {
    Site::new(label, "Fe", frac)
        .unwrap()
        .with_prop("magmom", json!(magmom))
}

fn bcc_iron(moment: f64) -> Structure {
    StructureBuilder::new("bcc-Fe", Lattice::cubic(2.87).unwrap())
        .site(site("Fe1", [0.0, 0.0, 0.0], moment))
        .site(site("Fe2", [0.5, 0.5, 0.5], moment))
        .build()
        .unwrap()
}

fn simple_cubic(moment: f64) -> Structure {
    StructureBuilder::new("sc-Fe", Lattice::cubic(2.87).unwrap())
        .site(site("Fe1", [0.0, 0.0, 0.0], moment))
        .build()
        .unwrap()
}

fn checkerboard() -> Structure {
    let supercell = simple_cubic(2.0).supercell(2, 2, 2).unwrap();
    let sites = supercell
        .sites()
        .iter()
        .map(|site| {
            let frac = site.frac();
            let parity =
                (2.0 * frac[0]).round() + (2.0 * frac[1]).round() + (2.0 * frac[2]).round();
            let sign = if parity % 2.0 == 0.0 { 1.0 } else { -1.0 };
            site.clone().with_prop("magmom", json!(2.0 * sign))
        })
        .collect::<Vec<_>>();
    let mut builder = StructureBuilder::new("sc-AFM", *supercell.lattice());
    for site in sites {
        builder = builder.site(site);
    }
    builder.build().unwrap()
}

fn ferrimagnet() -> Structure {
    let mut first = site("Fe1", [0.0, 0.0, 0.0], 5.0);
    first = first.with_prop("magmom", json!(5.0));
    let mut second = site("Fe2", [0.5, 0.5, 0.5], -3.0);
    second = second.with_prop("magmom", json!(-3.0));
    StructureBuilder::new("ferrimagnet", Lattice::cubic(2.87).unwrap())
        .site(first)
        .site(second)
        .build()
        .unwrap()
}

#[test]
fn ferromagnetic_bcc_iron_has_a_nonzero_magnetic_group() {
    let dataset = detect_magnetic_symmetry(&bcc_iron(2.2), 1e-3).unwrap();
    assert!(dataset.uni_number > 0);
    assert!(dataset.operations.len() > 1);
    assert_eq!(dataset.orbits, vec![0, 0]);
    assert!(
        dataset
            .operations
            .iter()
            .any(|operation| !operation.time_reversal)
    );
}

#[test]
fn ferromagnetic_bcc_iron_symmetrises_noisy_moments() {
    let structure = bcc_iron(2.2);
    let options = MagneticSymmetryOptions {
        mag_symprec: Some(0.3),
        ..MagneticSymmetryOptions::default()
    };
    let result = symmetrize_moments(&structure, &[2.2, 2.0], &options).unwrap();
    assert_eq!(result.moments.len(), 2);
    assert!(
        (result.moments[0] - 2.1).abs() < 1e-9,
        "{:?}",
        result.moments
    );
    assert!(
        (result.moments[1] - 2.1).abs() < 1e-9,
        "{:?}",
        result.moments
    );
    assert_eq!(result.orbits, vec![0, 0]);
    assert_eq!(result.stabilizers.len(), 2);
    assert!(result.stabilizers.iter().all(|count| *count > 0));
    assert!(result.operation_count >= 2);
    assert!(result.time_reversal_count <= result.operation_count);
    assert!(result.uni_number > 0);
}

#[test]
fn checkerboard_simple_cubic_is_antiferromagnetic_and_stable() {
    let structure = checkerboard();
    let moments = magmom_values(&structure).unwrap();
    assert_eq!(moments.len(), 8);
    let options = MagneticSymmetryOptions {
        symprec: 1e-4,
        mag_symprec: Some(0.1),
        rotate_basis: false,
    };
    let result = symmetrize_moments(&structure, &moments, &options).unwrap();
    for (before, after) in moments.iter().zip(&result.moments) {
        assert!((before - after).abs() < 1e-9, "{before} != {after}");
    }
    assert!(result.moments.iter().sum::<f64>().abs() < 1e-9);
    assert!(result.operation_count > 1);
    let orbit_magnitudes = result
        .orbits
        .iter()
        .map(|representative| result.moments[*representative].abs())
        .collect::<Vec<_>>();
    assert!(
        orbit_magnitudes
            .iter()
            .all(|value| (*value - 2.0).abs() < 1e-9)
    );
}

#[test]
fn two_sublattice_ferrimagnet_keeps_distinct_magnitudes() {
    let structure = ferrimagnet();
    let options = MagneticSymmetryOptions::default();
    let result = symmetrize_structure_moments(&structure, &options).unwrap();
    assert!(
        (result.moments[0] - 5.0).abs() < 1e-9,
        "{:?}",
        result.moments
    );
    assert!(
        (result.moments[1] + 3.0).abs() < 1e-9,
        "{:?}",
        result.moments
    );
    assert_eq!(result.stabilizers, vec![48, 48]);
    assert_eq!(result.operation_count, 48);
}

#[test]
fn detection_with_options_accepts_the_same_cell() {
    let structure = bcc_iron(2.2);
    let options = MagneticSymmetryOptions::default();
    let strict = detect_magnetic_symmetry_with_options(&structure, &options).unwrap();
    let legacy = detect_magnetic_symmetry(&structure, 1e-3).unwrap();
    assert!(strict.uni_number > 0);
    assert!(legacy.uni_number > 0);
    assert_eq!(strict.orbits, legacy.orbits);
}

#[test]
fn non_collinear_moments_are_rejected() {
    let structure = StructureBuilder::new("vector", Lattice::cubic(2.87).unwrap())
        .site(
            Site::new("Fe1", "Fe", [0.0, 0.0, 0.0])
                .unwrap()
                .with_prop("magmom", json!([0.0, 0.0, 2.0])),
        )
        .build()
        .unwrap();
    let error = magmom_values(&structure).unwrap_err();
    assert!(matches!(error, Error::InvalidInput(_)));
    let error =
        symmetrize_structure_moments(&structure, &MagneticSymmetryOptions::default()).unwrap_err();
    assert!(matches!(error, Error::InvalidInput(_)));
}

#[test]
fn moment_array_shape_must_match_the_p1_cell() {
    let structure = bcc_iron(2.2);
    let error =
        symmetrize_moments(&structure, &[2.2], &MagneticSymmetryOptions::default()).unwrap_err();
    match error {
        Error::InvalidInput(message) => {
            assert!(message.contains("expected 2 moments"), "{message}")
        }
        other => panic!("unexpected error: {other:?}"),
    }
}

#[test]
fn invalid_tolerances_are_rejected() {
    let structure = bcc_iron(2.2);
    let options = MagneticSymmetryOptions {
        symprec: 0.0,
        ..MagneticSymmetryOptions::default()
    };
    assert!(matches!(
        symmetrize_moments(&structure, &[2.2, 2.2], &options).unwrap_err(),
        Error::Symmetry(_)
    ));
    let options = MagneticSymmetryOptions {
        mag_symprec: Some(f64::NAN),
        ..MagneticSymmetryOptions::default()
    };
    assert!(matches!(
        symmetrize_moments(&structure, &[2.2, 2.2], &options).unwrap_err(),
        Error::Symmetry(_)
    ));
}

#[test]
fn magmom_values_and_format_match_the_vasp_writer() {
    let values = magmom_values(&checkerboard()).unwrap();
    assert_eq!(values.len(), 8);
    let text = format_magmom(&[2.0, -2.0, 0.0]).unwrap();
    assert_eq!(text, "2.0000 -2.0000 0.0000");
    assert!(format_magmom(&[]).is_err());
    assert!(format_magmom(&[f64::INFINITY]).is_err());
}

#[test]
fn typed_collinear_site_moments_are_supported() {
    let first = Site::new("Fe1", "Fe", [0.0, 0.0, 0.0])
        .unwrap()
        .with_collinear_moment(MagneticAxis::C, 2.2)
        .unwrap();
    let second = Site::new("Fe2", "Fe", [0.5, 0.5, 0.5])
        .unwrap()
        .with_collinear_moment(MagneticAxis::C, -2.2)
        .unwrap();
    let structure = StructureBuilder::new("typed", Lattice::cubic(2.87).unwrap())
        .site(first)
        .site(second)
        .build()
        .unwrap();
    assert_eq!(magmom_values(&structure).unwrap(), vec![2.2, -2.2]);
    let options = MagneticSymmetryOptions {
        symprec: 1e-4,
        mag_symprec: Some(0.1),
        rotate_basis: false,
    };
    let result = symmetrize_structure_moments(&structure, &options).unwrap();
    assert!((result.moments[0] - 2.2).abs() < 1e-9);
    assert!((result.moments[1] + 2.2).abs() < 1e-9);
}

#[test]
fn standardization_keeps_an_unambiguous_site_moment() {
    let structure = StructureBuilder::new("fe", Lattice::cubic(2.87).unwrap())
        .site(
            Site::new("Fe1", "Fe", [0.0, 0.0, 0.0])
                .unwrap()
                .with_collinear_moment(MagneticAxis::C, 2.2)
                .unwrap(),
        )
        .build()
        .unwrap();
    let result =
        xtalforge_symmetry::standardize(&structure, xtalforge_symmetry::SymmetryOptions::default())
            .unwrap();
    assert!(!result.structure.sites().is_empty());
    for site in result.structure.sites() {
        assert!(
            site.magnetic_moment().is_some(),
            "moment lost at {}",
            site.label()
        );
    }
}
