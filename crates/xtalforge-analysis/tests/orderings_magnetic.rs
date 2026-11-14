//! Collinear magnetic ordering enumeration (Task 19.1.4).

#![allow(
    clippy::unwrap_used,
    clippy::expect_used,
    clippy::panic,
    clippy::indexing_slicing
)]

use serde_json::json;
use xtalforge_analysis::{
    MagneticOrderingKind, MagneticOrderingOptions, magnetic_orderings,
    magnetic_orderings_in_supercell,
};
use xtalforge_core::{Error, Lattice, MagneticAxis, Site, Structure, StructureBuilder};
use xtalforge_symmetry::MagneticSymmetryOptions;

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

fn ferrimagnet() -> Structure {
    StructureBuilder::new("ferrimagnet", Lattice::cubic(2.87).unwrap())
        .site(site("Fe1", [0.0, 0.0, 0.0], 5.0))
        .site(site("Fe2", [0.5, 0.5, 0.5], 3.0))
        .build()
        .unwrap()
}

fn options(
    supercell: [u32; 3],
    max_candidates: usize,
    max_orderings: usize,
) -> MagneticOrderingOptions {
    MagneticOrderingOptions {
        supercell,
        max_orderings,
        max_candidates,
        symmetry: MagneticSymmetryOptions {
            symprec: 1e-4,
            mag_symprec: Some(0.1),
            rotate_basis: false,
        },
    }
}

#[test]
fn bcc_iron_1x1x1_gives_ferro_then_antiferro() {
    let structure = bcc_iron(2.2);
    let settings = options([1, 1, 1], 8, 8);
    let orderings = magnetic_orderings(&structure, &settings).unwrap();
    assert_eq!(orderings.len(), 2, "{orderings:#?}");
    assert_eq!(orderings[0].kind, MagneticOrderingKind::Ferromagnetic);
    assert_eq!(orderings[1].kind, MagneticOrderingKind::Antiferromagnetic);
    assert_eq!(orderings[0].moments, vec![2.2, 2.2]);
    assert_eq!(orderings[1].moments, vec![2.2, -2.2]);
    assert!((orderings[0].net_moment - 4.4).abs() < 1e-12);
    assert!(orderings[1].net_moment.abs() < 1e-12);
    assert_eq!(orderings[0].operation_count, orderings[1].operation_count);
    assert_ne!(orderings[0].uni_number, orderings[1].uni_number);
    assert!(orderings[0].uni_number.is_some());
    assert!(orderings[1].uni_number.is_some());
    assert_eq!(orderings[0].magmom, "2.2000 2.2000");
    assert_eq!(orderings[1].magmom, "2.2000 -2.2000");
}

#[test]
fn enumeration_is_deterministic_and_sorted_by_symmetry() {
    let structure = simple_cubic(2.0);
    let settings = options([2, 2, 1], 8, 8);
    let first = magnetic_orderings(&structure, &settings).unwrap();
    let second = magnetic_orderings(&structure, &settings).unwrap();
    assert_eq!(first, second);
    assert_eq!(first.len(), 8, "expected every sign pattern");
    assert_eq!(first[0].kind, MagneticOrderingKind::Ferromagnetic);
    for pair in first.windows(2) {
        assert!(
            pair[0].operation_count >= pair[1].operation_count,
            "{} < {}",
            pair[0].operation_count,
            pair[1].operation_count
        );
    }
    assert!(first.iter().any(
        |ordering| ordering.kind == MagneticOrderingKind::Antiferromagnetic
            && ordering.net_moment.abs() < 1e-12
    ));
    assert!(
        first
            .iter()
            .any(|ordering| ordering.kind == MagneticOrderingKind::Ferrimagnetic)
    );
    for ordering in &first {
        assert_eq!(ordering.moments.len(), 4);
        assert_eq!(ordering.structure.sites().len(), 4);
        assert_eq!(
            ordering.structure.p1().unwrap().atoms.len(),
            ordering.moments.len()
        );
    }
}

#[test]
fn ferrimagnet_has_ferro_and_ferri_orderings() {
    let structure = ferrimagnet();
    let orderings = magnetic_orderings(&structure, &options([1, 1, 1], 8, 8)).unwrap();
    assert_eq!(orderings.len(), 2);
    let ferro = orderings
        .iter()
        .find(|ordering| ordering.kind == MagneticOrderingKind::Ferromagnetic)
        .unwrap();
    let ferri = orderings
        .iter()
        .find(|ordering| ordering.kind == MagneticOrderingKind::Ferrimagnetic)
        .unwrap();
    assert!((ferro.net_moment - 8.0).abs() < 1e-12);
    assert!((ferri.net_moment - 2.0).abs() < 1e-12);
    assert_eq!(ferro.magmom, "5.0000 3.0000");
    assert_eq!(ferri.magmom, "5.0000 -3.0000");
    assert!(ferro.uni_number.is_some());
    assert!(ferri.uni_number.is_some());
}

#[test]
fn supplied_supercell_matches_the_expanded_cell() {
    let structure = simple_cubic(2.0);
    let supercell = structure.supercell(2, 1, 1).unwrap();
    let settings = options([2, 1, 1], 16, 16);
    let expanded = magnetic_orderings(&structure, &settings).unwrap();
    let supplied = magnetic_orderings_in_supercell(&supercell, &settings).unwrap();
    assert_eq!(expanded, supplied);
    assert_eq!(supplied.len(), 2);
    assert!(
        supplied
            .iter()
            .all(|ordering| ordering.structure.sites().len() == 2)
    );
}

#[test]
fn max_orderings_truncates_the_deterministic_order() {
    let structure = simple_cubic(2.0);
    let all = magnetic_orderings(&structure, &options([2, 2, 1], 8, 8)).unwrap();
    let truncated = magnetic_orderings(&structure, &options([2, 2, 1], 8, 5)).unwrap();
    assert_eq!(truncated.len(), 5);
    assert_eq!(truncated, all[..5]);
}

#[test]
fn typed_collinear_moments_enumerate_like_properties() {
    let typed_site = |label: &str, frac: [f64; 3], moment: f64| {
        Site::new(label, "Fe", frac)
            .unwrap()
            .with_collinear_moment(MagneticAxis::C, moment)
            .unwrap()
    };
    let structure = StructureBuilder::new("typed-bcc", Lattice::cubic(2.87).unwrap())
        .site(typed_site("Fe1", [0.0, 0.0, 0.0], 2.2))
        .site(typed_site("Fe2", [0.5, 0.5, 0.5], 2.2))
        .build()
        .unwrap();
    let orderings = magnetic_orderings(&structure, &options([1, 1, 1], 8, 8)).unwrap();
    assert_eq!(orderings.len(), 2);
    assert_eq!(orderings[0].moments, vec![2.2, 2.2]);
    assert_eq!(orderings[1].moments, vec![2.2, -2.2]);
    assert!(
        orderings
            .iter()
            .all(|ordering| ordering.uni_number.is_some())
    );
}

#[test]
fn bad_options_and_inputs_are_rejected() {
    let structure = simple_cubic(2.0);
    assert!(matches!(
        magnetic_orderings(&structure, &options([0, 1, 1], 8, 8)).unwrap_err(),
        Error::InvalidInput(_)
    ));
    let mut zero_orderings = options([1, 1, 1], 8, 0);
    zero_orderings.max_orderings = 0;
    assert!(matches!(
        magnetic_orderings(&structure, &zero_orderings).unwrap_err(),
        Error::InvalidInput(_)
    ));
    let mut zero_candidates = options([1, 1, 1], 0, 8);
    zero_candidates.max_candidates = 0;
    assert!(matches!(
        magnetic_orderings(&structure, &zero_candidates).unwrap_err(),
        Error::InvalidInput(_)
    ));
    let nonmagnetic = StructureBuilder::new("no-moments", Lattice::cubic(2.87).unwrap())
        .site(Site::new("Fe1", "Fe", [0.0, 0.0, 0.0]).unwrap())
        .build()
        .unwrap();
    assert!(matches!(
        magnetic_orderings(&nonmagnetic, &options([1, 1, 1], 8, 8)).unwrap_err(),
        Error::InvalidInput(_)
    ));
    let noncollinear = StructureBuilder::new("vector", Lattice::cubic(2.87).unwrap())
        .site(
            Site::new("Fe1", "Fe", [0.0, 0.0, 0.0])
                .unwrap()
                .with_prop("magmom", json!([0.0, 0.0, 2.0])),
        )
        .build()
        .unwrap();
    assert!(matches!(
        magnetic_orderings(&noncollinear, &options([1, 1, 1], 8, 8)).unwrap_err(),
        Error::InvalidInput(_)
    ));
}
