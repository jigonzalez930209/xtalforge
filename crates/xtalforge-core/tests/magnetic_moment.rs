//! Integration tests for magnetic moments as site data (Task 19.1.1).

#![allow(
    clippy::unwrap_used,
    clippy::expect_used,
    clippy::panic,
    clippy::indexing_slicing
)]

use serde_json::json;
use xtalforge_core::{
    Lattice, MagneticAxis, MagneticMoment, Site, Structure, StructureBuilder, StructureDocument,
    SymOp, make_slab,
};

fn magnetic_site() -> Site {
    Site::new("Mn1", "Mn", [0.25, 0.25, 0.25])
        .unwrap()
        .with_occupancy(0.75)
        .unwrap()
        .with_u_iso(0.01)
        .unwrap()
        .with_charge(2.0)
        .unwrap()
        .with_collinear_moment(MagneticAxis::C, 4.2)
        .unwrap()
        .with_prop("source", json!("test"))
}

fn magnetic_structure() -> Structure {
    StructureBuilder::new("Mn", Lattice::cubic(5.0).unwrap())
        .site(magnetic_site())
        .build()
        .unwrap()
}

#[test]
fn assignment_accessors_and_removal() {
    let plain = Site::new("C1", "C", [0.0, 0.0, 0.0]).unwrap();
    assert!(plain.magnetic_moment().is_none());
    assert!(!plain.is_partial());

    let site = plain.with_collinear_moment(MagneticAxis::B, -3.5).unwrap();
    let moment = site.magnetic_moment().unwrap();
    assert_eq!(moment.components(), [0.0, -3.5, 0.0]);
    assert!(moment.is_collinear());
    assert!(!moment.is_zero());
    assert!((moment.magnitude() - 3.5).abs() < 1e-12);

    let site = site.with_moment_vector([1.0, -2.0, 0.5]).unwrap();
    let moment = site.magnetic_moment().unwrap();
    assert!(!moment.is_collinear());
    assert_eq!(moment.components(), [1.0, -2.0, 0.5]);
    assert!((moment.magnitude() - 5.25_f64.sqrt()).abs() < 1e-12);

    let site = site.with_magnetic_moment(MagneticMoment::vector([0.0; 3]).unwrap());
    assert!(site.magnetic_moment().unwrap().is_zero());
    assert!(site.without_magnetic_moment().magnetic_moment().is_none());

    assert_eq!(
        MagneticMoment::vector([0.0, 0.0, 4.0]).unwrap(),
        MagneticMoment::collinear(MagneticAxis::C, 4.0).unwrap()
    );
    let zero = MagneticMoment::collinear(MagneticAxis::A, 0.0).unwrap();
    assert!(zero.is_zero());
    assert!(!zero.is_collinear());
}

#[test]
fn invalid_moments_and_labels_are_rejected() {
    assert!(MagneticMoment::collinear(MagneticAxis::A, f64::NAN).is_err());
    assert!(MagneticMoment::collinear(MagneticAxis::A, f64::INFINITY).is_err());
    assert!(MagneticMoment::vector([0.0, f64::INFINITY, 0.0]).is_err());
    let site = Site::new("C1", "C", [0.0, 0.0, 0.0]).unwrap();
    assert!(
        site.clone()
            .with_collinear_moment(MagneticAxis::A, f64::NAN)
            .is_err()
    );
    assert!(site.with_moment_vector([f64::NAN, 0.0, 0.0]).is_err());
    assert!(
        Site::new("C1", "C", [0.0, 0.0, 0.0])
            .unwrap()
            .with_label("   ")
            .is_err()
    );
}

#[test]
fn clone_and_builder_keep_the_moment() {
    let structure = magnetic_structure();
    let expected = Some(MagneticMoment::collinear(MagneticAxis::C, 4.2).unwrap());
    assert_eq!(structure.sites()[0].magnetic_moment(), expected);
    assert_eq!(structure.clone().sites()[0].magnetic_moment(), expected);
    assert_eq!(structure.sites()[0].clone().magnetic_moment(), expected);
}

#[test]
fn p1_expansion_preserves_moments_and_site_properties() {
    let inversion = SymOp::parse("-x,-y,-z").unwrap();
    let structure = StructureBuilder::new("Mn", Lattice::cubic(5.0).unwrap())
        .hall("P -1")
        .operation(SymOp::identity())
        .operation(inversion)
        .site(
            Site::new("Mn1", "Mn", [0.1, 0.2, 0.3])
                .unwrap()
                .with_occupancy(0.5)
                .unwrap()
                .with_charge(3.0)
                .unwrap()
                .with_prop("magmom", json!(2.0))
                .with_moment_vector([1.5, -0.5, 2.5])
                .unwrap(),
        )
        .build()
        .unwrap();
    let expanded = structure.as_p1().unwrap();
    assert_eq!(expanded.sites().len(), 2);
    for site in expanded.sites() {
        assert_eq!(
            site.magnetic_moment().unwrap().components(),
            [1.5, -0.5, 2.5]
        );
        assert!((site.occupancy() - 0.5).abs() < 1e-12);
        assert_eq!(site.charge(), Some(3.0));
        assert_eq!(site.props()["magmom"], json!(2.0));
    }
}

#[test]
fn supercell_preserves_moments_and_site_properties() {
    let structure = magnetic_structure();
    let cell = structure.supercell(2, 1, 2).unwrap();
    assert_eq!(cell.sites().len(), 4);
    for site in cell.sites() {
        assert_eq!(
            site.magnetic_moment(),
            Some(MagneticMoment::collinear(MagneticAxis::C, 4.2).unwrap())
        );
        assert!((site.occupancy() - 0.75).abs() < 1e-12);
        assert_eq!(site.u_iso(), Some(0.01));
        assert_eq!(site.charge(), Some(2.0));
        assert_eq!(site.props()["source"], json!("test"));
    }
}

#[test]
fn slab_preserves_moments_and_site_properties() {
    let structure = magnetic_structure();
    let slab = make_slab(&structure, [0, 0, 1], 1, 8.0, 0.0).unwrap();
    assert!(!slab.structure.sites().is_empty());
    for site in slab.structure.sites() {
        assert_eq!(
            site.magnetic_moment(),
            structure.sites()[0].magnetic_moment()
        );
        assert_eq!(site.occupancy(), structure.sites()[0].occupancy());
        assert_eq!(site.u_iso(), structure.sites()[0].u_iso());
        assert_eq!(site.charge(), structure.sites()[0].charge());
        assert_eq!(site.props(), structure.sites()[0].props());
    }
}

#[test]
fn document_round_trip_keeps_moments() {
    let structure = magnetic_structure();
    let document = structure.to_document();
    assert_eq!(document.schema_version, xtalforge_core::SCHEMA_VERSION);
    let text = document.to_json().unwrap();
    assert!(text.contains("\"magnetic_moment\""));
    assert!(text.contains("\"collinear\""));
    let back = StructureDocument::from_json(&text).unwrap();
    assert_eq!(back, document);
    assert_eq!(
        back.structure.sites()[0].magnetic_moment(),
        Some(MagneticMoment::collinear(MagneticAxis::C, 4.2).unwrap())
    );

    let plain = StructureBuilder::new("C", Lattice::cubic(4.0).unwrap())
        .site(Site::new("C1", "C", [0.0, 0.0, 0.0]).unwrap())
        .build()
        .unwrap();
    let text = plain.to_document().to_json().unwrap();
    assert!(!text.contains("magnetic_moment"));
}
