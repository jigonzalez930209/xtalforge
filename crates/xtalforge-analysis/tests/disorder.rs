//! Integration tests for partial occupancies and the scattering-weight hook.

#![allow(
    clippy::unwrap_used,
    clippy::expect_used,
    clippy::panic,
    clippy::indexing_slicing
)]

use xtalforge_analysis::{
    RadiationKind, disorder_composition, disorder_density, disorder_formula, disorder_groups,
    disorder_scattering_weights, neutron_scattering_length, xray_form_factor,
};
use xtalforge_core::{Lattice, Site, Structure, StructureBuilder};

fn zinc_oxide_alloy() -> Structure {
    StructureBuilder::new("ZnO alloy", Lattice::cubic(4.0).unwrap())
        .site(
            Site::new("Zn1", "Zn", [0.0, 0.0, 0.0])
                .unwrap()
                .with_occupancy(0.5)
                .unwrap()
                .with_disorder_group(7),
        )
        .site(
            Site::new("O1", "O", [0.0, 0.0, 0.0])
                .unwrap()
                .with_occupancy(0.5)
                .unwrap()
                .with_disorder_group(7),
        )
        .build()
        .unwrap()
}

#[test]
fn groups_merge_occupants_at_one_position() {
    let groups = disorder_groups(&zinc_oxide_alloy()).unwrap();
    assert_eq!(groups.len(), 1);
    let group = groups.first().unwrap();
    assert_eq!(group.disorder_group, Some(7));
    assert_eq!(group.positions.len(), 1);
    let position = group.positions.first().unwrap();
    assert_eq!(position.occupants.len(), 2);
    assert_eq!(position.occupants[0].element, "O");
    assert_eq!(position.occupants[1].element, "Zn");
    assert!((position.total_occupancy() - 1.0).abs() < 1e-12);
}

#[test]
fn composition_formula_and_density_are_occupancy_weighted() {
    let structure = zinc_oxide_alloy();
    let composition = disorder_composition(&structure).unwrap();
    assert_eq!(composition.elements.len(), 2);
    assert!((composition.elements[0].occupancy - 0.5).abs() < 1e-12);
    assert_eq!(disorder_formula(&structure).unwrap(), "Zn0.500 O0.500");
    let expected = structure.density().unwrap();
    let density = disorder_density(&structure).unwrap();
    assert!((density - expected).abs() < 1e-12);
    assert!(composition.density(0.0).is_err());
}

#[test]
fn a_vacancy_is_a_disorder_group() {
    let structure = StructureBuilder::new("vacancy", Lattice::cubic(4.0).unwrap())
        .site(
            Site::new("Zn1", "Zn", [0.0, 0.0, 0.0])
                .unwrap()
                .with_occupancy(0.4)
                .unwrap(),
        )
        .build()
        .unwrap();
    let groups = disorder_groups(&structure).unwrap();
    assert_eq!(groups.len(), 1);
    let group = groups.first().unwrap();
    assert_eq!(group.positions.len(), 1);
    assert!((group.positions[0].total_occupancy() - 0.4).abs() < 1e-12);
    assert!(group.positions[0].is_disordered());
}

#[test]
fn explicit_groups_can_span_two_split_positions() {
    let structure = StructureBuilder::new("split", Lattice::cubic(4.0).unwrap())
        .site(
            Site::new("A1", "Na", [0.0, 0.0, 0.0])
                .unwrap()
                .with_occupancy(0.5)
                .unwrap()
                .with_disorder_group(3),
        )
        .site(
            Site::new("B1", "K", [0.0, 0.0, 0.0])
                .unwrap()
                .with_occupancy(0.5)
                .unwrap()
                .with_disorder_group(3),
        )
        .site(
            Site::new("A2", "Na", [0.5, 0.5, 0.5])
                .unwrap()
                .with_occupancy(0.5)
                .unwrap()
                .with_disorder_group(3),
        )
        .site(
            Site::new("B2", "K", [0.5, 0.5, 0.5])
                .unwrap()
                .with_occupancy(0.5)
                .unwrap()
                .with_disorder_group(3),
        )
        .build()
        .unwrap();
    let groups = disorder_groups(&structure).unwrap();
    assert_eq!(groups.len(), 1);
    let group = groups.first().unwrap();
    assert_eq!(group.positions.len(), 2);
    let occupants = group.aggregate_occupants();
    assert_eq!(occupants.len(), 2);
    assert!((occupants[0].occupancy - 1.0).abs() < 1e-12);
    assert!((occupants[1].occupancy - 1.0).abs() < 1e-12);
}

#[test]
fn xray_weights_are_the_occupancy_weighted_amplitudes() {
    let structure = zinc_oxide_alloy();
    let weights = disorder_scattering_weights(&structure, RadiationKind::XRay, 0.2).unwrap();
    assert_eq!(weights.scatterers.len(), 1);
    let scatterer = weights.scatterers.first().unwrap();
    let expected =
        0.5 * (xray_form_factor("Zn", 0.2).unwrap() + xray_form_factor("O", 0.2).unwrap());
    assert!((scatterer.coherent_amplitude - expected).abs() < 1e-12);
    assert!(scatterer.incoherent_intensity > 0.0);
    assert!((weights.total_occupancy - 1.0).abs() < 1e-12);
}

#[test]
fn neutron_weights_use_scattering_lengths() {
    let structure = zinc_oxide_alloy();
    let weights = disorder_scattering_weights(&structure, RadiationKind::Neutron, 0.0).unwrap();
    let scatterer = weights.scatterers.first().unwrap();
    let expected =
        0.5 * (neutron_scattering_length("Zn").unwrap() + neutron_scattering_length("O").unwrap());
    assert!((scatterer.coherent_amplitude - expected).abs() < 1e-12);
    assert!(disorder_scattering_weights(&structure, RadiationKind::XRay, -1.0).is_err());
}
