//! Integration tests for atom commands (Task 3.10.1).

#![allow(
    clippy::unwrap_used,
    clippy::expect_used,
    clippy::panic,
    clippy::indexing_slicing
)]

use xtalforge_commands::atoms::{AddSites, DeleteSites, MoveSites, SetElement, SetOccupancy};
use xtalforge_commands::command::{Command, CommandStack, Origin};
use xtalforge_core::{ExplicitBond, Lattice, Site, Structure, StructureBuilder};

fn base() -> Structure {
    let oxygen = Site::new("O1", "O", [0.5, 0.5, 0.5])
        .unwrap()
        .with_occupancy(0.5)
        .unwrap()
        .with_charge(-1.0)
        .unwrap();
    StructureBuilder::new("atoms", Lattice::cubic(4.0).unwrap())
        .metadata("source", "test".into())
        .site(Site::new("C1", "C", [0.0, 0.0, 0.0]).unwrap())
        .site(oxygen)
        .explicit_bond(ExplicitBond {
            i: 0,
            j: 1,
            image: [0, 0, 0],
            order: 1.0,
            remove: false,
        })
        .build()
        .unwrap()
}

#[test]
fn add_sites_appends_and_undoes() {
    let mut structure = base();
    let mut stack = CommandStack::new();
    let command = AddSites::new(vec![Site::new("H1", "H", [0.25, 0.25, 0.25]).unwrap()]);
    assert_eq!(command.description(), "Add H");
    stack
        .push(&mut structure, Box::new(command), Origin::User)
        .unwrap();
    assert_eq!(structure.sites().len(), 3);
    assert_eq!(structure.sites()[2].label(), "H1");
    assert_eq!(
        structure.metadata().get("source").unwrap().as_str(),
        Some("test")
    );
    assert_eq!(structure.explicit_bonds().len(), 1);

    assert!(stack.undo(&mut structure));
    assert_eq!(structure.sites().len(), 2);
    assert!(stack.redo(&mut structure));
    assert_eq!(structure.sites().len(), 3);
}

#[test]
fn delete_sites_validates_ranges_and_remaps_bonds() {
    let structure = base();
    let error = DeleteSites::new(vec![7]).apply(&structure).unwrap_err();
    assert_eq!(error.code(), "invalid_input");

    let mut structure = base();
    let mut stack = CommandStack::new();
    let command = DeleteSites::new(vec![1, 1]);
    assert_eq!(command.indices(), &[1, 1]);
    stack
        .push(&mut structure, Box::new(command), Origin::Script)
        .unwrap();
    assert_eq!(structure.sites().len(), 1);
    assert!(structure.explicit_bonds().is_empty());
    assert_eq!(structure.sites()[0].label(), "C1");

    assert!(stack.undo(&mut structure));
    assert_eq!(structure.sites().len(), 2);
    assert_eq!(structure.explicit_bonds().len(), 1);
}

#[test]
fn set_element_keeps_site_properties_and_round_trips() {
    assert!(SetElement::new(vec![0], "Qq").is_err());
    let mut structure = base();
    let mut stack = CommandStack::new();
    let command = SetElement::new(vec![1], "S").unwrap();
    assert_eq!(command.description(), "Change element to S");
    stack
        .push(&mut structure, Box::new(command), Origin::User)
        .unwrap();
    assert_eq!(structure.sites()[1].element(), "S");
    assert_eq!(structure.sites()[1].occupancy(), 0.5);
    assert_eq!(structure.sites()[1].charge(), Some(-1.0));

    assert!(stack.undo(&mut structure));
    assert_eq!(structure.sites()[1].element(), "O");
    assert!(stack.redo(&mut structure));
    assert_eq!(structure.sites()[1].element(), "S");

    let range = SetElement::new(vec![3], "N").unwrap();
    assert!(range.apply(&structure).is_err());
}

#[test]
fn move_sites_round_trips_and_rejects_invalid_input() {
    let mut structure = base();
    let mut stack = CommandStack::new();
    let command = MoveSites::new(vec![(0, [0.25, 0.25, 0.25])]);
    assert_eq!(command.description(), "Move site");
    stack
        .push(&mut structure, Box::new(command), Origin::User)
        .unwrap();
    assert_eq!(structure.sites()[0].frac(), [0.25, 0.25, 0.25]);
    assert!(stack.undo(&mut structure));
    assert_eq!(structure.sites()[0].frac(), [0.0, 0.0, 0.0]);
    assert!(stack.redo(&mut structure));
    assert_eq!(structure.sites()[0].frac(), [0.25, 0.25, 0.25]);

    let out_of_range = MoveSites::new(vec![(5, [0.0, 0.0, 0.0])]);
    assert!(out_of_range.apply(&structure).is_err());
    let not_finite = MoveSites::new(vec![(0, [f64::NAN, 0.0, 0.0])]);
    assert!(not_finite.apply(&structure).is_err());
}

#[test]
fn set_occupancy_round_trips_and_rejects_invalid_values() {
    let mut structure = base();
    let mut stack = CommandStack::new();
    let command = SetOccupancy::new(0, 0.25);
    assert_eq!(command.index(), 0);
    assert_eq!(command.description(), "Set occupancy 0.25");
    stack
        .push(&mut structure, Box::new(command), Origin::Agent)
        .unwrap();
    assert_eq!(structure.sites()[0].occupancy(), 0.25);
    assert!(stack.undo(&mut structure));
    assert_eq!(structure.sites()[0].occupancy(), 1.0);
    assert!(stack.redo(&mut structure));
    assert_eq!(structure.sites()[0].occupancy(), 0.25);

    assert!(SetOccupancy::new(0, 1.5).apply(&structure).is_err());
    assert!(SetOccupancy::new(9, 0.5).apply(&structure).is_err());
}
