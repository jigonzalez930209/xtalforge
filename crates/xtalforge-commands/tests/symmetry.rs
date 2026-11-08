//! Integration tests for symmetry commands (Task 3.10.4).

#![allow(
    clippy::unwrap_used,
    clippy::expect_used,
    clippy::panic,
    clippy::indexing_slicing
)]

use xtalforge_commands::command::{Command, CommandStack, Origin};
use xtalforge_commands::symmetry::{FindSymmetry, ReduceToP1, SetSpaceGroup, Standardize};
use xtalforge_core::{Lattice, Site, Structure, StructureBuilder, SymOp};
use xtalforge_symmetry::SetMode;

fn nacl() -> Structure {
    let mut builder = StructureBuilder::new("NaCl", Lattice::cubic(5.64).unwrap());
    for position in [
        [0.0, 0.0, 0.0],
        [0.0, 0.5, 0.5],
        [0.5, 0.0, 0.5],
        [0.5, 0.5, 0.0],
    ] {
        builder = builder.site(Site::new("Na1", "Na", position).unwrap());
    }
    for position in [
        [0.5, 0.5, 0.5],
        [0.5, 0.0, 0.0],
        [0.0, 0.5, 0.0],
        [0.0, 0.0, 0.5],
    ] {
        builder = builder.site(Site::new("Cl1", "Cl", position).unwrap());
    }
    builder.build().unwrap()
}

fn inverted() -> Structure {
    StructureBuilder::new("inverted", Lattice::cubic(4.0).unwrap())
        .hall("P -1")
        .operation(SymOp::identity())
        .operation(SymOp::parse("-x,-y,-z").unwrap())
        .site(Site::new("C1", "C", [0.1, 0.2, 0.3]).unwrap())
        .build()
        .unwrap()
}

#[test]
fn reduce_to_p1_expands_orbits_and_undoes() {
    let mut structure = inverted();
    let mut stack = CommandStack::new();
    stack
        .push(&mut structure, Box::new(ReduceToP1), Origin::User)
        .unwrap();
    assert!(structure.is_p1());
    assert_eq!(structure.sites().len(), 2);
    assert!(structure.operations().is_empty());
    assert_eq!(ReduceToP1.description(), "Reduce to P1");

    assert!(stack.undo(&mut structure));
    assert_eq!(structure.sites().len(), 1);
    assert_eq!(structure.hall(), "P -1");
    assert!(stack.redo(&mut structure));
    assert_eq!(structure.sites().len(), 2);
}

#[test]
fn find_symmetry_reduces_the_nacl_cell() {
    let mut structure = nacl();
    let mut stack = CommandStack::new();
    stack
        .push(
            &mut structure,
            Box::new(FindSymmetry::new()),
            Origin::Script,
        )
        .unwrap();
    assert_eq!(structure.sites().len(), 2);
    assert_ne!(structure.hall(), "P 1");
    assert!(!structure.operations().is_empty());

    assert!(stack.undo(&mut structure));
    assert_eq!(structure.sites().len(), 8);
    assert_eq!(structure.hall(), "P 1");
}

#[test]
fn find_symmetry_honours_a_tolerance() {
    let command = FindSymmetry::with_tolerance(0.005);
    assert!((command.options().symprec - 0.005).abs() < 1e-12);
    assert_eq!(command.description(), "Find symmetry");
    let error = FindSymmetry::with_tolerance(0.0)
        .apply(&nacl())
        .unwrap_err();
    assert_eq!(error.code(), "symmetry");
}

#[test]
fn set_space_group_reinterprets_the_sites() {
    assert!(SetSpaceGroup::from_number(0, SetMode::Reinterpret).is_err());
    let structure = StructureBuilder::new("rock salt", Lattice::cubic(5.64).unwrap())
        .site(Site::new("Na1", "Na", [0.0, 0.0, 0.0]).unwrap())
        .site(Site::new("Cl1", "Cl", [0.5, 0.5, 0.5]).unwrap())
        .build()
        .unwrap();
    let mut structure = structure;
    let mut stack = CommandStack::new();
    let command = SetSpaceGroup::from_number(225, SetMode::Reinterpret).unwrap();
    assert_eq!(command.mode(), SetMode::Reinterpret);
    stack
        .push(&mut structure, Box::new(command), Origin::User)
        .unwrap();
    assert_ne!(structure.hall(), "P 1");
    assert!(!structure.operations().is_empty());
    assert_eq!(structure.sites().len(), 2);

    assert!(stack.undo(&mut structure));
    assert_eq!(structure.hall(), "P 1");
    assert!(stack.redo(&mut structure));
    assert_ne!(structure.hall(), "P 1");
}

#[test]
fn standardize_rebuilds_the_conventional_cell() {
    let mut structure = nacl();
    let mut stack = CommandStack::new();
    let command = Standardize::new();
    assert!(!command.primitive());
    stack
        .push(&mut structure, Box::new(command), Origin::Agent)
        .unwrap();
    assert_eq!(structure.sites().len(), 2);
    assert_ne!(structure.hall(), "P 1");
    assert!((structure.lattice().volume() - 5.64_f64.powi(3)).abs() < 1e-3);

    assert!(stack.undo(&mut structure));
    assert_eq!(structure.sites().len(), 8);
    assert_eq!(structure.hall(), "P 1");
}
