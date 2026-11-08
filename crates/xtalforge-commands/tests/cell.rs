//! Integration tests for cell commands (Task 3.10.3).

#![allow(
    clippy::unwrap_used,
    clippy::expect_used,
    clippy::panic,
    clippy::indexing_slicing
)]

use xtalforge_commands::cell::{PositionPolicy, SetLattice, ShiftOrigin, Supercell, WrapIntoCell};
use xtalforge_commands::command::{Command, CommandStack, Origin};
use xtalforge_core::{Lattice, Site, Structure, StructureBuilder, SymOp};

fn base() -> Structure {
    StructureBuilder::new("cell", Lattice::cubic(4.0).unwrap())
        .metadata("source", "test".into())
        .site(Site::new("C1", "C", [0.25, 0.5, 0.5]).unwrap())
        .site(Site::new("O1", "O", [0.5, 0.25, 0.5]).unwrap())
        .build()
        .unwrap()
}

fn assert_lengths(structure: &Structure, a: f64) {
    assert!((structure.lattice().lengths()[0] - a).abs() < 1e-9);
}

fn assert_frac(structure: &Structure, index: usize, expected: [f64; 3]) {
    let actual = structure.sites()[index].frac();
    for (actual, expected) in actual.iter().zip(expected.iter()) {
        assert!((actual - expected).abs() < 1e-12);
    }
}

#[test]
fn set_lattice_keeps_fractional_positions() {
    let mut structure = base();
    let mut stack = CommandStack::new();
    let new_lattice = Lattice::cubic(8.0).unwrap();
    let command = SetLattice::new(new_lattice, PositionPolicy::KeepFractional);
    assert_eq!(command.keep(), PositionPolicy::KeepFractional);
    stack
        .push(&mut structure, Box::new(command), Origin::User)
        .unwrap();
    assert_lengths(&structure, 8.0);
    assert_frac(&structure, 0, [0.25, 0.5, 0.5]);
    assert_eq!(
        structure.metadata().get("source").unwrap().as_str(),
        Some("test")
    );

    assert!(stack.undo(&mut structure));
    assert_lengths(&structure, 4.0);
    assert!(stack.redo(&mut structure));
    assert_lengths(&structure, 8.0);
}

#[test]
fn set_lattice_keeps_cartesian_positions() {
    let mut structure = base();
    let mut stack = CommandStack::new();
    let cart_before = structure.lattice().to_cart(structure.sites()[0].frac());
    let new_lattice = Lattice::cubic(8.0).unwrap();
    stack
        .push(
            &mut structure,
            Box::new(SetLattice::new(new_lattice, PositionPolicy::KeepCartesian)),
            Origin::Script,
        )
        .unwrap();
    assert_frac(&structure, 0, [0.125, 0.25, 0.25]);
    let cart_after = structure.lattice().to_cart(structure.sites()[0].frac());
    for (before, after) in cart_before.iter().zip(cart_after.iter()) {
        assert!((before - after).abs() < 1e-9);
    }

    assert!(stack.undo(&mut structure));
    assert_lengths(&structure, 4.0);
    assert_frac(&structure, 0, [0.25, 0.5, 0.5]);
}

#[test]
fn supercell_builds_p1_and_undoes() {
    let mut structure = base();
    let mut stack = CommandStack::new();
    let command = Supercell::new(2, 2, 2);
    assert_eq!(command.counts(), (2, 2, 2));
    stack
        .push(&mut structure, Box::new(command), Origin::User)
        .unwrap();
    assert!(structure.is_p1());
    assert_eq!(structure.sites().len(), 16);
    assert!((structure.lattice().volume() - 512.0).abs() < 1e-9);

    assert!(stack.undo(&mut structure));
    assert_eq!(structure.sites().len(), 2);
    assert!((structure.lattice().volume() - 64.0).abs() < 1e-9);
    assert!(stack.redo(&mut structure));
    assert_eq!(structure.sites().len(), 16);

    assert!(Supercell::new(0, 1, 1).apply(&structure).is_err());
}

#[test]
fn shift_origin_moves_and_wraps_sites() {
    let mut structure = base();
    let mut stack = CommandStack::new();
    let command = ShiftOrigin::new([0.5, 0.0, 0.0]);
    assert_eq!(command.delta(), [0.5, 0.0, 0.0]);
    stack
        .push(&mut structure, Box::new(command), Origin::User)
        .unwrap();
    assert_frac(&structure, 0, [0.75, 0.5, 0.5]);
    assert_frac(&structure, 1, [0.0, 0.25, 0.5]);
    assert_eq!(
        structure.metadata().get("source").unwrap().as_str(),
        Some("test")
    );

    assert!(stack.undo(&mut structure));
    assert_frac(&structure, 0, [0.25, 0.5, 0.5]);

    assert!(
        ShiftOrigin::new([f64::NAN, 0.0, 0.0])
            .apply(&structure)
            .is_err()
    );
}

#[test]
fn shift_origin_requires_p1() {
    let structure = StructureBuilder::new("inverted", Lattice::cubic(4.0).unwrap())
        .hall("P -1")
        .operation(SymOp::parse("-x,-y,-z").unwrap())
        .site(Site::new("C1", "C", [0.1, 0.1, 0.1]).unwrap())
        .build()
        .unwrap();
    assert!(ShiftOrigin::new([0.1, 0.0, 0.0]).apply(&structure).is_err());
}

#[test]
fn wrap_into_cell_folds_sites_and_undoes() {
    let mut structure = StructureBuilder::new("wrapped", Lattice::cubic(4.0).unwrap())
        .site(Site::new("C1", "C", [1.2, -0.2, 0.5]).unwrap())
        .build()
        .unwrap();
    let mut stack = CommandStack::new();
    stack
        .push(&mut structure, Box::new(WrapIntoCell), Origin::Script)
        .unwrap();
    assert_frac(&structure, 0, [0.2, 0.8, 0.5]);
    assert!(stack.undo(&mut structure));
    assert_frac(&structure, 0, [1.2, -0.2, 0.5]);
}
