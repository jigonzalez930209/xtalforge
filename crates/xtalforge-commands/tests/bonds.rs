//! Integration tests for bond commands (Task 3.10.2).

#![allow(
    clippy::unwrap_used,
    clippy::expect_used,
    clippy::panic,
    clippy::indexing_slicing
)]

use xtalforge_commands::bonds::{AddBond, RemoveBond, SetBondRules};
use xtalforge_commands::command::{Command, CommandStack, Origin};
use xtalforge_core::{ExplicitBond, Lattice, Site, Structure, StructureBuilder};

fn base() -> Structure {
    StructureBuilder::new("bonds", Lattice::cubic(4.0).unwrap())
        .site(Site::new("C1", "C", [0.0, 0.0, 0.0]).unwrap())
        .site(Site::new("O1", "O", [0.5, 0.5, 0.5]).unwrap())
        .build()
        .unwrap()
}

#[test]
fn add_bond_stores_an_explicit_bond_and_undoes() {
    let mut structure = base();
    let mut stack = CommandStack::new();
    let command = AddBond::between(0, 1);
    assert_eq!(command.bond().i, 0);
    stack
        .push(&mut structure, Box::new(command), Origin::User)
        .unwrap();
    assert_eq!(structure.explicit_bonds().len(), 1);
    assert_eq!(structure.explicit_bonds()[0].order, 1.0);
    assert!(!structure.explicit_bonds()[0].remove);

    assert!(stack.undo(&mut structure));
    assert!(structure.explicit_bonds().is_empty());
    assert!(stack.redo(&mut structure));
    assert_eq!(structure.explicit_bonds().len(), 1);
}

#[test]
fn add_bond_rejects_invalid_order_and_indexes() {
    let structure = base();
    let bad_order = AddBond::new(ExplicitBond {
        i: 0,
        j: 1,
        image: [0, 0, 0],
        order: 0.0,
        remove: false,
    });
    assert!(bad_order.apply(&structure).is_err());

    let out_of_range = AddBond::between(0, 8);
    assert!(out_of_range.apply(&structure).is_err());
}

#[test]
fn remove_bond_stores_a_suppression_and_add_replaces_it() {
    let mut structure = base();
    let mut stack = CommandStack::new();
    stack
        .push(
            &mut structure,
            Box::new(RemoveBond::new(0, 1, [0, 0, 0])),
            Origin::User,
        )
        .unwrap();
    assert_eq!(structure.explicit_bonds().len(), 1);
    assert!(structure.explicit_bonds()[0].remove);

    assert!(stack.undo(&mut structure));
    assert!(structure.explicit_bonds().is_empty());

    stack
        .push(
            &mut structure,
            Box::new(RemoveBond::new(0, 1, [0, 0, 0])),
            Origin::User,
        )
        .unwrap();
    stack
        .push(
            &mut structure,
            Box::new(AddBond::between(0, 1)),
            Origin::Agent,
        )
        .unwrap();
    assert_eq!(structure.explicit_bonds().len(), 1);
    assert!(!structure.explicit_bonds()[0].remove);

    assert!(stack.undo(&mut structure));
    assert!(structure.explicit_bonds()[0].remove);
}

#[test]
fn set_bond_rules_writes_metadata_and_undoes() {
    let mut structure = base();
    let mut stack = CommandStack::new();
    stack
        .push(
            &mut structure,
            Box::new(SetBondRules::new(1.35, 0.5)),
            Origin::Script,
        )
        .unwrap();
    let rules = structure.metadata().get("bond_rules").unwrap();
    assert_eq!(rules["scale"], 1.35);
    assert_eq!(rules["minimum_distance"], 0.5);

    assert!(stack.undo(&mut structure));
    assert!(!structure.metadata().contains_key("bond_rules"));
    assert!(stack.redo(&mut structure));
    assert!(structure.metadata().contains_key("bond_rules"));

    assert!(SetBondRules::new(0.0, 0.4).apply(&structure).is_err());
    assert!(SetBondRules::new(1.2, -1.0).apply(&structure).is_err());
}

#[test]
fn remove_bond_range_is_validated() {
    let structure = base();
    let command = RemoveBond::new(0, 9, [0, 0, 0]);
    assert_eq!(command.i(), 0);
    assert_eq!(command.j(), 9);
    let error = command.apply(&structure).unwrap_err();
    assert_eq!(error.code(), "invalid_input");
}
