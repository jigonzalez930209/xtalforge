//! Integration tests for force-field command wrappers (Task 14.3.3).

#![allow(
    clippy::unwrap_used,
    clippy::expect_used,
    clippy::panic,
    clippy::indexing_slicing
)]

use xtalforge_commands::command::{Command, CommandStack, Origin};
use xtalforge_commands::ff::{
    ATOM_TYPE_KEY, AddHydrogens, ApplyOptimizedGeometry, Prepare, PrepareStep, SetAtomTypes,
};
use xtalforge_core::{Lattice, Site, Structure, StructureBuilder, SymOp};

fn base() -> Structure {
    StructureBuilder::new("ff", Lattice::cubic(10.0).unwrap())
        .site(Site::new("C1", "C", [0.5, 0.5, 0.5]).unwrap())
        .site(Site::new("O1", "O", [0.6, 0.5, 0.5]).unwrap())
        .build()
        .unwrap()
}

#[test]
fn optimized_geometry_moves_sites_and_carries_the_cell() {
    let mut structure = base();
    let mut stack = CommandStack::new();
    let relaxed = Lattice::cubic(10.5).unwrap();
    let command = ApplyOptimizedGeometry::new(vec![[0.5, 0.5, 0.5], [0.62, 0.5, 0.5]])
        .with_lattice(relaxed)
        .with_description("Optimise with UFF");
    assert_eq!(command.fracs().len(), 2);
    assert!(command.lattice().is_some());
    assert_eq!(command.description(), "Optimise with UFF");
    assert!(command.max_displacement(&structure).unwrap() > 0.19);
    stack
        .push(&mut structure, Box::new(command), Origin::Agent)
        .unwrap();
    assert!((structure.sites()[1].frac()[0] - 0.62).abs() < 1e-12);
    assert!((structure.lattice().lengths()[0] - 10.5).abs() < 1e-12);

    assert!(stack.undo(&mut structure));
    assert!((structure.sites()[1].frac()[0] - 0.6).abs() < 1e-12);
    assert!((structure.lattice().lengths()[0] - 10.0).abs() < 1e-12);
    assert!(stack.redo(&mut structure));
    assert!((structure.lattice().lengths()[0] - 10.5).abs() < 1e-12);

    let wrong = ApplyOptimizedGeometry::new(vec![[0.5, 0.5, 0.5]]);
    assert!(wrong.apply(&structure).is_err());
    assert!(wrong.max_displacement(&structure).is_err());
}

#[test]
fn add_hydrogens_appends_sites_and_bonds_and_undoes() {
    let mut structure = base();
    let mut stack = CommandStack::new();
    let hydrogen = Site::new("H1", "H", [0.55, 0.6, 0.5]).unwrap();
    let command = AddHydrogens::new(vec![hydrogen]).with_bond(0, 0);
    assert_eq!(command.sites().len(), 1);
    assert_eq!(command.bonds(), &[(0, 0)]);
    stack
        .push(&mut structure, Box::new(command), Origin::Script)
        .unwrap();
    assert_eq!(structure.sites().len(), 3);
    assert_eq!(structure.sites()[2].element(), "H");
    assert_eq!(structure.explicit_bonds().len(), 1);
    assert_eq!(
        (
            structure.explicit_bonds()[0].i,
            structure.explicit_bonds()[0].j
        ),
        (0, 2)
    );

    assert!(stack.undo(&mut structure));
    assert_eq!(structure.sites().len(), 2);
    assert!(structure.explicit_bonds().is_empty());

    let oxygen = Site::new("O2", "O", [0.0, 0.0, 0.0]).unwrap();
    assert!(AddHydrogens::new(vec![oxygen]).apply(&structure).is_err());
    let hydrogen = Site::new("H1", "H", [0.0, 0.0, 0.0]).unwrap();
    assert!(AddHydrogens::new(Vec::new()).apply(&structure).is_err());
    assert!(
        AddHydrogens::new(vec![hydrogen.clone()])
            .with_bond(9, 0)
            .apply(&structure)
            .is_err()
    );
    assert!(
        AddHydrogens::new(vec![hydrogen])
            .with_bond(0, 4)
            .apply(&structure)
            .is_err()
    );
}

#[test]
fn add_hydrogens_requires_p1_for_bonds() {
    let structure = StructureBuilder::new("inverted", Lattice::cubic(10.0).unwrap())
        .hall("P -1")
        .operation(SymOp::identity())
        .operation(SymOp::parse("-x,-y,-z").unwrap())
        .site(Site::new("C1", "C", [0.1, 0.1, 0.1]).unwrap())
        .build()
        .unwrap();
    let hydrogen = Site::new("H1", "H", [0.2, 0.1, 0.1]).unwrap();
    assert!(
        AddHydrogens::new(vec![hydrogen.clone()])
            .with_bond(0, 0)
            .apply(&structure)
            .is_err()
    );
    assert!(AddHydrogens::new(vec![hydrogen]).apply(&structure).is_ok());
}

#[test]
fn prepare_installs_the_pipeline_result_as_one_step() {
    let original = base();
    let mut prepared = original.clone();
    prepared = prepared.with_site(Site::new("H2", "H", [0.7, 0.5, 0.5]).unwrap());
    let steps = vec![
        PrepareStep::new("minimize", true, "converged in 12 steps"),
        PrepareStep::new("charges", true, "QEq charges stored"),
        PrepareStep::new("ewald", false, "skipped: disabled"),
    ];
    let command = Prepare::new(prepared.clone(), steps.clone()).unwrap();
    assert_eq!(command.steps(), steps.as_slice());
    assert_eq!(command.prepared().sites().len(), 3);
    assert_eq!(command.description(), "Prepare for simulation");

    let mut structure = original.clone();
    let mut stack = CommandStack::new();
    stack
        .push(&mut structure, Box::new(command), Origin::Script)
        .unwrap();
    assert_eq!(structure.sites().len(), 3);
    assert!(stack.undo(&mut structure));
    assert_eq!(structure, original);
    assert!(stack.redo(&mut structure));
    assert_eq!(structure.sites().len(), 3);

    assert!(Prepare::new(prepared.clone(), Vec::new()).is_err());
    assert!(
        Prepare::new(
            prepared,
            vec![PrepareStep::new("minimize", false, "skipped")]
        )
        .is_err()
    );
}

#[test]
fn atom_types_are_stored_cleared_and_undone() {
    let mut structure = base();
    let mut stack = CommandStack::new();
    let command = SetAtomTypes::new(vec![1], Some("O_R"));
    assert_eq!(command.indices(), &[1]);
    assert_eq!(command.type_name(), Some("O_R"));
    assert_eq!(command.description(), "Set atom type to O_R");
    stack
        .push(&mut structure, Box::new(command), Origin::User)
        .unwrap();
    assert_eq!(
        structure.sites()[1]
            .props()
            .get(ATOM_TYPE_KEY)
            .unwrap()
            .as_str(),
        Some("O_R")
    );

    let clear = SetAtomTypes::new(vec![1], None);
    assert_eq!(clear.description(), "Clear atom type override");
    stack
        .push(&mut structure, Box::new(clear), Origin::User)
        .unwrap();
    assert!(!structure.sites()[1].props().contains_key(ATOM_TYPE_KEY));

    assert!(stack.undo(&mut structure));
    assert!(structure.sites()[1].props().contains_key(ATOM_TYPE_KEY));
    assert!(
        SetAtomTypes::new(vec![7], Some("C3"))
            .apply(&structure)
            .is_err()
    );
    assert_eq!(ATOM_TYPE_KEY, "uff_type");
}
