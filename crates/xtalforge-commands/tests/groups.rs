//! Integration tests for atom-group commands (Task 3.10.5).

#![allow(
    clippy::unwrap_used,
    clippy::expect_used,
    clippy::panic,
    clippy::indexing_slicing
)]

use xtalforge_commands::atoms::MoveSites;
use xtalforge_commands::command::Origin;
use xtalforge_commands::document::Document;
use xtalforge_commands::groups::{
    AddAtomGroup, AtomGroup, GroupCommand, GroupStore, RemoveAtomGroup, SetGroupColor,
    SetGroupHidden,
};
use xtalforge_core::{Lattice, Site, Structure, StructureBuilder};

fn base() -> Structure {
    StructureBuilder::new("groups", Lattice::cubic(4.0).unwrap())
        .site(Site::new("C1", "C", [0.0, 0.0, 0.0]).unwrap())
        .site(Site::new("O1", "O", [0.5, 0.5, 0.5]).unwrap())
        .build()
        .unwrap()
}

fn group() -> AtomGroup {
    AtomGroup::new("ligand", vec![0, 1]).unwrap()
}

#[test]
fn add_and_remove_groups_round_trip() {
    let mut document = Document::new(base());
    document
        .apply_group(Box::new(AddAtomGroup::new(group())), Origin::User)
        .unwrap();
    assert_eq!(document.groups().len(), 1);
    assert_eq!(document.groups()[0].name, "ligand");
    assert!(document.dirty());

    assert!(document.undo());
    assert!(document.groups().is_empty());
    assert!(!document.dirty());
    assert!(document.redo());
    assert_eq!(document.groups().len(), 1);

    document
        .apply_group(Box::new(RemoveAtomGroup::new(0)), Origin::Script)
        .unwrap();
    assert!(document.groups().is_empty());
    assert!(document.undo());
    assert_eq!(document.groups().len(), 1);
}

#[test]
fn group_commands_validate_names_atoms_and_rows() {
    let mut document = Document::new(base());
    let blank = AtomGroup {
        name: "  ".to_owned(),
        atoms: vec![0],
        color: None,
        hidden: false,
    };
    assert!(
        document
            .apply_group(Box::new(AddAtomGroup::new(blank)), Origin::User)
            .is_err()
    );

    let out_of_range = AtomGroup::new("bad", vec![7]).unwrap();
    assert!(
        document
            .apply_group(Box::new(AddAtomGroup::new(out_of_range)), Origin::User)
            .is_err()
    );

    let duplicate = AtomGroup::new("twice", vec![0, 0]).unwrap();
    assert!(
        document
            .apply_group(Box::new(AddAtomGroup::new(duplicate)), Origin::User)
            .is_err()
    );

    assert!(
        document
            .apply_group(Box::new(RemoveAtomGroup::new(0)), Origin::User)
            .is_err()
    );
    assert!(
        document
            .apply_group(
                Box::new(SetGroupColor::new(0, Some([1, 2, 3]))),
                Origin::User
            )
            .is_err()
    );
    assert!(!document.dirty());
}

#[test]
fn colour_and_visibility_round_trip() {
    let mut document = Document::new(base());
    document
        .apply_group(Box::new(AddAtomGroup::new(group())), Origin::User)
        .unwrap();
    document.mark_clean();
    assert!(!document.dirty());

    document
        .apply_group(
            Box::new(SetGroupColor::new(0, Some([255, 0, 0]))),
            Origin::User,
        )
        .unwrap();
    assert_eq!(document.groups()[0].color, Some([255, 0, 0]));
    assert!(document.dirty());
    assert!(document.undo());
    assert_eq!(document.groups()[0].color, None);
    assert!(!document.dirty());
    assert!(document.redo());
    assert_eq!(document.groups()[0].color, Some([255, 0, 0]));

    document
        .apply_group(Box::new(SetGroupHidden::new(0, true)), Origin::User)
        .unwrap();
    assert!(document.groups()[0].hidden);
    assert!(document.undo());
    assert!(!document.groups()[0].hidden);
    assert!(document.redo());
    assert!(document.groups()[0].hidden);
}

#[test]
fn group_and_structure_edits_share_one_timeline() {
    let mut document = Document::new(base());
    document
        .apply(
            Box::new(MoveSites::new(vec![(0, [0.1, 0.1, 0.1])])),
            Origin::User,
        )
        .unwrap();
    document
        .apply_group(Box::new(AddAtomGroup::new(group())), Origin::Agent)
        .unwrap();
    assert_eq!(document.groups().len(), 1);
    assert_eq!(document.structure().sites()[0].frac(), [0.1, 0.1, 0.1]);

    assert!(document.undo());
    assert!(document.groups().is_empty());
    assert_eq!(document.structure().sites()[0].frac(), [0.1, 0.1, 0.1]);
    assert!(document.undo());
    assert_eq!(document.structure().sites()[0].frac(), [0.0, 0.0, 0.0]);
    assert!(!document.undo());

    assert!(document.redo());
    assert_eq!(document.structure().sites()[0].frac(), [0.1, 0.1, 0.1]);
    assert!(document.redo());
    assert_eq!(document.groups().len(), 1);
    assert!(!document.redo());
}

#[test]
fn group_edit_after_undo_does_not_inherit_the_clean_mark() {
    let mut document = Document::new(base());
    document
        .apply_group(Box::new(AddAtomGroup::new(group())), Origin::User)
        .unwrap();
    document.mark_clean();
    assert!(document.undo());
    let other = AtomGroup::new("other", vec![1]).unwrap();
    document
        .apply_group(Box::new(AddAtomGroup::new(other)), Origin::User)
        .unwrap();
    assert!(document.dirty());
}

#[test]
fn group_store_reports_range_errors_directly() {
    let mut store = GroupStore::new();
    assert!(store.is_empty());
    AddAtomGroup::new(group()).apply(&mut store).unwrap();
    assert_eq!(store.len(), 1);
    assert!(RemoveAtomGroup::new(5).apply(&mut store).is_err());
    assert!(
        SetGroupColor::new(5, Some([0, 0, 0]))
            .apply(&mut store)
            .is_err()
    );
    assert!(SetGroupHidden::new(5, true).apply(&mut store).is_err());

    let description = AddAtomGroup::new(group()).description();
    assert_eq!(description, "Add atom group ligand");
    assert_eq!(
        SetGroupHidden::new(0, true).description(),
        "Hide atom group"
    );
}
