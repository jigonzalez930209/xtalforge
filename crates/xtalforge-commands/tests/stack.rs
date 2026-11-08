//! Integration tests for the command stack, helpers and document (Task 3.9).

#![allow(
    clippy::unwrap_used,
    clippy::expect_used,
    clippy::panic,
    clippy::indexing_slicing
)]

use xtalforge_commands::atoms::{AddSites, DeleteSites, MoveSites};
use xtalforge_commands::command::{
    Command, CommandRecord, CommandStack, MacroCommand, Origin, ReplaceStructure, SnapshotEdit,
    StructureOperation,
};
use xtalforge_commands::document::Document;
use xtalforge_core::{Error, Lattice, Site, Structure, StructureBuilder};

fn base() -> Structure {
    StructureBuilder::new("stack", Lattice::cubic(4.0).unwrap())
        .site(Site::new("C1", "C", [0.0, 0.0, 0.0]).unwrap())
        .site(Site::new("O1", "O", [0.5, 0.5, 0.5]).unwrap())
        .build()
        .unwrap()
}

fn move_first(target: [f64; 3]) -> MoveSites {
    MoveSites::new(vec![(0, target)])
}

#[test]
fn push_undo_redo_tracks_dirty_state_and_log() {
    let mut structure = base();
    let mut stack = CommandStack::new();
    assert!(!stack.can_undo());
    assert!(!stack.can_redo());
    assert!(!stack.is_dirty());

    stack
        .push(
            &mut structure,
            Box::new(move_first([0.25, 0.25, 0.25])),
            Origin::Script,
        )
        .unwrap();
    assert_eq!(structure.sites()[0].frac(), [0.25, 0.25, 0.25]);
    assert!(stack.can_undo());
    assert!(stack.is_dirty());
    assert_eq!(stack.undo_depth(), 1);
    assert_eq!(stack.log().len(), 1);
    assert_eq!(stack.log()[0].origin, Origin::Script);

    stack.mark_clean();
    assert!(!stack.is_dirty());

    assert!(stack.undo(&mut structure));
    assert_eq!(structure.sites()[0].frac(), [0.0, 0.0, 0.0]);
    assert!(stack.can_redo());
    assert!(stack.is_dirty());
    assert!(stack.log().is_empty());

    assert!(stack.redo(&mut structure));
    assert_eq!(structure.sites()[0].frac(), [0.25, 0.25, 0.25]);
    assert!(!stack.is_dirty());
}

#[test]
fn new_push_discards_the_redo_branch() {
    let mut structure = base();
    let mut stack = CommandStack::new();
    stack
        .push(
            &mut structure,
            Box::new(move_first([0.1, 0.0, 0.0])),
            Origin::User,
        )
        .unwrap();
    assert!(stack.undo(&mut structure));
    assert!(stack.can_redo());

    stack
        .push(
            &mut structure,
            Box::new(move_first([0.2, 0.0, 0.0])),
            Origin::User,
        )
        .unwrap();
    assert!(!stack.can_redo());
    assert!(!stack.redo(&mut structure));
    assert_eq!(structure.sites()[0].frac(), [0.2, 0.0, 0.0]);
}

#[test]
fn push_after_undo_does_not_inherit_the_clean_mark() {
    let mut structure = base();
    let mut stack = CommandStack::new();
    stack
        .push(
            &mut structure,
            Box::new(move_first([0.1, 0.0, 0.0])),
            Origin::User,
        )
        .unwrap();
    stack.mark_clean();
    assert!(stack.undo(&mut structure));

    stack
        .push(
            &mut structure,
            Box::new(move_first([0.2, 0.0, 0.0])),
            Origin::User,
        )
        .unwrap();
    assert!(stack.is_dirty());
    assert!(stack.undo(&mut structure));
    assert!(stack.is_dirty());
}

#[test]
fn max_depth_drops_the_oldest_step() {
    let mut structure = base();
    let mut stack = CommandStack::with_max_depth(2);
    assert_eq!(stack.max_depth(), 2);
    for x in [0.1, 0.2, 0.3] {
        stack
            .push(
                &mut structure,
                Box::new(move_first([x, 0.0, 0.0])),
                Origin::User,
            )
            .unwrap();
    }
    assert_eq!(stack.trimmed(), 1);
    assert_eq!(stack.undo_depth(), 2);
    assert_eq!(stack.log().len(), 2);
    assert!(stack.undo(&mut structure));
    assert!(stack.undo(&mut structure));
    assert!(!stack.undo(&mut structure));
    assert_eq!(structure.sites()[0].frac(), [0.1, 0.0, 0.0]);
}

#[test]
fn macro_applies_in_order_and_undoes_as_one_step() {
    let mut structure = base();
    let mut macro_command = MacroCommand::new("Move both sites");
    macro_command.push(Box::new(move_first([0.1, 0.0, 0.0])));
    macro_command.push(Box::new(MoveSites::new(vec![(1, [0.6, 0.5, 0.5])])));
    assert_eq!(macro_command.len(), 2);
    assert!(!macro_command.is_empty());

    let mut stack = CommandStack::new();
    stack
        .push(&mut structure, Box::new(macro_command), Origin::User)
        .unwrap();
    assert_eq!(structure.sites()[0].frac(), [0.1, 0.0, 0.0]);
    assert_eq!(structure.sites()[1].frac(), [0.6, 0.5, 0.5]);
    assert_eq!(stack.undo_depth(), 1);

    assert!(stack.undo(&mut structure));
    assert_eq!(structure.sites()[0].frac(), [0.0, 0.0, 0.0]);
    assert_eq!(structure.sites()[1].frac(), [0.5, 0.5, 0.5]);
}

#[test]
fn snapshot_operation_and_replace_round_trip() {
    let mut structure = base();
    let mut stack = CommandStack::new();

    let snapshot = SnapshotEdit::new("Append nitrogen", |current| {
        let site = Site::new("N1", "N", [0.25, 0.25, 0.25])?;
        Ok(current.clone().with_site(site))
    });
    stack
        .push(&mut structure, Box::new(snapshot), Origin::User)
        .unwrap();
    assert_eq!(structure.sites().len(), 3);
    assert!(stack.undo(&mut structure));
    assert_eq!(structure.sites().len(), 2);

    let operation = StructureOperation::new("Double a", |current| {
        let [a, b, c, alpha, beta, gamma] = current.lattice().parameters();
        let lattice = Lattice::from_parameters(2.0 * a, b, c, alpha, beta, gamma)?;
        Ok(current.clone().with_lattice(lattice))
    });
    let preview = operation.preview(&structure).unwrap();
    assert!((preview.lattice().lengths()[0] - 8.0).abs() < 1e-9);
    stack
        .push(&mut structure, Box::new(operation), Origin::User)
        .unwrap();
    assert!((structure.lattice().lengths()[0] - 8.0).abs() < 1e-9);
    assert!(stack.undo(&mut structure));
    assert!((structure.lattice().lengths()[0] - 4.0).abs() < 1e-9);

    let replacement = base().with_name("replacement");
    stack
        .push(
            &mut structure,
            Box::new(ReplaceStructure::new(replacement)),
            Origin::Agent,
        )
        .unwrap();
    assert_eq!(structure.name(), "replacement");
    assert!(stack.undo(&mut structure));
    assert_eq!(structure.name(), "stack");
}

#[test]
fn failed_command_leaves_the_stack_untouched() {
    let mut structure = base();
    let mut stack = CommandStack::new();
    let failing = SnapshotEdit::new("Fail", |_| {
        Err(Error::InvalidInput("deliberate failure".to_owned()))
    });
    let error = stack
        .push(&mut structure, Box::new(failing), Origin::User)
        .unwrap_err();
    assert_eq!(error.code(), "invalid_input");
    assert!(!stack.is_dirty());
    assert!(!stack.can_undo());
    assert!(stack.log().is_empty());
    assert_eq!(structure.sites().len(), 2);
}

#[test]
fn command_log_is_serializable() {
    let mut structure = base();
    let mut stack = CommandStack::new();
    stack
        .push(
            &mut structure,
            Box::new(AddSites::new(vec![
                Site::new("H1", "H", [0.1, 0.1, 0.1]).unwrap(),
            ])),
            Origin::Agent,
        )
        .unwrap();
    let json = serde_json::to_string(stack.log()).unwrap();
    let decoded: Vec<CommandRecord> = serde_json::from_str(&json).unwrap();
    assert_eq!(decoded, stack.log());
    assert!(json.contains("agent"));
}

#[test]
fn document_applies_undoes_and_tracks_selection() {
    let mut document = Document::new(base());
    assert_eq!(document.name(), "stack");
    document.set_name("edited");
    assert_eq!(document.name(), "edited");
    assert!(!document.dirty());

    document.set_selection(vec![1, 1, 0]).unwrap();
    assert_eq!(document.selection(), &[1, 0]);
    assert!(document.set_selection(vec![9]).is_err());
    assert_eq!(document.selection(), &[1, 0]);

    document
        .apply(Box::new(move_first([0.3, 0.0, 0.0])), Origin::Agent)
        .unwrap();
    assert!(document.dirty());
    assert_eq!(document.origin(), Origin::Agent);
    assert!(document.can_undo());
    assert!(document.undo());
    assert!(!document.dirty());
    assert!(document.can_redo());
    assert!(document.redo());
    assert!(document.dirty());

    document.mark_clean();
    assert!(!document.dirty());

    document.set_selection(vec![1]).unwrap();
    document
        .apply(Box::new(DeleteSites::new(vec![1])), Origin::User)
        .unwrap();
    assert_eq!(document.structure().sites().len(), 1);
    assert!(document.selection().is_empty());
    assert!(document.undo());
    assert_eq!(document.structure().sites().len(), 2);
}

#[test]
fn command_trait_reports_descriptions() {
    let command = move_first([0.2, 0.2, 0.2]);
    assert_eq!(command.description(), "Move site");
    let macro_command: Box<dyn Command> = Box::new(MacroCommand::new("Grouped"));
    assert_eq!(macro_command.description(), "Grouped");
}
