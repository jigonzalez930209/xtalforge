//! Document state: structure, history, selection and atom groups (Task 3.9.2).

use xtalforge_core::{Error, Structure};

use crate::command::{Command, CommandStack, Origin};
use crate::edit::check_index;
use crate::groups::{AtomGroup, GroupCommand, GroupHistory, GroupStore};

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum EditKind {
    Structure,
    Groups,
}

/// The editable state of one structure.
#[derive(Debug)]
pub struct Document {
    structure: Structure,
    stack: CommandStack,
    groups: GroupStore,
    group_history: GroupHistory,
    timeline: Vec<EditKind>,
    redo_timeline: Vec<EditKind>,
    selection: Vec<usize>,
    origin: Origin,
    name: String,
}

impl Document {
    /// Creates a clean document around one structure.
    #[must_use]
    pub fn new(structure: Structure) -> Self {
        let name = structure.name().to_owned();
        Self {
            structure,
            stack: CommandStack::new(),
            groups: GroupStore::new(),
            group_history: GroupHistory::default(),
            timeline: Vec::new(),
            redo_timeline: Vec::new(),
            selection: Vec::new(),
            origin: Origin::User,
            name,
        }
    }

    /// Applies a structure command and records it in the undo history.
    ///
    /// # Errors
    ///
    /// Returns the command error and leaves the document untouched when the
    /// edit is invalid.
    pub fn apply(&mut self, command: Box<dyn Command>, origin: Origin) -> Result<(), Error> {
        self.stack.push(&mut self.structure, command, origin)?;
        self.group_history.discard_redo();
        self.redo_timeline.clear();
        self.drop_timeline(EditKind::Structure, self.stack.trimmed());
        self.timeline.push(EditKind::Structure);
        self.origin = origin;
        self.prune_selection();
        Ok(())
    }

    /// Applies a document-state command over the atom groups.
    ///
    /// # Errors
    ///
    /// Returns the command error and leaves the document untouched when the
    /// edit is invalid.
    pub fn apply_group(
        &mut self,
        command: Box<dyn GroupCommand>,
        origin: Origin,
    ) -> Result<(), Error> {
        command.validate(&self.structure)?;
        self.group_history.push(&mut self.groups, command, origin)?;
        self.stack.discard_redo();
        self.redo_timeline.clear();
        self.drop_timeline(EditKind::Groups, self.group_history.trimmed());
        self.timeline.push(EditKind::Groups);
        self.origin = origin;
        Ok(())
    }

    /// Undoes the most recent edit of either kind.
    pub fn undo(&mut self) -> bool {
        let Some(kind) = self.timeline.pop() else {
            return false;
        };
        let undone = match kind {
            EditKind::Structure => self.stack.undo(&mut self.structure),
            EditKind::Groups => self.group_history.undo(&mut self.groups),
        };
        if !undone {
            self.timeline.push(kind);
            return false;
        }
        self.redo_timeline.push(kind);
        self.prune_selection();
        true
    }

    /// Redoes the most recently undone edit of either kind.
    pub fn redo(&mut self) -> bool {
        let Some(kind) = self.redo_timeline.pop() else {
            return false;
        };
        let redone = match kind {
            EditKind::Structure => self.stack.redo(&mut self.structure),
            EditKind::Groups => self.group_history.redo(&mut self.groups),
        };
        if !redone {
            self.redo_timeline.push(kind);
            return false;
        }
        self.timeline.push(kind);
        self.prune_selection();
        true
    }

    /// Whether there are unsaved edits.
    #[must_use]
    pub fn dirty(&self) -> bool {
        self.stack.is_dirty() || self.group_history.is_dirty()
    }

    /// Marks the current position as saved.
    pub fn mark_clean(&mut self) {
        self.stack.mark_clean();
        self.group_history.mark_clean();
    }

    /// Current structure.
    #[must_use]
    pub fn structure(&self) -> &Structure {
        &self.structure
    }

    /// Selected asymmetric-unit site indexes.
    #[must_use]
    pub fn selection(&self) -> &[usize] {
        &self.selection
    }

    /// Replaces the selection, validating every site index.
    ///
    /// Duplicates are dropped and order is preserved.
    ///
    /// # Errors
    ///
    /// Returns [`Error::InvalidInput`] when an index is out of range.
    pub fn set_selection(&mut self, selection: Vec<usize>) -> Result<(), Error> {
        for index in &selection {
            check_index(*index, self.structure.sites().len(), "site")?;
        }
        let mut cleaned = Vec::with_capacity(selection.len());
        for index in selection {
            if !cleaned.contains(&index) {
                cleaned.push(index);
            }
        }
        self.selection = cleaned;
        Ok(())
    }

    /// Document display name.
    #[must_use]
    pub fn name(&self) -> &str {
        &self.name
    }

    /// Sets the display name without touching the structure.
    pub fn set_name(&mut self, name: impl Into<String>) {
        self.name = name.into();
    }

    /// Origin of the most recent edit.
    #[must_use]
    pub fn origin(&self) -> Origin {
        self.origin
    }

    /// Atom groups in insertion order.
    #[must_use]
    pub fn groups(&self) -> &[AtomGroup] {
        self.groups.groups()
    }

    /// One atom group by row.
    #[must_use]
    pub fn group(&self, index: usize) -> Option<&AtomGroup> {
        self.groups.get(index)
    }

    /// Command log of the document.
    #[must_use]
    pub fn history(&self) -> &[crate::command::CommandRecord] {
        self.stack.log()
    }

    /// Whether at least one structure or group edit can be undone.
    #[must_use]
    pub fn can_undo(&self) -> bool {
        !self.timeline.is_empty()
    }

    /// Whether at least one structure or group edit can be redone.
    #[must_use]
    pub fn can_redo(&self) -> bool {
        !self.redo_timeline.is_empty()
    }

    fn prune_selection(&mut self) {
        let len = self.structure.sites().len();
        self.selection.retain(|index| *index < len);
    }

    fn drop_timeline(&mut self, kind: EditKind, count: usize) {
        for _ in 0..count {
            if let Some(position) = self.timeline.iter().position(|entry| *entry == kind) {
                self.timeline.remove(position);
            }
        }
    }
}
