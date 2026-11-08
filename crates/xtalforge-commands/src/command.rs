//! Command trait and undo stack (Task 3.9.1).

use std::fmt;

use serde::{Deserialize, Serialize};
use xtalforge_core::{Error, Structure};

/// Default maximum number of undo steps kept by [`CommandStack`].
pub const DEFAULT_MAX_DEPTH: usize = 200;

/// Who initiated a command.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default, Serialize, Deserialize)]
#[serde(rename_all = "snake_case")]
pub enum Origin {
    /// A person using a front end.
    #[default]
    User,
    /// A script or the CLI.
    Script,
    /// The agent layer.
    Agent,
}

/// One entry of the serializable command log.
#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
pub struct CommandRecord {
    /// Human-readable description of the command.
    pub description: String,
    /// Who initiated it.
    pub origin: Origin,
}

/// One undoable edit.
pub trait Command: fmt::Debug + Send {
    /// Applies the edit and returns the new structure.
    ///
    /// # Errors
    ///
    /// Returns a typed core error when the edit is invalid.
    fn apply(&self, structure: &Structure) -> Result<Structure, Error>;

    /// Human-readable description.
    fn description(&self) -> String;
}

/// The closure type used by [`SnapshotEdit`] and [`StructureOperation`].
pub type StructureFn = dyn Fn(&Structure) -> Result<Structure, Error> + Send + Sync;

/// A command that applies several commands in order.
#[derive(Debug)]
pub struct MacroCommand {
    description: String,
    commands: Vec<Box<dyn Command>>,
}

impl MacroCommand {
    /// Creates an empty macro with a description.
    #[must_use]
    pub fn new(description: impl Into<String>) -> Self {
        Self {
            description: description.into(),
            commands: Vec::new(),
        }
    }

    /// Creates a macro from already-built commands.
    #[must_use]
    pub fn from_commands(description: impl Into<String>, commands: Vec<Box<dyn Command>>) -> Self {
        Self {
            description: description.into(),
            commands,
        }
    }

    /// Appends one command.
    pub fn push(&mut self, command: Box<dyn Command>) {
        self.commands.push(command);
    }

    /// The commands in application order.
    #[must_use]
    pub fn commands(&self) -> &[Box<dyn Command>] {
        &self.commands
    }

    /// Number of commands in the macro.
    #[must_use]
    pub fn len(&self) -> usize {
        self.commands.len()
    }

    /// Whether the macro has no commands.
    #[must_use]
    pub fn is_empty(&self) -> bool {
        self.commands.is_empty()
    }
}

impl Command for MacroCommand {
    fn apply(&self, structure: &Structure) -> Result<Structure, Error> {
        let mut current = structure.clone();
        for command in &self.commands {
            current = command.apply(&current)?;
        }
        Ok(current)
    }

    fn description(&self) -> String {
        self.description.clone()
    }
}

/// A command that swaps the whole structure for a stored replacement.
#[derive(Debug, Clone)]
pub struct ReplaceStructure {
    replacement: Structure,
    description: String,
}

impl ReplaceStructure {
    /// Creates a replacement with the default description `Replace structure`.
    #[must_use]
    pub fn new(replacement: Structure) -> Self {
        Self {
            replacement,
            description: "Replace structure".to_owned(),
        }
    }

    /// Creates a replacement with an explicit description.
    #[must_use]
    pub fn with_description(replacement: Structure, description: impl Into<String>) -> Self {
        Self {
            replacement,
            description: description.into(),
        }
    }

    /// The structure this command installs.
    #[must_use]
    pub fn structure(&self) -> &Structure {
        &self.replacement
    }
}

impl Command for ReplaceStructure {
    fn apply(&self, _structure: &Structure) -> Result<Structure, Error> {
        Ok(self.replacement.clone())
    }

    fn description(&self) -> String {
        self.description.clone()
    }
}

/// An arbitrary edit made undoable by keeping before and after snapshots.
pub struct SnapshotEdit {
    description: String,
    edit: Box<StructureFn>,
}

impl SnapshotEdit {
    /// Creates an edit from a description and a pure structure function.
    #[must_use]
    pub fn new<F>(description: impl Into<String>, edit: F) -> Self
    where
        F: Fn(&Structure) -> Result<Structure, Error> + Send + Sync + 'static,
    {
        Self {
            description: description.into(),
            edit: Box::new(edit),
        }
    }
}

impl fmt::Debug for SnapshotEdit {
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        formatter
            .debug_struct("SnapshotEdit")
            .field("description", &self.description)
            .finish_non_exhaustive()
    }
}

impl Command for SnapshotEdit {
    fn apply(&self, structure: &Structure) -> Result<Structure, Error> {
        (self.edit)(structure)
    }

    fn description(&self) -> String {
        self.description.clone()
    }
}

/// A whole-structure rebuild that can be previewed before it is pushed.
pub struct StructureOperation {
    description: String,
    operation: Box<StructureFn>,
}

impl StructureOperation {
    /// Creates an operation from a description and a pure rebuild function.
    #[must_use]
    pub fn new<F>(description: impl Into<String>, operation: F) -> Self
    where
        F: Fn(&Structure) -> Result<Structure, Error> + Send + Sync + 'static,
    {
        Self {
            description: description.into(),
            operation: Box::new(operation),
        }
    }

    /// Computes the result without changing the input.
    ///
    /// # Errors
    ///
    /// Returns the error of the wrapped operation.
    pub fn preview(&self, structure: &Structure) -> Result<Structure, Error> {
        (self.operation)(structure)
    }
}

impl fmt::Debug for StructureOperation {
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        formatter
            .debug_struct("StructureOperation")
            .field("description", &self.description)
            .finish_non_exhaustive()
    }
}

impl Command for StructureOperation {
    fn apply(&self, structure: &Structure) -> Result<Structure, Error> {
        self.preview(structure)
    }

    fn description(&self) -> String {
        self.description.clone()
    }
}

#[derive(Debug)]
struct Entry {
    before: Structure,
    after: Structure,
    before_rev: u64,
    after_rev: u64,
    record: CommandRecord,
}

/// Bounded undo/redo history over structures.
#[derive(Debug)]
pub struct CommandStack {
    done: Vec<Entry>,
    undone: Vec<Entry>,
    log: Vec<CommandRecord>,
    clean_rev: Option<u64>,
    current_rev: u64,
    next_rev: u64,
    max_depth: usize,
    trimmed: usize,
}

impl CommandStack {
    /// Creates an empty stack with the default depth of 200.
    #[must_use]
    pub fn new() -> Self {
        Self::with_max_depth(DEFAULT_MAX_DEPTH)
    }

    /// Creates an empty stack with an explicit maximum depth.
    ///
    /// Values below one are raised to one, so a stack always keeps its most
    /// recent edit.
    #[must_use]
    pub fn with_max_depth(max_depth: usize) -> Self {
        Self {
            done: Vec::new(),
            undone: Vec::new(),
            log: Vec::new(),
            clean_rev: Some(0),
            current_rev: 0,
            next_rev: 1,
            max_depth: max_depth.max(1),
            trimmed: 0,
        }
    }

    /// Applies a command, records it and clears the redo branch.
    ///
    /// # Errors
    ///
    /// Returns the command error and leaves the structure and the history
    /// untouched when the command is invalid.
    pub fn push(
        &mut self,
        structure: &mut Structure,
        command: Box<dyn Command>,
        origin: Origin,
    ) -> Result<(), Error> {
        let before = structure.clone();
        let after = command.apply(&before)?;
        let record = CommandRecord {
            description: command.description(),
            origin,
        };
        *structure = after.clone();
        let before_rev = self.current_rev;
        let after_rev = self.next_rev;
        self.next_rev += 1;
        self.current_rev = after_rev;
        self.done.push(Entry {
            before,
            after,
            before_rev,
            after_rev,
            record: record.clone(),
        });
        self.log.push(record);
        self.undone.clear();
        self.trimmed = 0;
        self.trim_to_depth();
        Ok(())
    }

    /// Undoes the most recent step.
    ///
    /// Returns `true` when a step was undone.
    pub fn undo(&mut self, structure: &mut Structure) -> bool {
        let Some(entry) = self.done.pop() else {
            return false;
        };
        *structure = entry.before.clone();
        self.current_rev = entry.before_rev;
        self.log.pop();
        self.undone.push(entry);
        true
    }

    /// Redoes the most recently undone step.
    ///
    /// Returns `true` when a step was redone.
    pub fn redo(&mut self, structure: &mut Structure) -> bool {
        let Some(entry) = self.undone.pop() else {
            return false;
        };
        *structure = entry.after.clone();
        self.current_rev = entry.after_rev;
        self.log.push(entry.record.clone());
        self.done.push(entry);
        true
    }

    /// Drops every done and undone step.
    pub fn clear(&mut self) {
        self.done.clear();
        self.undone.clear();
        self.log.clear();
        self.trimmed = 0;
    }

    /// Discards the redo branch without touching the done history.
    pub fn discard_redo(&mut self) {
        self.undone.clear();
    }

    /// Whether at least one step can be undone.
    #[must_use]
    pub fn can_undo(&self) -> bool {
        !self.done.is_empty()
    }

    /// Whether at least one step can be redone.
    #[must_use]
    pub fn can_redo(&self) -> bool {
        !self.undone.is_empty()
    }

    /// Number of applied steps.
    #[must_use]
    pub fn undo_depth(&self) -> usize {
        self.done.len()
    }

    /// Number of undone steps available for redo.
    #[must_use]
    pub fn redo_depth(&self) -> usize {
        self.undone.len()
    }

    /// Whether the current position differs from the last clean mark.
    #[must_use]
    pub fn is_dirty(&self) -> bool {
        self.clean_rev != Some(self.current_rev)
    }

    /// Marks the current position as clean.
    pub fn mark_clean(&mut self) {
        self.clean_rev = Some(self.current_rev);
    }

    /// Maximum number of retained steps.
    #[must_use]
    pub fn max_depth(&self) -> usize {
        self.max_depth
    }

    /// Changes the maximum depth and drops the oldest steps if needed.
    pub fn set_max_depth(&mut self, max_depth: usize) {
        self.max_depth = max_depth.max(1);
        self.trimmed = 0;
        self.trim_to_depth();
    }

    /// Serialized log of the applied command records, oldest first.
    #[must_use]
    pub fn log(&self) -> &[CommandRecord] {
        &self.log
    }

    /// Number of oldest steps dropped by the most recent push or resize.
    #[must_use]
    pub fn trimmed(&self) -> usize {
        self.trimmed
    }

    fn trim_to_depth(&mut self) {
        while self.done.len() > self.max_depth {
            self.done.remove(0);
            self.log.remove(0);
            self.trimmed += 1;
        }
    }
}

impl Default for CommandStack {
    fn default() -> Self {
        Self::new()
    }
}
