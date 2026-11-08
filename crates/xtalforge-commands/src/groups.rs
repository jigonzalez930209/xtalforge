//! Atom groups: named display sets stored with the document (Task 3.10.5).

use serde::{Deserialize, Serialize};
use xtalforge_core::{Error, Structure};

use crate::command::{DEFAULT_MAX_DEPTH, Origin};
use crate::edit::check_index;

/// A named set of asymmetric-unit sites with display overrides.
#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
pub struct AtomGroup {
    /// Display name.
    pub name: String,
    /// Asymmetric-unit site indexes.
    pub atoms: Vec<usize>,
    /// Optional RGB colour; `None` uses the element colours.
    pub color: Option<[u8; 3]>,
    /// Whether the group's atoms are hidden.
    pub hidden: bool,
}

impl AtomGroup {
    /// Creates a visible group with element colours.
    ///
    /// # Errors
    ///
    /// Returns [`Error::InvalidInput`] for a blank name.
    pub fn new(name: impl Into<String>, atoms: Vec<usize>) -> Result<Self, Error> {
        let name = name.into();
        validate_name(&name)?;
        Ok(Self {
            name,
            atoms,
            color: None,
            hidden: false,
        })
    }
}

/// The document's atom groups.
#[derive(Debug, Clone, Default, PartialEq, Eq, Serialize, Deserialize)]
pub struct GroupStore {
    groups: Vec<AtomGroup>,
}

impl GroupStore {
    /// Creates an empty store.
    #[must_use]
    pub fn new() -> Self {
        Self::default()
    }

    /// Number of groups.
    #[must_use]
    pub fn len(&self) -> usize {
        self.groups.len()
    }

    /// Whether there are no groups.
    #[must_use]
    pub fn is_empty(&self) -> bool {
        self.groups.is_empty()
    }

    /// All groups in insertion order.
    #[must_use]
    pub fn groups(&self) -> &[AtomGroup] {
        &self.groups
    }

    /// One group by row.
    #[must_use]
    pub fn get(&self, index: usize) -> Option<&AtomGroup> {
        self.groups.get(index)
    }

    /// Appends a group.
    ///
    /// # Errors
    ///
    /// Returns [`Error::InvalidInput`] for a blank name.
    pub fn add(&mut self, group: AtomGroup) -> Result<(), Error> {
        validate_name(&group.name)?;
        self.groups.push(group);
        Ok(())
    }

    /// Removes a group by row.
    ///
    /// # Errors
    ///
    /// Returns [`Error::InvalidInput`] when the row is out of range.
    pub fn remove(&mut self, index: usize) -> Result<(), Error> {
        check_index(index, self.groups.len(), "group")?;
        self.groups.remove(index);
        Ok(())
    }

    /// Recolours a group.
    ///
    /// # Errors
    ///
    /// Returns [`Error::InvalidInput`] when the row is out of range.
    pub fn set_color(&mut self, index: usize, color: Option<[u8; 3]>) -> Result<(), Error> {
        check_index(index, self.groups.len(), "group")?;
        if let Some(group) = self.groups.get_mut(index) {
            group.color = color;
        }
        Ok(())
    }

    /// Shows or hides a group.
    ///
    /// # Errors
    ///
    /// Returns [`Error::InvalidInput`] when the row is out of range.
    pub fn set_hidden(&mut self, index: usize, hidden: bool) -> Result<(), Error> {
        check_index(index, self.groups.len(), "group")?;
        if let Some(group) = self.groups.get_mut(index) {
            group.hidden = hidden;
        }
        Ok(())
    }
}

/// A document-state command over [`GroupStore`].
pub trait GroupCommand: std::fmt::Debug + Send {
    /// Applies the edit to the group store.
    ///
    /// # Errors
    ///
    /// Returns a typed error when the edit is invalid.
    fn apply(&self, groups: &mut GroupStore) -> Result<(), Error>;

    /// Human-readable description.
    fn description(&self) -> String;

    /// Checks the command against the current structure.
    ///
    /// # Errors
    ///
    /// Returns a typed error when an atom index is out of range.
    fn validate(&self, _structure: &Structure) -> Result<(), Error> {
        Ok(())
    }
}

/// Adds an atom group.
#[derive(Debug, Clone)]
pub struct AddAtomGroup {
    group: AtomGroup,
}

impl AddAtomGroup {
    /// Creates the command.
    #[must_use]
    pub fn new(group: AtomGroup) -> Self {
        Self { group }
    }

    /// The group that will be added.
    #[must_use]
    pub fn group(&self) -> &AtomGroup {
        &self.group
    }
}

impl GroupCommand for AddAtomGroup {
    fn apply(&self, groups: &mut GroupStore) -> Result<(), Error> {
        let mut unique = self.group.atoms.clone();
        unique.sort_unstable();
        unique.dedup();
        if unique.len() != self.group.atoms.len() {
            return Err(Error::InvalidInput(
                "an atom group must not repeat an atom index".to_owned(),
            ));
        }
        groups.add(self.group.clone())
    }

    fn description(&self) -> String {
        format!("Add atom group {}", self.group.name)
    }

    fn validate(&self, structure: &Structure) -> Result<(), Error> {
        validate_name(&self.group.name)?;
        for index in &self.group.atoms {
            check_index(*index, structure.sites().len(), "site")?;
        }
        Ok(())
    }
}

/// Removes an atom group by row.
#[derive(Debug, Clone)]
pub struct RemoveAtomGroup {
    index: usize,
}

impl RemoveAtomGroup {
    /// Creates the command.
    #[must_use]
    pub fn new(index: usize) -> Self {
        Self { index }
    }

    /// Row the command removes.
    #[must_use]
    pub fn index(&self) -> usize {
        self.index
    }
}

impl GroupCommand for RemoveAtomGroup {
    fn apply(&self, groups: &mut GroupStore) -> Result<(), Error> {
        groups.remove(self.index)
    }

    fn description(&self) -> String {
        "Remove atom group".to_owned()
    }
}

/// Sets or clears the display colour of an atom group.
#[derive(Debug, Clone)]
pub struct SetGroupColor {
    index: usize,
    color: Option<[u8; 3]>,
}

impl SetGroupColor {
    /// Creates the command.
    #[must_use]
    pub fn new(index: usize, color: Option<[u8; 3]>) -> Self {
        Self { index, color }
    }

    /// Row the command recolours.
    #[must_use]
    pub fn index(&self) -> usize {
        self.index
    }
}

impl GroupCommand for SetGroupColor {
    fn apply(&self, groups: &mut GroupStore) -> Result<(), Error> {
        groups.set_color(self.index, self.color)
    }

    fn description(&self) -> String {
        "Colour atom group".to_owned()
    }
}

/// Shows or hides an atom group.
#[derive(Debug, Clone)]
pub struct SetGroupHidden {
    index: usize,
    hidden: bool,
}

impl SetGroupHidden {
    /// Creates the command.
    #[must_use]
    pub fn new(index: usize, hidden: bool) -> Self {
        Self { index, hidden }
    }

    /// Row the command toggles.
    #[must_use]
    pub fn index(&self) -> usize {
        self.index
    }
}

impl GroupCommand for SetGroupHidden {
    fn apply(&self, groups: &mut GroupStore) -> Result<(), Error> {
        groups.set_hidden(self.index, self.hidden)
    }

    fn description(&self) -> String {
        if self.hidden {
            "Hide atom group".to_owned()
        } else {
            "Show atom group".to_owned()
        }
    }
}

fn validate_name(name: &str) -> Result<(), Error> {
    if name.trim().is_empty() {
        return Err(Error::InvalidInput(
            "an atom group needs a non-empty name".to_owned(),
        ));
    }
    Ok(())
}

#[derive(Debug)]
struct GroupEntry {
    before: GroupStore,
    after: GroupStore,
    before_rev: u64,
    after_rev: u64,
}

/// Undo history for document-owned atom groups.
#[derive(Debug)]
pub(crate) struct GroupHistory {
    done: Vec<GroupEntry>,
    undone: Vec<GroupEntry>,
    clean_rev: Option<u64>,
    current_rev: u64,
    next_rev: u64,
    max_depth: usize,
    trimmed: usize,
}

impl GroupHistory {
    /// Applies a group command and records it.
    pub(crate) fn push(
        &mut self,
        groups: &mut GroupStore,
        command: Box<dyn GroupCommand>,
        _origin: Origin,
    ) -> Result<(), Error> {
        let before = groups.clone();
        let mut after = groups.clone();
        command.apply(&mut after)?;
        *groups = after.clone();
        let before_rev = self.current_rev;
        let after_rev = self.next_rev;
        self.next_rev += 1;
        self.current_rev = after_rev;
        self.done.push(GroupEntry {
            before,
            after,
            before_rev,
            after_rev,
        });
        self.undone.clear();
        self.trimmed = 0;
        self.trim_to_depth();
        Ok(())
    }

    pub(crate) fn undo(&mut self, groups: &mut GroupStore) -> bool {
        let Some(entry) = self.done.pop() else {
            return false;
        };
        *groups = entry.before.clone();
        self.current_rev = entry.before_rev;
        self.undone.push(entry);
        true
    }

    pub(crate) fn redo(&mut self, groups: &mut GroupStore) -> bool {
        let Some(entry) = self.undone.pop() else {
            return false;
        };
        *groups = entry.after.clone();
        self.current_rev = entry.after_rev;
        self.done.push(entry);
        true
    }

    pub(crate) fn discard_redo(&mut self) {
        self.undone.clear();
    }

    pub(crate) fn is_dirty(&self) -> bool {
        self.clean_rev != Some(self.current_rev)
    }

    pub(crate) fn mark_clean(&mut self) {
        self.clean_rev = Some(self.current_rev);
    }

    pub(crate) fn trimmed(&self) -> usize {
        self.trimmed
    }

    fn trim_to_depth(&mut self) {
        while self.done.len() > self.max_depth {
            self.done.remove(0);
            self.trimmed += 1;
        }
    }
}

impl Default for GroupHistory {
    fn default() -> Self {
        Self {
            done: Vec::new(),
            undone: Vec::new(),
            clean_rev: Some(0),
            current_rev: 0,
            next_rev: 1,
            max_depth: DEFAULT_MAX_DEPTH,
            trimmed: 0,
        }
    }
}
