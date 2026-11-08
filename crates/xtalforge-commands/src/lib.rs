//! Undoable commands and document state for xtalforge (Phase 3).
//!
//! Every mutation goes through a [`Command`](command::Command); the
//! [`Document`](document::Document) owns the structure, the undo stack, the
//! selection and the atom groups.

pub mod atoms;
pub mod bonds;
pub mod cell;
pub mod clipboard;
pub mod command;
pub mod document;
mod edit;
pub mod groups;
pub mod later;
pub mod symmetry;

pub use atoms::{AddSites, DeleteSites, MoveSites, SetElement, SetOccupancy};
pub use bonds::{AddBond, RemoveBond, SetBondRules};
pub use cell::{PositionPolicy, SetLattice, ShiftOrigin, Supercell, WrapIntoCell};
pub use clipboard::{Fragment, PasteFragment, copy_selection};
pub use command::{
    Command, CommandRecord, CommandStack, MacroCommand, Origin, ReplaceStructure, SnapshotEdit,
    StructureOperation,
};
pub use document::Document;
pub use groups::{AtomGroup, GroupStore};
pub use later::{Interpenetrate, MarkConnectionPoints, SetCharges};
