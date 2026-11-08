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
pub mod symmetry;
