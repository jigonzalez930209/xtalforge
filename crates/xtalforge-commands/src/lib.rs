//! Undoable commands and document state for xtalforge (Phase 3).
//!
//! Every mutation goes through a [`Command`]; the [`Document`] owns the
//! structure, the undo stack, the selection and the atom groups.

pub mod atoms;
pub mod bonds;
pub mod cell;
pub mod clipboard;
pub mod command;
pub mod document;
mod edit;
pub mod ff;
pub mod groups;
pub mod later;
pub mod symmetry;

pub use atoms::{
    AddSites, DeleteSites, MoveSites, ReplaceSpecies, SetElement, SetOccupancy, SortSites,
};
pub use bonds::{AddBond, RemoveBond, SetBondRules};
pub use cell::{
    MakeSlab, PositionPolicy, SetLattice, ShiftOrigin, Supercell, TransformCell, WrapIntoCell,
};
pub use clipboard::{Fragment, InsertMolecules, PasteFragment, copy_selection};
pub use command::{
    Command, CommandRecord, CommandStack, MacroCommand, Origin, ReplaceStructure, SnapshotEdit,
    StructureOperation,
};
pub use document::Document;
pub use ff::{
    ATOM_TYPE_KEY, AddHydrogens, ApplyOptimizedGeometry, Prepare, PrepareStep, SetAtomTypes,
};
pub use groups::{AtomGroup, GroupStore};
pub use later::{Interpenetrate, MarkConnectionPoints, SetCharges};
pub use symmetry::{AssignWyckoff, DescendToSubgroup, Invert, MergeDuplicates};
