//! SHELX `.res`/`.ins` reader and writer (Task 3.3.5).
//!
//! Ported from Crystal-Builder `xtal/io/shelx.py` (MIT). The reader takes the
//! cell, `LATT`, `SYMM`, `SFAC`, `FVAR` and the atom records, resolves SHELX
//! free-variable fields and keeps each atom's `PART` as its disorder group.
//! `SYMM` and `LATT` are recorded as metadata but not applied: the atoms are
//! read in P1. The writer emits a P1 file with no free variables.

#![allow(clippy::indexing_slicing)] // fixed-format parser, lengths checked above
mod expr;
mod read;
mod write;

pub use read::read_shelx;
pub use write::write_shelx;
