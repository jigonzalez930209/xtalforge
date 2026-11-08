//! LAMMPS data reader and writer (Task 3.3.4).
//!
//! Ported from Crystal-Builder `xtal/io/lammps.py` (MIT). The format covers
//! `atom_style full`: the restricted triclinic box, its masses, the atoms with
//! molecule id and charge, and the bonds. Bonds are read and written as
//! explicit-bond records; nothing is perceived here.

#![allow(clippy::indexing_slicing)] // fixed-format parser, lengths checked above
mod read;
mod write;

pub use read::read_lammps;
pub use write::write_lammps;
