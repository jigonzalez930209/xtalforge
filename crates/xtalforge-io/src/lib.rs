//! File-format readers and writers for xtalforge.
//!
//! Phase 1 provides CIF, POSCAR and XYZ; Phase 3 adds the remaining formats,
//! the format registry, project files and trajectories. Every reader and
//! writer reports failures with a stable error code and never panics on user
//! input.

pub mod cif;
pub mod dispatch;
pub mod poscar;
pub mod text;
pub mod xyz;

pub use cif::{CifBlock, CifDocument, CifLoop, read_cif, read_cif_all, read_cif_block, write_cif};
pub use dispatch::{read_structure, write_structure};
pub use poscar::{read_poscar, write_poscar};
pub use text::{decode_text, read_text};
pub use xyz::{read_xyz, write_xyz};
