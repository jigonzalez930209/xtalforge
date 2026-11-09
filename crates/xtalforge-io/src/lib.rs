//! File-format readers and writers for xtalforge.
//!
//! Phase 1 provides CIF, POSCAR and XYZ; Phase 3 adds the remaining formats,
//! the format registry, project files and trajectories. Every reader and
//! writer reports failures with a stable error code and never panics on user
//! input.

pub mod atomic;
pub mod cgd;
pub mod cif;
pub mod cssr;
pub mod cube;
pub mod dispatch;
pub mod gencodec;
pub mod grid;
pub mod lammps;
pub mod pdb;
pub mod pmg;
pub mod poscar;
pub mod project;
pub mod registry;
pub mod shelx;
pub mod text;
pub mod trajectory;
pub mod xy;
pub mod xyz;

pub use cgd::{CgdEntry, CgdFile, CgdNode, read_cgd, write_cgd};
pub use cif::{CifBlock, CifDocument, CifLoop, read_cif, read_cif_all, read_cif_block, write_cif};
pub use cssr::{read_cssr, write_cssr};
pub use cube::{read_cube, write_cube};
pub use dispatch::{read_structure, write_structure};
pub use gencodec::{read_gen, write_gen};
pub use grid::Grid3;
pub use lammps::{read_lammps, write_lammps};
pub use pdb::write_pdb;
pub use pmg::{read_pmg_json, write_pmg_json};
pub use poscar::{read_poscar, write_poscar};
pub use project::{FORMAT_VERSION, ProjectParts, read_project, write_project};
pub use shelx::{read_shelx, write_shelx};
pub use text::{decode_text, read_text};
pub use trajectory::{
    Frame, Trajectory, XyzIndex, read_trajectory, read_trajectory_xyz, read_xdatcar,
    write_trajectory_xyz,
};
pub use xy::{PatternPoint, read_xy, write_xy};
pub use xyz::{read_xyz, write_xyz};
