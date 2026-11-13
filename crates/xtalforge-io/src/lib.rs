//! File-format readers and writers for xtalforge.
//!
//! Phase 1 provides CIF, POSCAR and XYZ; Phase 3 adds the remaining formats,
//! the format registry, project files and trajectories. Phase 13 adds the
//! plane-wave code formats, mmCIF/PDB, MOL2/SDF, the VASP output readers and
//! arbitrary extended-XYZ properties. Every reader and writer reports failures
//! with a stable error code and never panics on user input.

pub mod abinit;
pub mod aims;
pub mod atomic;
pub mod castep;
pub mod cgd;
pub mod cif;
pub mod cp2k;
pub mod crystal;
pub mod cssr;
pub mod cube;
pub mod dispatch;
pub mod gencodec;
pub mod grid;
pub mod lammps;
pub mod mmcif;
pub mod mol2;
mod molecule;
pub mod optimade;
pub mod pdb;
pub mod pmg;
pub mod poscar;
pub mod project;
pub mod pwscf;
pub mod registry;
pub mod sdf;
pub mod shelx;
pub mod stream;
pub mod text;
pub mod trajectory;
pub mod vasp;
pub mod xy;
pub mod xyz;

pub use abinit::{read_abinit, write_abinit};
pub use aims::{read_aims, write_aims};
pub use castep::{read_castep, write_castep};
pub use cgd::{CgdEntry, CgdFile, CgdNode, read_cgd, write_cgd};
pub use cif::{
    CifBlock, CifDocument, CifLoop, read_cif, read_cif_all, read_cif_block, read_cif_stream,
    write_cif,
};
pub use cp2k::{read_cp2k, write_cp2k};
pub use crystal::{read_crystal, write_crystal};
pub use cssr::{read_cssr, write_cssr};
pub use cube::{read_cube, write_cube};
pub use dispatch::{read_structure, write_structure};
pub use gencodec::{read_gen, write_gen};
pub use grid::Grid3;
pub use lammps::{read_lammps, write_lammps};
pub use mmcif::{read_mmcif, read_mmcif_block};
pub use mol2::{read_mol2, write_mol2};
pub use optimade::write_optimade;
pub use pdb::{read_pdb, write_pdb};
pub use pmg::{read_pmg_json, write_pmg_json};
pub use poscar::{read_poscar, write_poscar};
pub use project::{FORMAT_VERSION, ProjectParts, read_project, write_project};
pub use pwscf::{read_pwscf, write_pwscf};
pub use sdf::{read_sdf, read_sdf_all, write_sdf};
pub use shelx::{read_shelx, write_shelx};
pub use stream::{LammpsStream, StreamFrame, XyzStream};
pub use text::{decode_text, read_text};
pub use trajectory::{
    Frame, Trajectory, XyzIndex, read_geo_end, read_lammps_dump, read_lammps_dump_with_types,
    read_trajectory, read_trajectory_xyz, read_xdatcar, write_trajectory_xyz,
};
pub use vasp::{read_outcar, read_outcar_structure, read_vasprun, read_vasprun_structure};
pub use xy::{PatternPoint, read_xy, write_xy};
pub use xyz::{read_xyz, write_xyz};
