//! The built-in format list (Task 3.1, extended by Task 13.2).

use super::{Capabilities, FormatSpec};
use crate::{
    read_abinit, read_aims, read_castep, read_cif, read_cp2k, read_crystal, read_cssr, read_gen,
    read_lammps, read_mmcif, read_mol2, read_outcar_structure, read_pdb, read_pmg_json,
    read_poscar, read_pwscf, read_sdf, read_shelx, read_vasprun_structure, read_xyz, write_abinit,
    write_aims, write_castep, write_cif, write_cp2k, write_crystal, write_cssr, write_gen,
    write_lammps, write_mol2, write_optimade, write_pdb, write_pmg_json, write_poscar, write_pwscf,
    write_sdf, write_shelx, write_xyz,
};

const EMPTY: Capabilities = Capabilities {
    cell: false,
    symmetry: false,
    occupancy: false,
    anisotropic_adp: false,
    charges: false,
    bonds: false,
    labels: false,
    frames: false,
};

/// The built-in registry.
pub const FORMATS: &[FormatSpec] = &[
    FormatSpec {
        name: "cif",
        extensions: &[".cif", ".mcif"],
        read: Some(read_cif),
        write: Some(write_cif),
        capabilities: Capabilities {
            cell: true,
            symmetry: true,
            occupancy: true,
            anisotropic_adp: true,
            charges: true,
            bonds: false,
            labels: true,
            frames: false,
        },
    },
    FormatSpec {
        name: "mmcif",
        extensions: &[".mmcif", ".pdbx"],
        read: Some(read_mmcif),
        write: None,
        capabilities: Capabilities {
            cell: true,
            symmetry: true,
            occupancy: true,
            labels: true,
            ..EMPTY
        },
    },
    FormatSpec {
        name: "poscar",
        extensions: &[".poscar", ".vasp", ".contcar"],
        read: Some(read_poscar),
        write: Some(write_poscar),
        capabilities: Capabilities {
            cell: true,
            labels: true,
            ..EMPTY
        },
    },
    FormatSpec {
        name: "pmg-json",
        extensions: &[".json"],
        read: Some(read_pmg_json),
        write: Some(write_pmg_json),
        capabilities: Capabilities {
            cell: true,
            occupancy: true,
            labels: true,
            ..EMPTY
        },
    },
    FormatSpec {
        name: "cssr",
        extensions: &[".cssr"],
        read: Some(read_cssr),
        write: Some(write_cssr),
        capabilities: Capabilities {
            cell: true,
            ..EMPTY
        },
    },
    FormatSpec {
        name: "gen",
        extensions: &[".gen"],
        read: Some(read_gen),
        write: Some(write_gen),
        capabilities: Capabilities {
            cell: true,
            labels: true,
            ..EMPTY
        },
    },
    FormatSpec {
        name: "lammps-data",
        extensions: &[".data", ".lmp"],
        read: Some(read_lammps),
        write: Some(write_lammps),
        capabilities: Capabilities {
            cell: true,
            charges: true,
            bonds: true,
            ..EMPTY
        },
    },
    FormatSpec {
        name: "shelx",
        extensions: &[".res", ".ins"],
        read: Some(read_shelx),
        write: Some(write_shelx),
        capabilities: Capabilities {
            cell: true,
            occupancy: true,
            labels: true,
            ..EMPTY
        },
    },
    FormatSpec {
        name: "pdb",
        extensions: &[".pdb", ".ent"],
        read: Some(read_pdb),
        write: Some(write_pdb),
        capabilities: Capabilities {
            cell: true,
            symmetry: true,
            occupancy: true,
            bonds: true,
            labels: true,
            ..EMPTY
        },
    },
    FormatSpec {
        name: "xyz",
        extensions: &[".xyz"],
        read: Some(read_xyz),
        write: Some(write_xyz),
        capabilities: Capabilities {
            cell: true,
            ..EMPTY
        },
    },
    FormatSpec {
        name: "extxyz",
        extensions: &[".extxyz"],
        read: Some(read_xyz),
        write: Some(write_xyz),
        capabilities: Capabilities {
            cell: true,
            charges: true,
            labels: true,
            ..EMPTY
        },
    },
    FormatSpec {
        name: "pwscf",
        extensions: &[".pwi", ".pwo", ".pwscf"],
        read: Some(read_pwscf),
        write: Some(write_pwscf),
        capabilities: Capabilities {
            cell: true,
            ..EMPTY
        },
    },
    FormatSpec {
        name: "castep",
        extensions: &[".cell"],
        read: Some(read_castep),
        write: Some(write_castep),
        capabilities: Capabilities {
            cell: true,
            occupancy: true,
            ..EMPTY
        },
    },
    FormatSpec {
        name: "fhi-aims",
        extensions: &[".geometry.in", ".aims"],
        read: Some(read_aims),
        write: Some(write_aims),
        capabilities: Capabilities {
            cell: true,
            ..EMPTY
        },
    },
    FormatSpec {
        name: "cp2k",
        extensions: &[".cp2k", ".inp"],
        read: Some(read_cp2k),
        write: Some(write_cp2k),
        capabilities: Capabilities {
            cell: true,
            ..EMPTY
        },
    },
    FormatSpec {
        name: "abinit",
        extensions: &[".abi", ".abinit"],
        read: Some(read_abinit),
        write: Some(write_abinit),
        capabilities: Capabilities {
            cell: true,
            ..EMPTY
        },
    },
    FormatSpec {
        name: "crystal",
        extensions: &[".crystal", ".d12"],
        read: Some(read_crystal),
        write: Some(write_crystal),
        capabilities: Capabilities {
            cell: true,
            symmetry: true,
            ..EMPTY
        },
    },
    FormatSpec {
        name: "mol2",
        extensions: &[".mol2"],
        read: Some(read_mol2),
        write: Some(write_mol2),
        capabilities: Capabilities {
            charges: true,
            bonds: true,
            labels: true,
            ..EMPTY
        },
    },
    FormatSpec {
        name: "sdf",
        extensions: &[".sdf", ".mol", ".mdl"],
        read: Some(read_sdf),
        write: Some(write_sdf),
        capabilities: Capabilities {
            charges: true,
            bonds: true,
            ..EMPTY
        },
    },
    FormatSpec {
        name: "optimade",
        extensions: &[".optimade"],
        read: None,
        write: Some(write_optimade),
        capabilities: Capabilities {
            cell: true,
            labels: true,
            ..EMPTY
        },
    },
    FormatSpec {
        name: "outcar",
        extensions: &[".outcar"],
        read: Some(read_outcar_structure),
        write: None,
        capabilities: Capabilities {
            cell: true,
            frames: true,
            ..EMPTY
        },
    },
    FormatSpec {
        name: "vasprun",
        extensions: &[".xml", ".vasprun.xml"],
        read: Some(read_vasprun_structure),
        write: None,
        capabilities: Capabilities {
            cell: true,
            frames: true,
            ..EMPTY
        },
    },
];
