//! Structure formats accepted by `convert`, `supercell` and `slab`.

use clap::ValueEnum;

/// A structure format the CLI can read and write.
#[derive(Clone, Copy, Debug, PartialEq, Eq, ValueEnum)]
pub(crate) enum StructureFormat {
    /// Crystallographic Information File.
    Cif,
    /// Magnetic CIF.
    Mcif,
    /// mmCIF/PDBx.
    Mmcif,
    /// VASP POSCAR/CONTCAR.
    Poscar,
    /// Plain or extended XYZ.
    Xyz,
    /// pymatgen Structure JSON.
    PmgJson,
    /// CSSR, as read by Zeo++.
    Cssr,
    /// DFTB+ GEN geometry.
    Gen,
    /// LAMMPS data file.
    Lammps,
    /// SHELX .res/.ins.
    Shelx,
    /// PDB.
    Pdb,
    /// Extended XYZ with arbitrary properties.
    Extxyz,
    /// Quantum ESPRESSO pw.x input or output.
    Pwscf,
    /// CASTEP .cell.
    Castep,
    /// FHI-aims geometry.in.
    FhiAims,
    /// CP2K input.
    Cp2k,
    /// ABINIT input.
    Abinit,
    /// CRYSTAL input.
    Crystal,
    /// Tripos MOL2.
    Mol2,
    /// SDF/MOL.
    Sdf,
    /// VASP OUTCAR.
    Outcar,
    /// VASP vasprun.xml.
    Vasprun,
    /// OPTIMADE structures entry.
    Optimade,
}

impl StructureFormat {
    /// Format name understood by `xtalforge_io`.
    pub(crate) const fn name(self) -> &'static str {
        match self {
            Self::Cif => "cif",
            Self::Mcif => "mcif",
            Self::Mmcif => "mmcif",
            Self::Poscar => "poscar",
            Self::Xyz => "xyz",
            Self::PmgJson => "pmg-json",
            Self::Cssr => "cssr",
            Self::Gen => "gen",
            Self::Lammps => "lammps-data",
            Self::Shelx => "shelx",
            Self::Pdb => "pdb",
            Self::Extxyz => "extxyz",
            Self::Pwscf => "pwscf",
            Self::Castep => "castep",
            Self::FhiAims => "fhi-aims",
            Self::Cp2k => "cp2k",
            Self::Abinit => "abinit",
            Self::Crystal => "crystal",
            Self::Mol2 => "mol2",
            Self::Sdf => "sdf",
            Self::Outcar => "outcar",
            Self::Vasprun => "vasprun",
            Self::Optimade => "optimade",
        }
    }

    /// Format inferred from a file extension.
    pub(crate) fn from_path(path: &str) -> Option<Self> {
        let extension = std::path::Path::new(path)
            .extension()
            .and_then(std::ffi::OsStr::to_str)?
            .to_ascii_lowercase();
        match extension.as_str() {
            "cif" => Some(Self::Cif),
            "mcif" => Some(Self::Mcif),
            "mmcif" | "pdbx" => Some(Self::Mmcif),
            "poscar" | "contcar" | "vasp" => Some(Self::Poscar),
            "xyz" => Some(Self::Xyz),
            "extxyz" => Some(Self::Extxyz),
            "json" => Some(Self::PmgJson),
            "cssr" => Some(Self::Cssr),
            "gen" => Some(Self::Gen),
            "data" | "lmp" => Some(Self::Lammps),
            "res" | "ins" => Some(Self::Shelx),
            "pdb" | "ent" => Some(Self::Pdb),
            "pwi" | "pwo" | "pwscf" => Some(Self::Pwscf),
            "cell" => Some(Self::Castep),
            "aims" | "in" => Some(Self::FhiAims),
            "cp2k" | "inp" => Some(Self::Cp2k),
            "abi" | "abinit" => Some(Self::Abinit),
            "crystal" | "d12" => Some(Self::Crystal),
            "mol2" => Some(Self::Mol2),
            "sdf" | "mol" | "mdl" => Some(Self::Sdf),
            "outcar" => Some(Self::Outcar),
            "vasprun" | "xml" => Some(Self::Vasprun),
            "optimade" => Some(Self::Optimade),
            _ => None,
        }
    }
}

/// A picture format the CLI can write.
#[derive(Clone, Copy, Debug, PartialEq, Eq, ValueEnum)]
pub(crate) enum PictureFormat {
    /// Vector document.
    Svg,
    /// Raster image.
    Png,
}

/// A mesh format the CLI can write.
#[derive(Clone, Copy, Debug, PartialEq, Eq, ValueEnum)]
pub(crate) enum MeshFormat {
    /// ASCII STL.
    Stl,
    /// Wavefront OBJ.
    Obj,
    /// glTF 2.0 JSON.
    Gltf,
}
