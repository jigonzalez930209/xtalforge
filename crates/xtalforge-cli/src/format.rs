//! Structure formats accepted by `convert`, `supercell` and `slab`.

use clap::ValueEnum;

/// A structure format the CLI can read and write.
#[derive(Clone, Copy, Debug, PartialEq, Eq, ValueEnum)]
pub(crate) enum StructureFormat {
    /// Crystallographic Information File.
    Cif,
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
    /// PDB, write-only.
    Pdb,
}

impl StructureFormat {
    /// Format name understood by `xtalforge_io`.
    pub(crate) const fn name(self) -> &'static str {
        match self {
            Self::Cif => "cif",
            Self::Poscar => "poscar",
            Self::Xyz => "xyz",
            Self::PmgJson => "pmg-json",
            Self::Cssr => "cssr",
            Self::Gen => "gen",
            Self::Lammps => "lammps-data",
            Self::Shelx => "shelx",
            Self::Pdb => "pdb",
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
            "poscar" | "contcar" | "vasp" => Some(Self::Poscar),
            "xyz" | "extxyz" => Some(Self::Xyz),
            "json" => Some(Self::PmgJson),
            "cssr" => Some(Self::Cssr),
            "gen" => Some(Self::Gen),
            "data" | "lmp" => Some(Self::Lammps),
            "res" | "ins" => Some(Self::Shelx),
            "pdb" => Some(Self::Pdb),
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
