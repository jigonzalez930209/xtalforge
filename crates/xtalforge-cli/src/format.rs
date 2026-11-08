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
}

impl StructureFormat {
    /// Format name understood by `xtalforge_io`.
    pub(crate) const fn name(self) -> &'static str {
        match self {
            Self::Cif => "cif",
            Self::Poscar => "poscar",
            Self::Xyz => "xyz",
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
            "xyz" => Some(Self::Xyz),
            _ => None,
        }
    }
}
