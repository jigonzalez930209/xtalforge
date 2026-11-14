//! `xtalforge` builds, manipulates, analyses and exports crystal structures.
//!
//! The facade re-exports the domain crates under stable module names:
//!
//! - [`core`]: elements, lattices, sites, structures, P1, supercells, slabs,
//!   local environments, surfaces, edits, EOS and hulls.
//! - [`symmetry`]: space groups, Wyckoff positions, standard cells, reduction,
//!   subgroups, point groups, Brillouin zones, k-paths and the matcher.
//! - [`io`]: every registered reader and writer, project files and grids.
//! - [`commands`]: undoable edits and documents.
//! - [`view`]: scene data, colormaps and STL.
//! - [`ff`]: UFF/UFF4MOF, QEq/EQeq/EEQ charges, optimizers, vibrations,
//!   thermochemistry and elastic tensors.
//! - [`analysis`]: PXRD, porosity, nets, defects, orderings and interfaces.
//! - [`build`]: molecules, carbon, polymers and the SMILES embedder.
//! - [`mof`]: the MOF builder and the Phase 10 platform.
//! - [`engines`]: simulation input sets, output parsers and electronic structure.
//! - [`modules`] and [`agent`]: workspaces, jobs and the agent/MCP surface.
//!
//! # Features
//!
//! | Feature | Default | Effect |
//! |---|---|---|
//! | `parallel` | yes | Enables the Rayon paths in core, IO and symmetry. |
//! | `analysis` | yes | Re-exports `xtalforge-analysis` as [`analysis`]. |
//! | `build` | yes | Re-exports `xtalforge-build` as [`build`]; implies `analysis`. |
//! | `ff` | yes | Re-exports `xtalforge-ff` as [`ff`]; implies `build`. |
//! | `mof` | yes | Re-exports `xtalforge-mof` as [`mof`]; implies `ff`. |
//! | `engines` | yes | Re-exports `xtalforge-engines` as [`engines`]; implies `mof`. |
//! | `modules` | yes | Re-exports `xtalforge-modules` as [`modules`]. |
//! | `agent` | yes | Re-exports `xtalforge-agent` as [`agent`]; implies `mof` and `modules`. |
//! | `full` | yes | Shorthand for every surface tier above. |
//! | `sci-form` | no | Links the published `sci-form` crate for in-process GFN0/GFN1/GFN2. |
//! | `ml` | no | Reserved for the Python machine-learning bridge (Task 8.4). |
//! | `external` | no | Reserved for the subprocess engines (Task 4.9). |
//! | `net` | no | Blocking HTTP for OPTIMADE, COD and Materials Project fetches (Task 13.3). |
//! | `mof-data` | no | Loads the full PORMAKE 0.2.3 bundle (Task 17.1.3). |
//!
//! The base build always carries `core`, `io`, `symmetry`, `view` and
//! `commands`. Disable default features with `default-features = false` to
//! build only that base plus the tiers you name; a tier-free base is the
//! smallest build and keeps the same results for the modules it carries
//! (AGENTS.md rule 6).
//!
//! # Example
//!
//! ```
//! use xtalforge::core::Structure;
//! use xtalforge::symmetry::{SymmetryOptions, detect_symmetry};
//!
//! let cif = "data_Si\n_cell_length_a 5.43\n_cell_length_b 5.43\n_cell_length_c 5.43\n_cell_angle_alpha 90\n_cell_angle_beta 90\n_cell_angle_gamma 90\n_symmetry_space_group_name_H-M 'F d -3 m'\nloop_\n_atom_site_label\n_atom_site_type_symbol\n_atom_site_fract_x\n_atom_site_fract_y\n_atom_site_fract_z\nSi1 Si 0 0 0\n";
//! let structure = xtalforge::io::read_cif(cif).unwrap();
//! let dataset = detect_symmetry(&structure, SymmetryOptions::default()).unwrap();
//! assert_eq!(dataset.number, 227);
//! let supercell = structure.supercell(2, 2, 2).unwrap();
//! let poscar = xtalforge::io::write_structure(&supercell, "poscar").unwrap();
//! assert!(poscar.contains("Si"));
//! let _: Structure = supercell;
//! ```
#[cfg(feature = "agent")]
pub use xtalforge_agent as agent;
#[cfg(feature = "analysis")]
pub use xtalforge_analysis as analysis;
#[cfg(feature = "build")]
pub use xtalforge_build as build;
pub use xtalforge_commands as commands;
pub use xtalforge_core as core;
#[cfg(feature = "engines")]
pub use xtalforge_engines as engines;
#[cfg(feature = "ff")]
pub use xtalforge_ff as ff;
pub use xtalforge_io as io;
#[cfg(feature = "modules")]
pub use xtalforge_modules as modules;
#[cfg(feature = "mof")]
pub use xtalforge_mof as mof;
pub use xtalforge_symmetry as symmetry;
pub use xtalforge_view as view;

pub mod remote;

/// A structure read from disk or text together with its warnings.
#[derive(Debug, Clone, PartialEq)]
pub struct ReadResult {
    /// Parsed structure.
    pub structure: xtalforge_core::Structure,
    /// Structured reader warnings.
    pub warnings: Vec<xtalforge_core::StructureWarning>,
}

/// Reads a structure from a path, detecting the format from its extension or
/// content.
///
/// # Errors
///
/// Returns [`xtalforge_core::Error::InvalidInput`] when the file cannot be
/// read or the content is not a registered structure format.
pub fn read(path: impl AsRef<std::path::Path>) -> Result<ReadResult, xtalforge_core::Error> {
    let path = path.as_ref();
    let text = std::fs::read_to_string(path).map_err(|error| {
        xtalforge_core::Error::InvalidInput(format!("could not read {}: {error}", path.display()))
    })?;
    let hint = path
        .extension()
        .and_then(std::ffi::OsStr::to_str)
        .unwrap_or("");
    read_str(&text, hint)
}

/// Reads a structure from text and returns its warnings separately.
///
/// # Errors
///
/// Returns [`xtalforge_core::Error::InvalidInput`] for an unknown format or a
/// reader error.
pub fn read_str(text: &str, hint: &str) -> Result<ReadResult, xtalforge_core::Error> {
    let structure = xtalforge_io::read_structure(text, hint)?;
    Ok(ReadResult {
        warnings: structure.warnings(),
        structure,
    })
}

/// Writes a structure to a path, choosing the format from the extension.
///
/// # Errors
///
/// Returns [`xtalforge_core::Error::InvalidInput`] when the extension is
/// unknown or the writer fails.
pub fn write(
    path: impl AsRef<std::path::Path>,
    structure: &xtalforge_core::Structure,
) -> Result<(), xtalforge_core::Error> {
    let path = path.as_ref();
    let format = path
        .extension()
        .and_then(std::ffi::OsStr::to_str)
        .ok_or_else(|| {
            xtalforge_core::Error::InvalidInput("the output path has no extension".to_owned())
        })?;
    let text = xtalforge_io::write_structure(structure, format)?;
    xtalforge_io::atomic::write_text(path, &text)
}

/// Common imports for typical xtalforge programs.
pub mod prelude {
    #[cfg(feature = "analysis")]
    pub use crate::analysis::calculate_pxrd;
    #[cfg(feature = "build")]
    pub use crate::build::{
        build_carbon, build_polymer, embed_smiles, fullerene, graphene, nanotube,
    };
    pub use crate::core::{
        Element, Error, Lattice, ParallelOptions, Site, Slab, Structure, StructureBuilder,
        StructureWarning, SupercellChoice, cell_volume, element_info, lattice_info,
        lattice_minimum_image, structure_info, version,
    };
    #[cfg(feature = "ff")]
    pub use crate::ff::{Calculator, UffCalculator, minimize};
    pub use crate::io::{read_structure, write_structure};
    pub use crate::symmetry::{SymmetryDataset, SymmetryOptions, detect_symmetry};
    pub use crate::{ReadResult, read, read_str, write};
}
