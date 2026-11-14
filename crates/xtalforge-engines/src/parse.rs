//! Engine output parsers (Task 18.2).
//!
//! Every parser returns a [`CalculationResult`](crate::result::CalculationResult)
//! with the same units: Angstrom, eV, eV/Angstrom and GPa. VASP keeps one
//! entry point per output file; Quantum ESPRESSO separates the text output,
//! the schema XML, the bands file and `projwfc`; CP2K, LAMMPS and GULP have
//! one entry point each, with the LAMMPS log and dump available separately.

mod common;
mod cp2k;
mod gulp;
mod lammps;
mod qe_bands;
mod qe_projwfc;
mod qe_text;
mod qe_xml;
mod vasp_doscar;
mod vasp_eigenval;
mod vasp_oszicar;
mod vasp_outcar;
mod vasp_procar;
mod vasp_vasprun;
mod xml;

pub use cp2k::parse as parse_cp2k;
pub use gulp::parse as parse_gulp;
pub use lammps::parse_dump_with_types as parse_lammps_dump_with_types;
pub use lammps::parse_log as parse_lammps_log;
pub use lammps::{parse as parse_lammps, parse_dump as parse_lammps_dump};
pub use qe_bands::parse as parse_qe_bands;
pub use qe_projwfc::parse as parse_qe_projwfc;
pub use qe_text::parse as parse_qe;
pub use qe_xml::parse as parse_qe_xml;
pub use vasp_doscar::parse as parse_vasp_doscar;
pub use vasp_eigenval::parse as parse_vasp_eigenval;
pub use vasp_oszicar::parse as parse_vasp_oszicar;
pub use vasp_outcar::parse as parse_vasp_outcar;
pub use vasp_procar::parse as parse_vasp_procar;
pub use vasp_vasprun::parse as parse_vasp_vasprun;
