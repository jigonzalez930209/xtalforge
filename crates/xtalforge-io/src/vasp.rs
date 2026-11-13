//! VASP output readers: `OUTCAR` and `vasprun.xml` (Task 13.2.4).
//!
//! Both readers produce a [`crate::trajectory::Trajectory`] with one frame per
//! ionic step, energies in eV, forces in eV/Angstrom and the stress tensor in
//! GPa. The structure entry points return the last frame so the format
//! registry can expose the files as single structures. VASP outputs are
//! read-only.

mod outcar;
mod vasprun;
mod xml;

pub use outcar::{read_outcar, read_outcar_structure};
pub use vasprun::{read_vasprun, read_vasprun_structure};

/// Whether the text looks like a `vasprun.xml`.
#[must_use]
pub fn looks_like_vasprun(text: &str) -> bool {
    text.contains("<modeling") || text.contains("<vasprun")
}

/// Whether the text looks like an `OUTCAR`.
#[must_use]
pub fn looks_like_outcar(text: &str) -> bool {
    text.contains("vasp.")
        || (text.contains("POSITION") && text.contains("TOTAL-FORCE"))
        || (text.contains("TITEL") && text.contains("ions per type"))
}
