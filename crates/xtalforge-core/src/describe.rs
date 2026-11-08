//! One-line atom descriptions for tooltips, the CLI and the agent layer.

use std::fmt::Write;

use crate::{Error, Structure};

/// Describes one P1 atom: element, label, site, coordination and charge.
///
/// # Errors
///
/// Returns [`Error::InvalidInput`] for an out-of-range index or a failed
/// expansion.
pub fn describe_atom(structure: &Structure, index: usize) -> Result<String, Error> {
    let p1 = structure.p1()?;
    let atom = p1
        .atoms
        .get(index)
        .ok_or_else(|| Error::InvalidInput(format!("atom index {index} is out of range")))?;
    let graph = structure.bond_graph()?;
    let coordination = graph.degree(index);
    let mut description = String::new();
    let _ = write!(
        description,
        "{} {} ({}) site={} coordination={}",
        atom.element, atom.label, atom.frac[0], atom.site_idx, coordination
    );
    let _ = write!(
        description,
        " frac=({:.4}, {:.4}, {:.4})",
        atom.frac[0], atom.frac[1], atom.frac[2]
    );
    if let Some(site) = structure.sites().get(atom.site_idx) {
        if let Some(charge) = site.charge() {
            let _ = write!(description, " charge={charge:+.2}");
        }
        if site.occupancy() < 1.0 {
            let _ = write!(description, " occupancy={:.2}", site.occupancy());
        }
    }
    Ok(description)
}
