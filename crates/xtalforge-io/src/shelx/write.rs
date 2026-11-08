//! SHELX instruction-file writer (Task 3.3.5).
//!
//! Ported from Crystal-Builder `xtal/io/shelx.py` (MIT). The writer emits P1
//! with `LATT 1` and no free variables.

use xtalforge_core::elements::element_by_symbol;
use xtalforge_core::{Error, Site, Structure};

/// Writes a P1 SHELX instruction file.
///
/// # Errors
///
/// Returns [`Error::InvalidInput`] for an empty structure or a failed P1
/// expansion.
pub fn write_shelx(structure: &Structure) -> Result<String, Error> {
    let p1 = structure.p1()?;
    if p1.atoms.is_empty() {
        return Err(Error::InvalidInput(
            "a SHELX file cannot hold a structure with no atoms".to_owned(),
        ));
    }
    let [a, b, c, alpha, beta, gamma] = structure.lattice().parameters();
    let mut species: Vec<String> = Vec::new();
    for atom in &p1.atoms {
        let symbol = element_by_symbol(&atom.element)?.symbol.clone();
        if !species.iter().any(|known| known == &symbol) {
            species.push(symbol);
        }
    }
    let mut out = format!("TITL {}\n", structure.name());
    out.push_str(&format!(
        "CELL 0.71073 {a:.6} {b:.6} {c:.6} {alpha:.6} {beta:.6} {gamma:.6}\n"
    ));
    out.push_str("LATT 1\n");
    out.push_str(&format!("SFAC {}\n", species.join(" ")));
    for (index, atom) in p1.atoms.iter().enumerate() {
        let symbol = element_by_symbol(&atom.element)?.symbol.clone();
        let sfac = species
            .iter()
            .position(|known| known == &symbol)
            .map_or(1, |position| position + 1);
        let site = structure.sites().get(atom.site_idx);
        let u_iso = site.and_then(Site::u_iso).unwrap_or(0.02);
        let group = site.and_then(Site::disorder_group).unwrap_or(0);
        let [x, y, z] = atom.frac;
        let sof = 10.0 + atom.occupancy;
        out.push_str(&format!(
            "{}{} {sfac} {x:.6} {y:.6} {z:.6} {sof:.5} {u_iso:.5} {group}\n",
            symbol,
            index + 1
        ));
    }
    out.push_str("END\n");
    Ok(out)
}
