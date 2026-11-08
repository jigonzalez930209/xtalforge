//! PDB writer for finite atom sets (Task 3.3.6).
//!
//! Ported from Crystal-Builder `xtal/io/pdb.py` (MIT). The writer expands the
//! structure to P1, writes every atom as a `HETATM` record with its element in
//! columns 77-78, and writes one `CONECT` record per explicit bond. Residues
//! are not modelled: every atom sits in a single `UNL` residue.

use std::collections::BTreeMap;

use xtalforge_core::{Error, Structure};

/// The widest serial a PDB five-column field holds.
const MAX_ATOMS: usize = 99_999;

/// Partners per `CONECT` line; longer lists continue on another line.
const PARTNERS_PER_LINE: usize = 4;

/// Writes a structure as a PDB file.
///
/// # Errors
///
/// Returns [`Error::InvalidInput`] for an empty structure, more than 99999
/// atoms or a failed P1 expansion.
pub fn write_pdb(structure: &Structure) -> Result<String, Error> {
    let p1 = structure.p1()?;
    if p1.atoms.is_empty() {
        return Err(Error::InvalidInput(
            "a PDB file cannot hold a structure with no atoms".to_owned(),
        ));
    }
    if p1.atoms.len() > MAX_ATOMS {
        return Err(Error::InvalidInput(format!(
            "{} atoms is more than a PDB can number ({MAX_ATOMS})",
            p1.atoms.len()
        )));
    }

    let mut out = String::new();
    for (index, atom) in p1.atoms.iter().enumerate() {
        let [x, y, z] = structure.lattice().to_cart(atom.frac);
        out.push_str(&hetatm_record(
            index + 1,
            &atom.element,
            [x, y, z],
            atom.occupancy,
        ));
        out.push('\n');
    }

    let mut site_atom: Vec<Option<usize>> = vec![None; structure.sites().len()];
    for (index, atom) in p1.atoms.iter().enumerate() {
        if let Some(slot) = site_atom.get_mut(atom.site_idx)
            && slot.is_none()
        {
            *slot = Some(index);
        }
    }
    let mut partners: BTreeMap<usize, Vec<usize>> = BTreeMap::new();
    for bond in structure.explicit_bonds() {
        if bond.remove {
            continue;
        }
        let Some(i) = site_atom.get(bond.i).copied().flatten() else {
            continue;
        };
        let Some(j) = site_atom.get(bond.j).copied().flatten() else {
            continue;
        };
        if i == j {
            continue;
        }
        partners.entry(i + 1).or_default().push(j + 1);
        partners.entry(j + 1).or_default().push(i + 1);
    }
    for (serial, mut bonded) in partners {
        bonded.sort_unstable();
        bonded.dedup();
        for chunk in bonded.chunks(PARTNERS_PER_LINE) {
            out.push_str(&format!("CONECT{serial:5}"));
            for partner in chunk {
                out.push_str(&format!("{partner:5}"));
            }
            out.push('\n');
        }
    }
    out.push_str("END\n");
    Ok(out)
}

fn hetatm_record(serial: usize, element: &str, cart: [f64; 3], occupancy: f64) -> String {
    let element = element.to_ascii_uppercase();
    let name = format!("{element:>2}  ");
    let [x, y, z] = cart;
    format!(
        "HETATM{serial:5} {name:4} UNL A   1    \
         {x:8.3}{y:8.3}{z:8.3}{occupancy:6.2}{:6.2}          {element:>2}",
        0.0_f64
    )
}
