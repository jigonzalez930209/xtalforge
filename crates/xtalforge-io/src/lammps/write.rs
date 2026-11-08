//! LAMMPS data writer (Task 3.3.4).
//!
//! Ported from Crystal-Builder `xtal/io/lammps.py` (MIT). The writer expands
//! the structure to P1 and takes the bonds from the stored explicit bonds.

use xtalforge_core::elements::element_by_symbol;
use xtalforge_core::{Error, Site, Structure};

/// Writes a P1 LAMMPS data file with `atom_style full`.
///
/// # Errors
///
/// Returns [`Error::InvalidInput`] for an empty structure or a failed P1
/// expansion.
pub fn write_lammps(structure: &Structure) -> Result<String, Error> {
    let p1 = structure.p1()?;
    if p1.atoms.is_empty() {
        return Err(Error::InvalidInput(
            "a LAMMPS data file cannot hold a structure with no atoms".to_owned(),
        ));
    }
    let lattice = structure.lattice().standard_orientation()?;
    let [row_a, row_b, row_c] = lattice.matrix();
    let [xhi, _, _] = row_a;
    let [xy, yhi, _] = row_b;
    let [xz, yz, zhi] = row_c;
    if xhi <= 0.0 || yhi <= 0.0 || zhi <= 0.0 {
        return Err(Error::InvalidInput(
            "the LAMMPS box needs positive x, y and z lengths".to_owned(),
        ));
    }

    let mut species: Vec<String> = Vec::new();
    for atom in &p1.atoms {
        let symbol = element_by_symbol(&atom.element)?.symbol.clone();
        if !species.iter().any(|known| known == &symbol) {
            species.push(symbol);
        }
    }
    let type_of = |symbol: &str| -> usize {
        species
            .iter()
            .position(|known| known == symbol)
            .map_or(1, |index| index + 1)
    };

    let mut site_atom: Vec<Option<usize>> = vec![None; structure.sites().len()];
    for (index, atom) in p1.atoms.iter().enumerate() {
        if let Some(slot) = site_atom.get_mut(atom.site_idx)
            && slot.is_none()
        {
            *slot = Some(index);
        }
    }
    let mut bonds: Vec<(usize, usize)> = Vec::new();
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
        if i != j && !bonds.contains(&(i, j)) && !bonds.contains(&(j, i)) {
            bonds.push((i, j));
        }
    }
    let mut pair_types: Vec<(String, String)> = Vec::new();
    for (i, j) in &bonds {
        let first = p1.atoms.get(*i).map_or("X", |atom| atom.element.as_str());
        let second = p1.atoms.get(*j).map_or("X", |atom| atom.element.as_str());
        let pair = if first <= second {
            (first.to_owned(), second.to_owned())
        } else {
            (second.to_owned(), first.to_owned())
        };
        if !pair_types.contains(&pair) {
            pair_types.push(pair);
        }
    }

    let mut out = String::from("LAMMPS data file written by xtalforge\n\n");
    out.push_str(&format!(
        "{} atoms\n{} bonds\n\n",
        p1.atoms.len(),
        bonds.len()
    ));
    out.push_str(&format!(
        "{} atom types\n{} bond types\n\n",
        species.len(),
        pair_types.len()
    ));
    out.push_str(&format!("0.0 {xhi:.10} xlo xhi\n"));
    out.push_str(&format!("0.0 {yhi:.10} ylo yhi\n"));
    out.push_str(&format!("0.0 {zhi:.10} zlo zhi\n"));
    if xy.abs() > 1e-10 || xz.abs() > 1e-10 || yz.abs() > 1e-10 {
        out.push_str(&format!("{xy:.10} {xz:.10} {yz:.10} xy xz yz\n"));
    }
    out.push_str("\nMasses\n\n");
    for (index, symbol) in species.iter().enumerate() {
        let mass = element_by_symbol(symbol)?.mass;
        out.push_str(&format!("{} {mass:.6}  # {symbol}\n", index + 1));
    }
    out.push_str("\nAtoms  # full\n\n");
    for (index, atom) in p1.atoms.iter().enumerate() {
        let [x, y, z] = lattice.to_cart(atom.frac);
        let charge = structure
            .sites()
            .get(atom.site_idx)
            .and_then(Site::charge)
            .unwrap_or(0.0);
        out.push_str(&format!(
            "{} 1 {} {charge:.6} {x:.10} {y:.10} {z:.10}\n",
            index + 1,
            type_of(&atom.element)
        ));
    }
    if !bonds.is_empty() {
        out.push_str("\nBonds\n\n");
        for (index, (i, j)) in bonds.iter().enumerate() {
            let first = p1.atoms.get(*i).map_or("X", |atom| atom.element.as_str());
            let second = p1.atoms.get(*j).map_or("X", |atom| atom.element.as_str());
            let pair = if first <= second {
                (first.to_owned(), second.to_owned())
            } else {
                (second.to_owned(), first.to_owned())
            };
            let bond_type = pair_types
                .iter()
                .position(|known| known == &pair)
                .map_or(1, |position| position + 1);
            out.push_str(&format!("{} {bond_type} {} {}\n", index + 1, i + 1, j + 1));
        }
    }
    Ok(out)
}
