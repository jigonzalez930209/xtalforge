//! SDF and MDL MOL writer (Task 13.2.5).
//!
//! Writes one V2000 record. Charges use `M  CHG`, bonds use the stored
//! explicit bonds and `sdf_data` metadata becomes data fields before `$$$$`.

use serde_json::Value;
use xtalforge_core::{Error, Structure};

/// Writes a structure as one SDF record with a V2000 molfile.
///
/// # Errors
///
/// Returns an error for an empty structure or a failed P1 expansion.
pub fn write_sdf(structure: &Structure) -> Result<String, Error> {
    let p1 = structure.p1()?;
    if p1.atoms.is_empty() {
        return Err(Error::InvalidInput(
            "an SDF record cannot hold a structure with no atoms".to_owned(),
        ));
    }
    let mut site_atom: Vec<Option<usize>> = vec![None; structure.sites().len()];
    for (index, atom) in p1.atoms.iter().enumerate() {
        if let Some(slot) = site_atom.get_mut(atom.site_idx)
            && slot.is_none()
        {
            *slot = Some(index);
        }
    }
    let bonds: Vec<(usize, usize, f64)> = structure
        .explicit_bonds()
        .iter()
        .filter(|bond| !bond.remove)
        .filter_map(|bond| {
            let i = site_atom.get(bond.i).copied().flatten()?;
            let j = site_atom.get(bond.j).copied().flatten()?;
            (i != j).then_some((i, j, bond.order))
        })
        .collect();
    let mut out = String::new();
    out.push_str(&format!("{}\n", structure.name().replace('\n', " ")));
    out.push_str("  xtalforge 1.0\n");
    out.push_str("structure written by xtalforge\n");
    out.push_str(&format!(
        "{:>3}{:>3}  0  0  0  0  0  0  0  0999 V2000\n",
        p1.atoms.len(),
        bonds.len()
    ));
    for atom in &p1.atoms {
        let [x, y, z] = structure.lattice().to_cart(atom.frac);
        let charge = structure
            .sites()
            .get(atom.site_idx)
            .and_then(|site| site.charge())
            .unwrap_or(0.0);
        out.push_str(&format!(
            "{x:>10.4}{y:>10.4}{z:>10.4} {:<3} 0{:>3}  0  0  0  0  0  0  0  0\n",
            atom.element,
            charge_code(charge)
        ));
    }
    for (i, j, order) in &bonds {
        out.push_str(&format!(
            "{:>3}{:>3}{:>3}  0  0  0  0\n",
            i + 1,
            j + 1,
            bond_code(*order)
        ));
    }
    let charged: Vec<(usize, i32)> = p1
        .atoms
        .iter()
        .enumerate()
        .filter_map(|(index, atom)| {
            let charge = structure
                .sites()
                .get(atom.site_idx)
                .and_then(|site| site.charge())
                .unwrap_or(0.0);
            let rounded = charge.round() as i32;
            (rounded != 0).then_some((index + 1, rounded))
        })
        .collect();
    for chunk in charged.chunks(8) {
        out.push_str(&format!("M  CHG{:>3}", chunk.len()));
        for (atom, charge) in chunk {
            out.push_str(&format!("{atom:>4}{charge:>4}"));
        }
        out.push('\n');
    }
    out.push_str("M  END\n");
    if let Some(Value::Object(data)) = structure.metadata().get("sdf_data") {
        for (key, value) in data {
            if let Some(text) = value.as_str() {
                out.push_str(&format!(">  <{key}>\n{text}\n\n"));
            }
        }
    }
    out.push_str("$$$$\n");
    Ok(out)
}

fn charge_code(charge: f64) -> i32 {
    match charge.round() as i32 {
        3 => 1,
        2 => 2,
        1 => 3,
        -1 => 5,
        -2 => 6,
        -3 => 7,
        _ => 0,
    }
}

fn bond_code(order: f64) -> i32 {
    if (order - 1.5).abs() < 0.2 {
        4
    } else if (order - 3.0).abs() < 0.2 {
        3
    } else if (order - 2.0).abs() < 0.2 {
        2
    } else {
        1
    }
}
