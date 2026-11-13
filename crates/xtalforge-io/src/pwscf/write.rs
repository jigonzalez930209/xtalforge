//! Quantum ESPRESSO `pw.x` input writer (Task 13.2.1).

use xtalforge_core::elements::element_by_symbol;
use xtalforge_core::{Error, Structure};

pub(super) fn write(structure: &Structure) -> Result<String, Error> {
    let p1 = structure.p1()?;
    if p1.atoms.is_empty() {
        return Err(Error::InvalidInput(
            "a pw.x input deck cannot hold a structure with no atoms".to_owned(),
        ));
    }
    let mut species: Vec<String> = Vec::new();
    for atom in &p1.atoms {
        let symbol = element_by_symbol(&atom.element)?.symbol.clone();
        if !species.iter().any(|known| known == &symbol) {
            species.push(symbol);
        }
    }
    let prefix = sanitize_name(structure.name());
    let mut out = String::from("&CONTROL\n");
    out.push_str("  calculation = 'scf'\n");
    out.push_str(&format!("  prefix = '{prefix}'\n"));
    out.push_str("  pseudo_dir = './'\n");
    out.push_str("  outdir = './tmp'\n");
    out.push_str("/\n&SYSTEM\n  ibrav = 0\n");
    out.push_str(&format!("  nat = {}\n", p1.atoms.len()));
    out.push_str(&format!("  ntyp = {}\n", species.len()));
    out.push_str("  ecutwfc = 50.0\n");
    out.push_str("/\n&ELECTRONS\n  conv_thr = 1.0d-8\n/\n");
    out.push_str("ATOMIC_SPECIES\n");
    for symbol in &species {
        let mass = element_by_symbol(symbol)?.mass;
        out.push_str(&format!("  {symbol} {mass:.4} {symbol}.UPF\n"));
    }
    out.push_str("CELL_PARAMETERS angstrom\n");
    for row in structure.lattice().matrix() {
        out.push_str(&format!("  {:.10} {:.10} {:.10}\n", row[0], row[1], row[2]));
    }
    out.push_str("ATOMIC_POSITIONS crystal\n");
    for atom in &p1.atoms {
        out.push_str(&format!(
            "  {} {:.10} {:.10} {:.10}\n",
            atom.element, atom.frac[0], atom.frac[1], atom.frac[2]
        ));
    }
    out.push_str("K_POINTS gamma\n");
    Ok(out)
}

fn sanitize_name(name: &str) -> String {
    let cleaned: String = name
        .chars()
        .map(|character| {
            if character.is_ascii_alphanumeric() || character == '_' {
                character
            } else {
                '_'
            }
        })
        .collect();
    if cleaned.is_empty() {
        "xtalforge".to_owned()
    } else {
        cleaned
    }
}
