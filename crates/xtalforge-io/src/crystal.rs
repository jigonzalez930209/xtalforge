//! CRYSTAL input reader and writer (Task 13.2.2).
//!
//! Reads the three-dimensional `CRYSTAL` geometry block: the optional
//! `IFLAG IFHR IFSO` flags, the space-group number, six cell parameters, the
//! number of atoms and one `Z x y z` line per atom in fractional coordinates.
//! Slabs, polymers and molecules (`IFLAG` 1-3) are outside this subset and
//! return a stable error. The writer emits a P1 block with `END`.

#![allow(clippy::indexing_slicing)] // fixed-format parser, lengths checked

use xtalforge_core::elements::{atomic_number, element_by_z};
use xtalforge_core::{Error, Lattice, Site, Structure, StructureBuilder};

use crate::molecule::parse_number;

/// Reads a CRYSTAL geometry block.
///
/// # Errors
///
/// Returns [`Error::InvalidInput`] for a missing `CRYSTAL` keyword, a
/// non-3D `IFLAG`, an unknown space group or a malformed atom line.
pub fn read_crystal(text: &str) -> Result<Structure, Error> {
    let lines: Vec<&str> = text.lines().collect();
    let title = lines.first().map_or("crystal", |line| line.trim());
    let keyword = lines
        .iter()
        .position(|line| line.trim().eq_ignore_ascii_case("CRYSTAL"))
        .ok_or_else(|| Error::InvalidInput("CRYSTAL input has no CRYSTAL keyword".to_owned()))?;
    let mut index = keyword + 1;
    skip_empty(&lines, &mut index);
    if let Some(flags) = lines.get(index).map(|line| line.trim())
        && flags.split_whitespace().count() == 3
        && flags
            .split_whitespace()
            .all(|token| token.parse::<i32>().is_ok())
    {
        let iflag = flags
            .split_whitespace()
            .next()
            .and_then(|token| token.parse::<i32>().ok())
            .unwrap_or(0);
        if iflag != 0 {
            return Err(Error::InvalidInput(format!(
                "CRYSTAL IFLAG={iflag} is not supported; only 3D crystals (IFLAG=0) are read"
            )));
        }
        index += 1;
    }
    skip_empty(&lines, &mut index);
    let group_line = lines.get(index).map(|line| line.trim()).ok_or_else(|| {
        Error::InvalidInput("CRYSTAL input ends before the space group".to_owned())
    })?;
    index += 1;
    let (group, hall) = match group_line.parse::<u32>() {
        Ok(0) => (None, "P 1".to_owned()),
        Ok(number) => {
            let group = xtalforge_symmetry::SpaceGroup::from_number(number)?;
            let hall = group.hall_symbol.clone();
            (Some(group.operations), hall)
        }
        Err(_) => {
            let group = xtalforge_symmetry::SpaceGroup::parse(group_line)?;
            (Some(group.operations), group.hall_symbol.clone())
        }
    };
    skip_empty(&lines, &mut index);
    let cell_line = lines
        .get(index)
        .map(|line| line.trim())
        .ok_or_else(|| Error::InvalidInput("CRYSTAL input ends before the cell".to_owned()))?;
    index += 1;
    let cell: Vec<f64> = cell_line
        .split_whitespace()
        .filter_map(|token| token.parse().ok())
        .collect();
    if cell.len() < 6 {
        return Err(Error::InvalidInput(
            "CRYSTAL cell needs a b c alpha beta gamma".to_owned(),
        ));
    }
    let lattice = Lattice::from_parameters(cell[0], cell[1], cell[2], cell[3], cell[4], cell[5])?;
    skip_empty(&lines, &mut index);
    let count: usize = lines
        .get(index)
        .map(|line| line.trim())
        .and_then(|line| line.parse().ok())
        .ok_or_else(|| Error::InvalidInput("CRYSTAL input has no atom count".to_owned()))?;
    index += 1;
    let mut builder = StructureBuilder::new(title, lattice).hall(hall);
    if let Some(operations) = group {
        for operation in operations {
            builder = builder.operation(operation);
        }
    }
    for atom_index in 0..count {
        let line = lines
            .get(index)
            .ok_or_else(|| Error::InvalidInput("CRYSTAL ended before all atoms".to_owned()))?;
        index += 1;
        let fields: Vec<&str> = line.split_whitespace().collect();
        if fields.len() < 4 {
            return Err(Error::InvalidInput(
                "a CRYSTAL atom line needs Z x y z".to_owned(),
            ));
        }
        let z: u8 = fields[0]
            .parse()
            .map_err(|_| Error::InvalidInput(format!("bad CRYSTAL atomic number: {fields:?}")))?;
        let element = element_by_z(z)?.symbol.clone();
        let frac = [
            parse_number(fields[1], "CRYSTAL x")?,
            parse_number(fields[2], "CRYSTAL y")?,
            parse_number(fields[3], "CRYSTAL z")?,
        ];
        builder = builder.site(Site::new(
            format!("{element}{}", atom_index + 1),
            &element,
            frac,
        )?);
    }
    builder.build()
}

/// Writes a P1 structure as a CRYSTAL geometry block.
///
/// # Errors
///
/// Returns an error for an empty structure, a failed P1 expansion or an
/// element without an atomic number.
pub fn write_crystal(structure: &Structure) -> Result<String, Error> {
    let p1 = structure.p1()?;
    if p1.atoms.is_empty() {
        return Err(Error::InvalidInput(
            "a CRYSTAL input cannot hold a structure with no atoms".to_owned(),
        ));
    }
    let (number, sites): (u32, Vec<(String, [f64; 3])>) =
        match xtalforge_symmetry::SpaceGroup::parse(structure.hall()) {
            Ok(group) if !structure.sites().is_empty() => (
                group.number,
                structure
                    .sites()
                    .iter()
                    .map(|site| (site.element().to_owned(), site.frac()))
                    .collect(),
            ),
            _ => (
                1,
                p1.atoms
                    .iter()
                    .map(|atom| (atom.element.clone(), atom.frac))
                    .collect(),
            ),
        };
    let [a, b, c, alpha, beta, gamma] = structure.lattice().parameters();
    let mut out = format!("{}\nCRYSTAL\n0 0 0\n{number}\n", structure.name());
    out.push_str(&format!(
        "{a:.10} {b:.10} {c:.10} {alpha:.6} {beta:.6} {gamma:.6}\n"
    ));
    out.push_str(&format!("{}\n", sites.len()));
    for (element, frac) in &sites {
        let z = atomic_number(element)?;
        out.push_str(&format!(
            "{} {:.10} {:.10} {:.10}\n",
            z, frac[0], frac[1], frac[2]
        ));
    }
    out.push_str("END\n");
    Ok(out)
}

fn skip_empty(lines: &[&str], index: &mut usize) {
    while lines.get(*index).is_some_and(|line| line.trim().is_empty()) {
        *index += 1;
    }
}
