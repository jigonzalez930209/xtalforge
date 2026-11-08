//! CSSR reader and writer for Zeo++ (Task 3.3.2).

#![allow(clippy::indexing_slicing)] // fixed-format parser, lengths checked above
use xtalforge_core::{Error, Lattice, Site, Structure, StructureBuilder};

/// Reads CSSR text.
///
/// # Errors
///
/// Returns [`Error::InvalidInput`] for a short header, bad cell parameters or
/// malformed atom lines.
pub fn read_cssr(text: &str) -> Result<Structure, Error> {
    let lines: Vec<&str> = text.lines().collect();
    if lines.len() < 4 {
        return Err(Error::InvalidInput(
            "a CSSR file has four header lines".to_owned(),
        ));
    }
    let cell: Vec<f64> = lines
        .first()
        .unwrap_or(&"")
        .split_whitespace()
        .filter_map(|value| value.parse().ok())
        .collect();
    let angles: Vec<f64> = lines
        .get(1)
        .unwrap_or(&"")
        .split_whitespace()
        .filter_map(|value| value.parse().ok())
        .collect();
    if cell.len() < 3 || angles.len() < 3 {
        return Err(Error::InvalidInput(
            "the first two CSSR lines are a b c and alpha beta gamma".to_owned(),
        ));
    }
    let lattice =
        Lattice::from_parameters(cell[0], cell[1], cell[2], angles[0], angles[1], angles[2])?;
    let counted: Vec<&str> = lines.get(2).unwrap_or(&"").split_whitespace().collect();
    let cartesian = counted.get(1).is_some_and(|flag| *flag == "1");
    let count: usize = counted
        .first()
        .and_then(|value| value.parse().ok())
        .or_else(|| {
            lines
                .iter()
                .skip(4)
                .filter(|line| line.split_whitespace().count() >= 5)
                .count()
                .checked_sub(0)
        })
        .unwrap_or(0);
    let mut builder = StructureBuilder::new("cssr", lattice);
    for index in 0..count {
        let line = lines
            .get(4 + index)
            .ok_or_else(|| Error::InvalidInput("CSSR ended before all atoms".to_owned()))?;
        let parts: Vec<&str> = line.split_whitespace().collect();
        if parts.len() < 5 {
            return Err(Error::InvalidInput(
                "a CSSR atom line needs index, element and three coordinates".to_owned(),
            ));
        }
        let element = parts[1];
        let coordinates = [
            parts[2]
                .parse::<f64>()
                .map_err(|_| Error::InvalidInput("bad CSSR x".to_owned()))?,
            parts[3]
                .parse::<f64>()
                .map_err(|_| Error::InvalidInput("bad CSSR y".to_owned()))?,
            parts[4]
                .parse::<f64>()
                .map_err(|_| Error::InvalidInput("bad CSSR z".to_owned()))?,
        ];
        let frac = if cartesian {
            lattice.to_frac(coordinates)
        } else {
            coordinates
        };
        builder = builder.site(Site::new(format!("{element}{}", index + 1), element, frac)?);
    }
    builder.build()
}

/// Writes a P1 CSSR for Zeo++.
///
/// # Errors
///
/// Returns the errors of the P1 expansion or an empty structure.
pub fn write_cssr(structure: &Structure) -> Result<String, Error> {
    let p1 = structure.p1()?;
    if p1.atoms.is_empty() {
        return Err(Error::InvalidInput(
            "a CSSR file cannot hold a structure with no atoms".to_owned(),
        ));
    }
    let [a, b, c, alpha, beta, gamma] = structure.lattice().parameters();
    let mut out = String::new();
    out.push_str(&format!("\t\t\t\t{a:.4}  {b:.4}  {c:.4}\n"));
    out.push_str(&format!(
        "\t\t{alpha:.3}  {beta:.3}  {gamma:.3}  SPGR =  1 P 1\t\t OPT = 1\n"
    ));
    out.push_str(&format!("{}   0\n", p1.atoms.len()));
    out.push_str("0 structure\n");
    for (index, atom) in p1.atoms.iter().enumerate() {
        let [x, y, z] = atom.frac;
        out.push_str(&format!(
            " {} {:<3} {x:9.5} {y:9.5} {z:9.5} 0 0 0 0 0 0 0 0 {:.4}\n",
            index + 1,
            atom.element,
            0.0
        ));
    }
    Ok(out)
}
