//! XYZ and extended XYZ reader and writer (Task 1.16.2).

#![allow(clippy::indexing_slicing)] // fixed-format parsers
use xtalforge_core::{Error, Lattice, Site, Structure, StructureBuilder};

/// Reads plain or extended XYZ.
///
/// Extended XYZ carries `Lattice="..."`; a plain file is placed in a padded
/// P1 box.
///
/// # Errors
///
/// Returns [`Error::InvalidInput`] for a bad atom count, a malformed line or
/// an invalid lattice.
pub fn read_xyz(text: &str) -> Result<Structure, Error> {
    let mut lines = text.lines();
    let count: usize = lines
        .next()
        .and_then(|line| line.trim().parse().ok())
        .ok_or_else(|| Error::InvalidInput("XYZ atom count is missing".to_owned()))?;
    let comment = lines.next().unwrap_or("").trim();
    let lattice = parse_lattice(comment)?;
    let mut points = Vec::new();
    for _ in 0..count {
        let line = lines
            .next()
            .ok_or_else(|| Error::InvalidInput("XYZ ended before all atoms".to_owned()))?;
        let values: Vec<&str> = line.split_whitespace().collect();
        if values.len() < 4 {
            return Err(Error::InvalidInput(
                "XYZ atom line needs element and three coordinates".to_owned(),
            ));
        }
        let element = values[0];
        let coordinates = [
            values[1]
                .parse()
                .map_err(|_| Error::InvalidInput("bad XYZ x".to_owned()))?,
            values[2]
                .parse()
                .map_err(|_| Error::InvalidInput("bad XYZ y".to_owned()))?,
            values[3]
                .parse()
                .map_err(|_| Error::InvalidInput("bad XYZ z".to_owned()))?,
        ];
        points.push((element.to_owned(), coordinates));
    }
    match lattice {
        Some(lattice) => {
            let mut builder = StructureBuilder::new(comment, lattice);
            for (index, (element, point)) in points.iter().enumerate() {
                builder = builder.site(Site::new(
                    format!("{element}{}", index + 1),
                    element,
                    lattice.to_frac(*point),
                )?);
            }
            builder.build()
        }
        None => Structure::molecule("xyz", &points, 5.0),
    }
}

/// Writes an extended XYZ with `Lattice` and `Properties`.
///
/// # Errors
///
/// Returns an error if the structure cannot be expanded.
pub fn write_xyz(structure: &Structure) -> Result<String, Error> {
    let p1 = structure.p1()?;
    let matrix = structure.lattice().matrix();
    let lattice = matrix
        .iter()
        .flat_map(|row| row.iter())
        .map(|value| format!("{value:.10}"))
        .collect::<Vec<_>>()
        .join(" ");
    let mut out = format!(
        "{}\nLattice=\"{}\" Properties=species:S:1:pos:R:3\n",
        p1.atoms.len(),
        lattice
    );
    for atom in &p1.atoms {
        let cart = structure.lattice().to_cart(atom.frac);
        out.push_str(&format!(
            "{} {:.10} {:.10} {:.10}\n",
            atom.element, cart[0], cart[1], cart[2]
        ));
    }
    Ok(out)
}

fn parse_lattice(comment: &str) -> Result<Option<Lattice>, Error> {
    let Some(start) = comment.find("Lattice=\"") else {
        return Ok(None);
    };
    let rest = &comment[start + 9..];
    let Some(end) = rest.find('"') else {
        return Err(Error::InvalidInput(
            "XYZ Lattice attribute is not closed".to_owned(),
        ));
    };
    let values: Vec<f64> = rest[..end]
        .split_whitespace()
        .filter_map(|value| value.parse().ok())
        .collect();
    if values.len() != 9 {
        return Err(Error::InvalidInput(
            "XYZ Lattice needs nine numbers".to_owned(),
        ));
    }
    Ok(Some(Lattice::new([
        [values[0], values[1], values[2]],
        [values[3], values[4], values[5]],
        [values[6], values[7], values[8]],
    ])?))
}
