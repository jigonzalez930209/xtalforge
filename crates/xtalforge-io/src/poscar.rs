//! POSCAR/CONTCAR reader and writer (Task 1.16.1).

#![allow(clippy::indexing_slicing)] // fixed-format parsers
use xtalforge_core::{Error, Lattice, Site, Structure, StructureBuilder};

/// Reads a VASP POSCAR or CONTCAR.
///
/// # Errors
///
/// Returns [`Error::InvalidInput`] for a missing species line, a bad scale,
/// a malformed lattice or a coordinate count mismatch.
pub fn read_poscar(text: &str) -> Result<Structure, Error> {
    let lines: Vec<&str> = text.lines().collect();
    if lines.len() < 8 {
        return Err(Error::InvalidInput("POSCAR is too short".to_owned()));
    }
    let comment = lines[0].trim();
    let scale: f64 = lines[1]
        .trim()
        .parse()
        .map_err(|_| Error::InvalidInput("POSCAR scale is not a number".to_owned()))?;
    let mut raw = [[0.0; 3]; 3];
    for (row, line) in lines[2..5].iter().enumerate() {
        let values: Vec<f64> = line
            .split_whitespace()
            .filter_map(|value| value.parse().ok())
            .collect();
        if values.len() < 3 {
            return Err(Error::InvalidInput(
                "POSCAR lattice row has fewer than three values".to_owned(),
            ));
        }
        raw[row] = [values[0], values[1], values[2]];
    }
    let determinant = raw[0][0] * (raw[1][1] * raw[2][2] - raw[1][2] * raw[2][1])
        - raw[0][1] * (raw[1][0] * raw[2][2] - raw[1][2] * raw[2][0])
        + raw[0][2] * (raw[1][0] * raw[2][1] - raw[1][1] * raw[2][0]);
    let factor = if scale < 0.0 {
        (-scale / determinant.abs()).cbrt()
    } else {
        scale
    };
    let matrix = raw.map(|row| row.map(|value| value * factor));
    let species_line = lines[5].trim();
    if species_line.is_empty()
        || species_line
            .chars()
            .all(|character| character.is_ascii_digit())
    {
        return Err(Error::InvalidInput(
            "POSCAR has no species line; VASP 4 files are not supported".to_owned(),
        ));
    }
    let species: Vec<&str> = species_line.split_whitespace().collect();
    let counts: Vec<usize> = lines[6]
        .split_whitespace()
        .filter_map(|value| value.parse().ok())
        .collect();
    if species.len() != counts.len() {
        return Err(Error::InvalidInput(
            "POSCAR species and count lines disagree".to_owned(),
        ));
    }
    let mut index = 7;
    if lines
        .get(index)
        .is_some_and(|line| line.trim().to_ascii_lowercase().starts_with('s'))
    {
        index += 1;
    }
    let mode = lines
        .get(index)
        .map(|line| line.trim().to_ascii_lowercase())
        .unwrap_or_default();
    index += 1;
    let direct = mode.starts_with('d');
    let lattice = Lattice::new(matrix)?;
    let mut builder = StructureBuilder::new(comment, lattice);
    let mut atom = 0;
    for (species, count) in species.iter().zip(counts) {
        for _ in 0..count {
            let line = lines.get(index).ok_or_else(|| {
                Error::InvalidInput("POSCAR ended before all coordinates".to_owned())
            })?;
            let values: Vec<f64> = line
                .split_whitespace()
                .take(3)
                .filter_map(|value| value.parse().ok())
                .collect();
            if values.len() < 3 {
                return Err(Error::InvalidInput(
                    "POSCAR coordinate has fewer than three values".to_owned(),
                ));
            }
            let frac = if direct {
                [values[0], values[1], values[2]]
            } else {
                lattice.to_frac([values[0] * scale, values[1] * scale, values[2] * scale])
            };
            builder = builder.site(Site::new(
                format!("{}{}", species, atom + 1),
                species,
                frac,
            )?);
            atom += 1;
            index += 1;
        }
    }
    builder.build()
}

/// Writes a P1 structure as a POSCAR.
///
/// # Errors
///
/// Returns an error if a site element cannot be represented.
pub fn write_poscar(structure: &Structure) -> Result<String, Error> {
    let p1 = structure.p1()?;
    let mut species: Vec<String> = Vec::new();
    let mut counts: Vec<usize> = Vec::new();
    for atom in &p1.atoms {
        if let Some(index) = species.iter().position(|element| element == &atom.element) {
            if let Some(count) = counts.get_mut(index) {
                *count += 1;
            }
        } else {
            species.push(atom.element.clone());
            counts.push(1);
        }
    }
    let mut out = String::new();
    out.push_str(&format!("{}\n1.0\n", structure.name()));
    for row in structure.lattice().matrix() {
        out.push_str(&format!("{:.10} {:.10} {:.10}\n", row[0], row[1], row[2]));
    }
    out.push_str(&format!("{}\n", species.join(" ")));
    out.push_str(&format!(
        "{}\nDirect\n",
        counts
            .iter()
            .map(usize::to_string)
            .collect::<Vec<_>>()
            .join(" ")
    ));
    for atom in &p1.atoms {
        out.push_str(&format!(
            "{:.10} {:.10} {:.10}\n",
            atom.frac[0], atom.frac[1], atom.frac[2]
        ));
    }
    Ok(out)
}
