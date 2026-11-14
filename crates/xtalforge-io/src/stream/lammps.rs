//! Streaming LAMMPS text-dump reader (Task 15.2.2).

use std::collections::BTreeMap;
use std::io::BufRead;

use xtalforge_core::{Error, Lattice};

use super::{StreamFrame, read_line};

/// Streaming LAMMPS text-dump reader.
///
/// Frames use the same header layout as [`crate::read_lammps_dump`] and accept
/// an `element` column or a numeric `type` column mapped through `types` (type
/// number minus one).
pub struct LammpsStream<R: BufRead> {
    reader: R,
    line: String,
    types: Vec<String>,
}

impl<R: BufRead> LammpsStream<R> {
    /// Wraps a buffered reader without a type mapping.
    #[must_use]
    pub fn new(reader: R) -> Self {
        Self {
            reader,
            line: String::new(),
            types: Vec::new(),
        }
    }

    /// Wraps a buffered reader with a LAMMPS type-to-element mapping.
    #[must_use]
    pub fn with_types(reader: R, types: &[&str]) -> Self {
        Self {
            reader,
            line: String::new(),
            types: types.iter().map(|value| (*value).to_owned()).collect(),
        }
    }

    /// Reads the next frame, or `None` at end of file.
    ///
    /// # Errors
    ///
    /// Returns [`Error::InvalidInput`] for a truncated frame, a bad box, an
    /// unknown type or a missing element/type column.
    pub fn next_frame(&mut self) -> Result<Option<StreamFrame>, Error> {
        if !self.seek_timestep()? {
            return Ok(None);
        }
        let step = self
            .required("LAMMPS timestep")?
            .trim()
            .parse::<u64>()
            .map_err(|_| Error::InvalidInput("LAMMPS timestep is not an integer".to_owned()))?;
        self.expect("ITEM: NUMBER OF ATOMS")?;
        let count = self
            .required("LAMMPS atom count")?
            .trim()
            .parse::<usize>()
            .map_err(|_| Error::InvalidInput("LAMMPS atom count is not an integer".to_owned()))?;
        let box_header = self.required("box header")?;
        if !box_header.starts_with("ITEM: BOX BOUNDS") {
            return Err(Error::InvalidInput(
                "LAMMPS dump is missing ITEM: BOX BOUNDS".to_owned(),
            ));
        }
        let mut bounds = [[0.0; 3]; 3];
        for bound in &mut bounds {
            let line = self.required("box bounds")?;
            let values: Vec<f64> = line
                .split_whitespace()
                .map(|value| value.parse::<f64>())
                .collect::<Result<_, _>>()
                .map_err(|_| Error::InvalidInput("LAMMPS box has a bad number".to_owned()))?;
            if values.len() < 2 {
                return Err(Error::InvalidInput(
                    "LAMMPS box row needs two bounds".to_owned(),
                ));
            }
            let first = values.first().copied().unwrap_or(0.0);
            let second = values.get(1).copied().unwrap_or(0.0);
            *bound = [first, second, values.get(2).copied().unwrap_or(0.0)];
        }
        let cell = lammps_lattice(bounds)?;
        let atom_header = self.required("atoms header")?;
        if !atom_header.starts_with("ITEM: ATOMS") {
            return Err(Error::InvalidInput(
                "LAMMPS dump is missing ITEM: ATOMS".to_owned(),
            ));
        }
        let columns: Vec<&str> = atom_header.split_whitespace().skip(2).collect();
        let element_column = columns.iter().position(|column| *column == "element");
        let type_column = columns.iter().position(|column| *column == "type");
        let (mode, x, y, z) = coordinate_columns(&columns)?;
        let mut elements = Vec::with_capacity(count);
        let mut frac = Vec::with_capacity(count);
        for _ in 0..count {
            let line = self.required("atom row")?;
            let values: Vec<&str> = line.split_whitespace().collect();
            let element = if let Some(column) = element_column {
                values
                    .get(column)
                    .ok_or_else(|| truncated("element column"))?
                    .to_string()
            } else if let Some(column) = type_column {
                let kind = values
                    .get(column)
                    .ok_or_else(|| truncated("type column"))?
                    .parse::<usize>()
                    .map_err(|_| Error::InvalidInput("LAMMPS type is not an integer".to_owned()))?;
                self.types
                    .get(kind.saturating_sub(1))
                    .cloned()
                    .ok_or_else(|| {
                        Error::InvalidInput(format!("LAMMPS type {kind} has no element mapping"))
                    })?
            } else {
                return Err(Error::InvalidInput(
                    "LAMMPS dump needs an element column or a type mapping".to_owned(),
                ));
            };
            let point = [
                field_number(&values, x, "coordinate column")?,
                field_number(&values, y, "coordinate column")?,
                field_number(&values, z, "coordinate column")?,
            ];
            elements.push(element);
            frac.push(match mode {
                CoordinateMode::Fractional => point,
                CoordinateMode::Cartesian => cell.to_frac(point),
            });
        }
        Ok(Some(StreamFrame {
            elements,
            frac,
            lattice: Some(cell.matrix()),
            energy: None,
            step: Some(step),
            info: BTreeMap::new(),
            forces: None,
            columns: BTreeMap::new(),
        }))
    }

    fn seek_timestep(&mut self) -> Result<bool, Error> {
        loop {
            let Some(line) = read_line(&mut self.reader, &mut self.line)? else {
                return Ok(false);
            };
            if line.trim() == "ITEM: TIMESTEP" {
                return Ok(true);
            }
        }
    }

    fn required(&mut self, part: &str) -> Result<String, Error> {
        read_line(&mut self.reader, &mut self.line)?.ok_or_else(|| truncated(part))
    }

    fn expect(&mut self, expected: &str) -> Result<(), Error> {
        let line = self.required(expected)?;
        if line.trim() == expected {
            Ok(())
        } else {
            Err(Error::InvalidInput(format!(
                "LAMMPS dump is missing {expected}"
            )))
        }
    }
}

fn field_number(values: &[&str], index: usize, part: &str) -> Result<f64, Error> {
    values
        .get(index)
        .ok_or_else(|| truncated(part))?
        .parse()
        .map_err(|_| Error::InvalidInput("a LAMMPS coordinate is not a number".to_owned()))
}

#[derive(Clone, Copy)]
enum CoordinateMode {
    Fractional,
    Cartesian,
}

fn coordinate_columns(columns: &[&str]) -> Result<(CoordinateMode, usize, usize, usize), Error> {
    let find = |names: &[&str]| {
        names
            .iter()
            .find_map(|name| columns.iter().position(|column| column == name))
    };
    if let (Some(x), Some(y), Some(z)) = (find(&["xs"]), find(&["ys"]), find(&["zs"])) {
        return Ok((CoordinateMode::Fractional, x, y, z));
    }
    if let (Some(x), Some(y), Some(z)) = (find(&["xsu"]), find(&["ysu"]), find(&["zsu"])) {
        return Ok((CoordinateMode::Fractional, x, y, z));
    }
    if let (Some(x), Some(y), Some(z)) =
        (find(&["x", "xu"]), find(&["y", "yu"]), find(&["z", "zu"]))
    {
        return Ok((CoordinateMode::Cartesian, x, y, z));
    }
    Err(Error::InvalidInput(
        "LAMMPS ATOMS section needs x/y/z, xu/yu/zu or xs/ys/zs".to_owned(),
    ))
}

fn lammps_lattice(bounds: [[f64; 3]; 3]) -> Result<Lattice, Error> {
    let [xlo_b, xhi_b, xy] = bounds[0];
    let [ylo_b, yhi_b, xz] = bounds[1];
    let [zlo, zhi, yz] = bounds[2];
    let xlo = xlo_b - 0.0_f64.min(xy).min(xz).min(xy + xz);
    let xhi = xhi_b - 0.0_f64.max(xy).max(xz).max(xy + xz);
    let ylo = ylo_b - 0.0_f64.min(yz);
    let yhi = yhi_b - 0.0_f64.max(yz);
    Lattice::new([
        [xhi - xlo, 0.0, 0.0],
        [xy, yhi - ylo, 0.0],
        [xz, yz, zhi - zlo],
    ])
}

fn truncated(part: &str) -> Error {
    Error::InvalidInput(format!("LAMMPS dump ended before {part}"))
}
