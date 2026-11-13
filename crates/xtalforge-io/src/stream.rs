//! Streaming frame readers for very large XYZ and LAMMPS text (Task 15.2.2).
//!
//! The string-based readers in [`crate::xyz`] and [`crate::trajectory`] keep
//! the whole file in memory and, for trajectories, every parsed frame. The
//! readers here borrow a [`BufRead`] and return one [`StreamFrame`] at a time,
//! so peak memory is one frame plus one line, independent of the frame count.
//!
//! The CIF stream is [`crate::cif::read_cif_stream`]; it has a different shape
//! because a CIF is a single structure rather than a frame sequence.

use std::collections::BTreeMap;
use std::io::BufRead;

use xtalforge_core::{Error, Lattice, Site, Structure, StructureBuilder};

use crate::trajectory::{parse_info, parse_lattice};

/// One frame read from a stream.
///
/// `frac` is fractional when `lattice` is present and Cartesian otherwise,
/// matching [`crate::Trajectory`]. Elements travel with the frame because a
/// stream does not assume that every frame names the same atoms.
#[derive(Debug, Clone, PartialEq)]
pub struct StreamFrame {
    /// Element symbols in atom order.
    pub elements: Vec<String>,
    /// Coordinates in the frame's basis.
    pub frac: Vec<[f64; 3]>,
    /// Cell of this frame, when the frame carries one.
    pub lattice: Option<[[f64; 3]; 3]>,
    /// Potential energy in eV, when the comment line carried one.
    pub energy: Option<f64>,
    /// Step number, when present.
    pub step: Option<u64>,
    /// Comment-line key/value pairs.
    pub info: BTreeMap<String, String>,
}

impl StreamFrame {
    /// Number of atoms in the frame.
    #[must_use]
    pub fn len(&self) -> usize {
        self.frac.len()
    }

    /// Whether the frame has no atoms.
    #[must_use]
    pub fn is_empty(&self) -> bool {
        self.frac.is_empty()
    }

    /// Builds a structure from this frame.
    ///
    /// # Errors
    ///
    /// Returns [`Error::InvalidInput`] when the frame has no lattice or the
    /// coordinates do not build sites.
    pub fn structure(&self, name: &str) -> Result<Structure, Error> {
        let matrix = self.lattice.ok_or_else(|| {
            Error::InvalidInput(
                "a streaming frame needs a lattice to become a structure".to_owned(),
            )
        })?;
        let lattice = Lattice::new(matrix)?;
        let mut builder = StructureBuilder::new(name, lattice);
        for (index, (element, frac)) in self.elements.iter().zip(self.frac.iter()).enumerate() {
            builder = builder.site(Site::new(
                format!("{element}{}", index + 1),
                element,
                *frac,
            )?);
        }
        builder.build()
    }
}

/// Streaming extended-XYZ reader.
///
/// The reader accepts the same comment-line `Lattice=`, `energy=` and `step=`
/// keys as [`crate::read_trajectory_xyz`]; extra per-atom property columns are
/// ignored because a frame only carries geometry.
pub struct XyzStream<R: BufRead> {
    reader: R,
    line: String,
}

impl<R: BufRead> XyzStream<R> {
    /// Wraps a buffered reader.
    #[must_use]
    pub fn new(reader: R) -> Self {
        Self {
            reader,
            line: String::new(),
        }
    }

    /// Reads the next frame, or `None` at end of file.
    ///
    /// # Errors
    ///
    /// Returns [`Error::InvalidInput`] for a truncated frame, a bad count, a
    /// malformed atom line or a bad lattice.
    pub fn next_frame(&mut self) -> Result<Option<StreamFrame>, Error> {
        let Some(count_line) = self.next_non_empty()? else {
            return Ok(None);
        };
        let count: usize = count_line
            .trim()
            .parse()
            .map_err(|_| Error::InvalidInput("bad XYZ atom count".to_owned()))?;
        let comment = self.read_line()?.unwrap_or_default();
        let info = parse_info(&comment);
        let lattice = info.get("Lattice").and_then(|value| parse_lattice(value));
        let mut elements = Vec::with_capacity(count);
        let mut frac = Vec::with_capacity(count);
        for _ in 0..count {
            let line = self.read_line()?.ok_or_else(|| {
                Error::InvalidInput("trajectory ended before all atoms".to_owned())
            })?;
            let parts: Vec<&str> = line.split_whitespace().collect();
            if parts.len() < 4 {
                return Err(Error::InvalidInput(
                    "an XYZ atom line needs element and three coordinates".to_owned(),
                ));
            }
            let Some(element) = parts.first() else {
                return Err(Error::InvalidInput("empty XYZ atom line".to_owned()));
            };
            let cart = [
                parse_number(parts.get(1).copied(), "bad XYZ x")?,
                parse_number(parts.get(2).copied(), "bad XYZ y")?,
                parse_number(parts.get(3).copied(), "bad XYZ z")?,
            ];
            elements.push((*element).to_owned());
            frac.push(match lattice {
                Some(matrix) => Lattice::new(matrix)?.to_frac(cart),
                None => cart,
            });
        }
        Ok(Some(StreamFrame {
            elements,
            frac,
            lattice,
            energy: info.get("energy").and_then(|value| value.parse().ok()),
            step: info.get("step").and_then(|value| value.parse().ok()),
            info,
        }))
    }

    fn next_non_empty(&mut self) -> Result<Option<String>, Error> {
        loop {
            match self.read_line()? {
                None => return Ok(None),
                Some(line) if line.trim().is_empty() => continue,
                Some(line) => return Ok(Some(line)),
            }
        }
    }

    fn read_line(&mut self) -> Result<Option<String>, Error> {
        read_line(&mut self.reader, &mut self.line)
    }
}

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

fn read_line<R: BufRead>(reader: &mut R, buffer: &mut String) -> Result<Option<String>, Error> {
    buffer.clear();
    let read = reader
        .read_line(buffer)
        .map_err(|error| Error::InvalidInput(format!("could not read the stream: {error}")))?;
    if read == 0 {
        return Ok(None);
    }
    Ok(Some(buffer.trim_end_matches(['\n', '\r']).to_owned()))
}

fn parse_number(value: Option<&str>, label: &str) -> Result<f64, Error> {
    value
        .ok_or_else(|| truncated("coordinate column"))?
        .parse()
        .map_err(|_| Error::InvalidInput(label.to_owned()))
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
