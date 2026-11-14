//! Streaming frame readers and writers for very large XYZ and LAMMPS text
//! (Tasks 15.2.2 and 17.2.2).
//!
//! The string-based readers in [`crate::xyz`] and [`crate::trajectory`] keep
//! the whole file in memory and, for trajectories, every parsed frame. The
//! readers here borrow a [`BufRead`] and return one [`StreamFrame`] at a time,
//! so peak memory is one frame plus one line, independent of the frame count.
//!
//! The writers in this module mirror that contract: [`XyzWriter`] and
//! [`LammpsWriter`] accept one frame at a time and write it straight to a
//! [`std::io::Write`], so a trajectory iterator never becomes a `Vec` of
//! frames. A frame carries element symbols, coordinates, an optional cell,
//! energy, step, forces and arbitrary custom per-atom columns.
//!
//! The CIF stream is [`crate::cif::read_cif_stream`]; it has a different shape
//! because a CIF is a single structure rather than a frame sequence. Its
//! writer is [`crate::cif::write_cif_stream`].

use std::collections::BTreeMap;
use std::io::BufRead;

use serde_json::Value;
use xtalforge_core::{Error, Lattice, Site, Structure, StructureBuilder};

use crate::trajectory::{parse_info, parse_lattice};
use crate::xyz::properties::{parse_atom, parse_properties, split_fields};

/// Per-atom extended-XYZ columns accumulated while reading a frame.
type Columns = BTreeMap<String, Vec<Value>>;

/// Forces and custom columns separated from the raw property map.
struct ParsedColumns {
    forces: Option<Vec<[f64; 3]>>,
    columns: Columns,
}

mod lammps;
mod lammps_write;
mod write;

pub use lammps::LammpsStream;
pub use lammps_write::{LammpsWriter, write_lammps_frame};
pub use write::{XyzWriter, write_xyz_frame};

/// One frame read from a stream or handed to a writer.
///
/// `frac` is fractional when `lattice` is present and Cartesian otherwise,
/// matching [`crate::Trajectory`]. Elements travel with the frame because a
/// stream does not assume that every frame names the same atoms.
///
/// The default value is an empty frame; it exists so callers can build a
/// frame with `..StreamFrame::default()` when only some fields are set.
#[derive(Debug, Clone, PartialEq, Default)]
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
    /// Forces in eV/Angstrom in atom order, when present.
    ///
    /// A reader fills this from a `forces` column, or from separate `fx`,
    /// `fy` and `fz` columns. A writer emits a single `forces:R:3` column.
    pub forces: Option<Vec<[f64; 3]>>,
    /// Extra per-atom extended-XYZ columns, excluding `species`, `pos` and
    /// `forces`.
    ///
    /// A scalar value writes one column; an array writes several values under
    /// one declared name. Every vector has one entry per atom.
    pub columns: BTreeMap<String, Vec<Value>>,
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
/// keys as [`crate::read_trajectory_xyz`], and reads the `Properties=`
/// declaration: a `forces` column or separate `fx`/`fy`/`fz` columns become
/// [`StreamFrame::forces`], and every other non-position column becomes
/// [`StreamFrame::columns`]. Without a `Properties=` declaration the reader
/// takes the historical element plus three Cartesian coordinates layout and
/// ignores any extra atom-line fields.
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
    /// malformed atom line, a bad lattice or a column that does not match its
    /// `Properties=` declaration.
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
        let spec = parse_properties(&comment)?;
        let mut elements = Vec::with_capacity(count);
        let mut frac = Vec::with_capacity(count);
        let mut columns: Columns = BTreeMap::new();
        for _ in 0..count {
            let line = self.read_line()?.ok_or_else(|| {
                Error::InvalidInput("trajectory ended before all atoms".to_owned())
            })?;
            let (element, cart, props) = match &spec {
                Some(spec) => {
                    let fields = split_fields(&line);
                    let record = parse_atom(&fields, spec)?;
                    (record.element, record.cart, Some(record.props))
                }
                None => {
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
                    ((*element).to_owned(), cart, None)
                }
            };
            if let Some(props) = props {
                for (key, value) in props {
                    columns.entry(key).or_default().push(value);
                }
            }
            elements.push(element);
            frac.push(match lattice {
                Some(matrix) => Lattice::new(matrix)?.to_frac(cart),
                None => cart,
            });
        }
        let parsed = split_forces(columns)?;
        Ok(Some(StreamFrame {
            elements,
            frac,
            lattice,
            energy: info.get("energy").and_then(|value| value.parse().ok()),
            step: info.get("step").and_then(|value| value.parse().ok()),
            info,
            forces: parsed.forces,
            columns: parsed.columns,
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

/// Moves a `forces` column or separate `fx`/`fy`/`fz` columns into the
/// dedicated field.
fn split_forces(mut columns: Columns) -> Result<ParsedColumns, Error> {
    if let Some(values) = columns.remove("forces") {
        return Ok(ParsedColumns {
            forces: Some(force_vectors(&values)?),
            columns,
        });
    }
    if ["fx", "fy", "fz"]
        .iter()
        .all(|key| columns.contains_key(*key))
    {
        let fx = columns.remove("fx").unwrap_or_default();
        let fy = columns.remove("fy").unwrap_or_default();
        let fz = columns.remove("fz").unwrap_or_default();
        if fx.len() != fy.len() || fx.len() != fz.len() {
            return Err(Error::InvalidInput(
                "the fx/fy/fz columns disagree on atom count".to_owned(),
            ));
        }
        let mut forces = Vec::with_capacity(fx.len());
        for index in 0..fx.len() {
            forces.push([
                column_real(fx.get(index), "fx")?,
                column_real(fy.get(index), "fy")?,
                column_real(fz.get(index), "fz")?,
            ]);
        }
        return Ok(ParsedColumns {
            forces: Some(forces),
            columns,
        });
    }
    Ok(ParsedColumns {
        forces: None,
        columns,
    })
}

fn force_vectors(values: &[Value]) -> Result<Vec<[f64; 3]>, Error> {
    let mut forces = Vec::with_capacity(values.len());
    for value in values {
        let Some(items) = value.as_array() else {
            return Err(Error::InvalidInput(
                "the forces column needs three numbers per atom".to_owned(),
            ));
        };
        if items.len() != 3 {
            return Err(Error::InvalidInput(
                "the forces column needs three numbers per atom".to_owned(),
            ));
        }
        forces.push([
            column_real(items.first(), "forces")?,
            column_real(items.get(1), "forces")?,
            column_real(items.get(2), "forces")?,
        ]);
    }
    Ok(forces)
}

fn column_real(value: Option<&Value>, name: &str) -> Result<f64, Error> {
    value
        .and_then(Value::as_f64)
        .filter(|value| value.is_finite())
        .ok_or_else(|| Error::InvalidInput(format!("the {name} column needs finite numbers")))
}

fn parse_number(value: Option<&str>, label: &str) -> Result<f64, Error> {
    value
        .ok_or_else(|| Error::InvalidInput("trajectory ended before all atoms".to_owned()))?
        .parse()
        .map_err(|_| Error::InvalidInput(label.to_owned()))
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
