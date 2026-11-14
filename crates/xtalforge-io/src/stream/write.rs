//! Streaming XYZ and LAMMPS dump writers (Task 17.2.2).
//!
//! Both writers take one [`StreamFrame`] at a time. The iterating helpers in
//! [`XyzWriter`] and [`super::LammpsWriter`] consume a lazy `IntoIterator`, so
//! a trajectory never has to exist as a `Vec` of frames.

use std::borrow::Borrow;
use std::collections::BTreeMap;
use std::io::Write;

use serde_json::Value;
use xtalforge_core::elements::element_by_symbol;
use xtalforge_core::{Error, Lattice};

use super::StreamFrame;

/// Column names that a custom column may not shadow.
const RESERVED: [&str; 5] = ["species", "pos", "forces", "Z", "label"];

/// Writes one frame as extended XYZ.
///
/// The comment line carries `Lattice=` when the frame has a cell, the
/// `Properties=` declaration, `energy=` and `step=` when present, and every
/// other [`StreamFrame::info`] pair in sorted order. Atom lines carry the
/// element, three Cartesian coordinates, an optional `forces` column and the
/// custom columns. Coordinates use ten decimal places, so the output is
/// deterministic.
///
/// # Errors
///
/// Returns [`Error::InvalidInput`] for an empty frame, an atom-count mismatch,
/// a non-finite coordinate, force or energy, an unknown element symbol, a bad
/// lattice, or a custom column with mixed types or a reserved name.
pub fn write_xyz_frame<W: Write>(writer: &mut W, frame: &StreamFrame) -> Result<(), Error> {
    let checked = CheckedFrame::new(frame)?;
    let mut header = String::new();
    if let Some(lattice) = checked.lattice {
        header.push_str(&format!("Lattice=\"{}\"", matrix_text(&lattice)));
    }
    header.push_str(" Properties=species:S:1:pos:R:3");
    if frame.forces.is_some() {
        header.push_str(":forces:R:3");
    }
    for column in &checked.columns {
        header.push_str(&format!(":{}:{}:{}", column.key, column.kind, column.count));
    }
    if let Some(energy) = frame.energy {
        header.push_str(&format!(" energy={energy:.10}"));
    }
    if let Some(step) = frame.step {
        header.push_str(&format!(" step={step}"));
    }
    append_info(&mut header, &frame.info)?;
    writeln!(writer, "{}", checked.symbols.len()).map_err(io_error)?;
    writeln!(writer, "{header}").map_err(io_error)?;
    for (index, symbol) in checked.symbols.iter().enumerate() {
        let frac = frame.frac.get(index).copied().ok_or_else(|| {
            Error::InvalidInput("the frame has fewer coordinates than elements".to_owned())
        })?;
        let cart = checked
            .lattice
            .map_or(frac, |lattice| lattice.to_cart(frac));
        if !cart.iter().all(|value| value.is_finite()) {
            return Err(Error::InvalidInput(
                "the frame has a non-finite Cartesian coordinate".to_owned(),
            ));
        }
        write!(
            writer,
            "{symbol} {:.10} {:.10} {:.10}",
            cart[0], cart[1], cart[2]
        )
        .map_err(io_error)?;
        if let Some(forces) = &frame.forces {
            let force = forces.get(index).copied().ok_or_else(|| {
                Error::InvalidInput("the frame has fewer forces than atoms".to_owned())
            })?;
            write!(
                writer,
                " {:.10} {:.10} {:.10}",
                force[0], force[1], force[2]
            )
            .map_err(io_error)?;
        }
        for column in &checked.columns {
            let values = frame.columns.get(&column.key).ok_or_else(|| {
                Error::InvalidInput(format!("the column {:?} disappeared", column.key))
            })?;
            let value = values.get(index).ok_or_else(|| {
                Error::InvalidInput(format!(
                    "the column {:?} is shorter than the frame",
                    column.key
                ))
            })?;
            write_column_value(writer, value)?;
        }
        writeln!(writer).map_err(io_error)?;
    }
    Ok(())
}

/// Streaming extended-XYZ writer.
pub struct XyzWriter<W: Write> {
    writer: W,
}

impl<W: Write> XyzWriter<W> {
    /// Wraps a writer.
    pub fn new(writer: W) -> Self {
        Self { writer }
    }

    /// Writes one frame.
    ///
    /// # Errors
    ///
    /// Returns the errors of [`write_xyz_frame`].
    pub fn write_frame(&mut self, frame: &StreamFrame) -> Result<(), Error> {
        write_xyz_frame(&mut self.writer, frame)
    }

    /// Writes every frame from a lazy iterator and returns the frame count.
    ///
    /// The iterator may yield owned frames or references; each frame is
    /// written before the next one is requested.
    ///
    /// # Errors
    ///
    /// Returns the first error from a frame or from the writer.
    pub fn write_frames<I, F>(&mut self, frames: I) -> Result<usize, Error>
    where
        I: IntoIterator<Item = F>,
        F: Borrow<StreamFrame>,
    {
        let mut count = 0_usize;
        for frame in frames {
            self.write_frame(frame.borrow())?;
            count += 1;
        }
        Ok(count)
    }

    /// Writes frames that arrive as `Result`s, such as a stream reader.
    ///
    /// # Errors
    ///
    /// Returns the first error from the source or from the writer.
    pub fn write_results<I>(&mut self, frames: I) -> Result<usize, Error>
    where
        I: IntoIterator<Item = Result<StreamFrame, Error>>,
    {
        let mut count = 0_usize;
        for frame in frames {
            self.write_frame(&frame?)?;
            count += 1;
        }
        Ok(count)
    }

    /// Returns the wrapped writer.
    pub fn into_inner(self) -> W {
        self.writer
    }
}

struct ColumnSpec {
    key: String,
    kind: char,
    count: usize,
}

impl ColumnSpec {
    fn new(key: &str, values: &[Value]) -> Result<Self, Error> {
        check_token(key, "column name")?;
        if RESERVED.contains(&key) {
            return Err(Error::InvalidInput(format!(
                "the custom column {key:?} uses a reserved name"
            )));
        }
        let Some(first) = values.first() else {
            return Err(Error::InvalidInput(format!(
                "the column {key:?} has no values"
            )));
        };
        let (kind, count) = column_shape(first)?;
        for value in values {
            if column_shape(value)? != (kind, count) {
                return Err(Error::InvalidInput(format!(
                    "the column {key:?} mixes value types"
                )));
            }
        }
        Ok(Self {
            key: key.to_owned(),
            kind,
            count,
        })
    }
}

pub(super) struct CheckedFrame {
    pub(super) symbols: Vec<&'static str>,
    pub(super) lattice: Option<Lattice>,
    columns: Vec<ColumnSpec>,
}

impl CheckedFrame {
    pub(super) fn new(frame: &StreamFrame) -> Result<Self, Error> {
        let count = frame.elements.len();
        if count == 0 {
            return Err(Error::InvalidInput(
                "a streaming frame needs at least one atom".to_owned(),
            ));
        }
        if frame.frac.len() != count {
            return Err(Error::InvalidInput(format!(
                "the frame has {count} elements but {} coordinates",
                frame.frac.len()
            )));
        }
        if !frame.frac.iter().flatten().all(|value| value.is_finite()) {
            return Err(Error::InvalidInput(
                "the frame has a non-finite coordinate".to_owned(),
            ));
        }
        if frame.energy.is_some_and(|value| !value.is_finite()) {
            return Err(Error::InvalidInput(
                "the frame energy is not finite".to_owned(),
            ));
        }
        let symbols = frame
            .elements
            .iter()
            .map(|element| element_by_symbol(element).map(|element| element.symbol.as_str()))
            .collect::<Result<Vec<_>, _>>()?;
        let lattice = frame.lattice.map(Lattice::new).transpose()?;
        if let Some(forces) = &frame.forces {
            if forces.len() != count {
                return Err(Error::InvalidInput(format!(
                    "the frame has {count} atoms but {} forces",
                    forces.len()
                )));
            }
            if !forces.iter().flatten().all(|value| value.is_finite()) {
                return Err(Error::InvalidInput(
                    "the frame has a non-finite force".to_owned(),
                ));
            }
        }
        let mut columns = Vec::with_capacity(frame.columns.len());
        for (key, values) in &frame.columns {
            if values.len() != count {
                return Err(Error::InvalidInput(format!(
                    "the column {key:?} has {} values for {count} atoms",
                    values.len()
                )));
            }
            columns.push(ColumnSpec::new(key, values)?);
        }
        Ok(Self {
            symbols,
            lattice,
            columns,
        })
    }
}

fn column_shape(value: &Value) -> Result<(char, usize), Error> {
    match value {
        Value::Number(number) if number.is_i64() || number.is_u64() => Ok(('I', 1)),
        Value::Number(number) => {
            if number.as_f64().is_some_and(f64::is_finite) {
                Ok(('R', 1))
            } else {
                Err(Error::InvalidInput(
                    "a custom column has a non-finite number".to_owned(),
                ))
            }
        }
        Value::Bool(_) => Ok(('L', 1)),
        Value::String(text) => {
            check_column_text(text)?;
            Ok(('S', 1))
        }
        Value::Array(values) => {
            if values.is_empty() {
                return Err(Error::InvalidInput(
                    "a custom column value cannot be an empty array".to_owned(),
                ));
            }
            let Some(first) = values.first() else {
                return Err(Error::InvalidInput(
                    "a custom column value cannot be an empty array".to_owned(),
                ));
            };
            let (kind, _) = column_shape(first)?;
            for value in values {
                let (current, width) = column_shape(value)?;
                if current != kind || width != 1 {
                    return Err(Error::InvalidInput(
                        "a custom column array mixes value types".to_owned(),
                    ));
                }
            }
            Ok((kind, values.len()))
        }
        _ => Err(Error::InvalidInput(
            "a custom column needs a number, a string, a logical or an array".to_owned(),
        )),
    }
}

fn append_info(header: &mut String, info: &BTreeMap<String, String>) -> Result<(), Error> {
    for (key, value) in info {
        if matches!(key.as_str(), "Lattice" | "Properties" | "energy" | "step") {
            continue;
        }
        check_token(key, "comment key")?;
        if value.contains('"') || value.chars().any(char::is_control) {
            return Err(Error::InvalidInput(format!(
                "the comment value for {key:?} cannot be quoted"
            )));
        }
        header.push(' ');
        header.push_str(key);
        header.push('=');
        if value.is_empty() || value.chars().any(char::is_whitespace) {
            header.push('"');
            header.push_str(value);
            header.push('"');
        } else {
            header.push_str(value);
        }
    }
    Ok(())
}

fn check_token(text: &str, what: &str) -> Result<(), Error> {
    if text.is_empty()
        || text.chars().any(|character| {
            character.is_whitespace() || matches!(character, ':' | '=' | '"' | '\'')
        })
    {
        return Err(Error::InvalidInput(format!(
            "{what} {text:?} is not a valid extended-XYZ token"
        )));
    }
    Ok(())
}

fn check_column_text(text: &str) -> Result<(), Error> {
    if text.is_empty()
        || text.contains('"')
        || text.contains('\'')
        || text.chars().any(char::is_control)
    {
        return Err(Error::InvalidInput(
            "a text column needs a non-empty value without quotes or control characters".to_owned(),
        ));
    }
    Ok(())
}

fn matrix_text(lattice: &Lattice) -> String {
    lattice
        .matrix()
        .iter()
        .flatten()
        .map(|value| format!("{value:.10}"))
        .collect::<Vec<_>>()
        .join(" ")
}

fn write_column_value<W: Write>(writer: &mut W, value: &Value) -> Result<(), Error> {
    writer.write_all(b" ").map_err(io_error)?;
    if let Value::Array(values) = value {
        for (index, value) in values.iter().enumerate() {
            if index > 0 {
                writer.write_all(b" ").map_err(io_error)?;
            }
            write_scalar(writer, value)?;
        }
        return Ok(());
    }
    write_scalar(writer, value)
}

fn write_scalar<W: Write>(writer: &mut W, value: &Value) -> Result<(), Error> {
    match value {
        Value::Number(number) => write!(writer, "{number}").map_err(io_error),
        Value::Bool(flag) => {
            let text = if *flag { "T" } else { "F" };
            write!(writer, "{text}").map_err(io_error)
        }
        Value::String(text) if text.chars().any(char::is_whitespace) => {
            write!(writer, "\"{text}\"").map_err(io_error)
        }
        Value::String(text) => write!(writer, "{text}").map_err(io_error),
        _ => Err(Error::InvalidInput(
            "a custom column value changed type after validation".to_owned(),
        )),
    }
}

pub(super) fn io_error(error: std::io::Error) -> Error {
    Error::InvalidInput(format!("could not write the stream: {error}"))
}
