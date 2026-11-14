//! Streaming LAMMPS dump writer (Task 17.2.2).

use std::borrow::Borrow;
use std::io::Write;

use xtalforge_core::Error;

use super::StreamFrame;
use super::write::{CheckedFrame, io_error};

/// Writes one frame as a LAMMPS text dump.
///
/// The periodic box is written in the restricted triclinic form LAMMPS
/// expects (`xy xz yz`), the timestep comes from [`StreamFrame::step`] or is
/// zero when the frame has none, and the atom rows carry `id`, `element` and
/// Cartesian `x`, `y`, `z`. The output is deterministic.
///
/// # Errors
///
/// Returns [`Error::InvalidInput`] for an empty frame, an atom-count mismatch,
/// a non-finite coordinate, an unknown element symbol, a missing lattice or a
/// cell that cannot become a positive LAMMPS box.
pub fn write_lammps_frame<W: Write>(writer: &mut W, frame: &StreamFrame) -> Result<(), Error> {
    let checked = CheckedFrame::new(frame)?;
    let Some(lattice) = checked.lattice else {
        return Err(Error::InvalidInput(
            "a LAMMPS dump needs a lattice".to_owned(),
        ));
    };
    let cell = lattice.standard_orientation()?;
    let [[ax, _, _], [bx, by, _], [cx, cy, cz]] = cell.matrix();
    if ax <= 0.0 || by <= 0.0 || cz <= 0.0 {
        return Err(Error::InvalidInput(
            "the LAMMPS box needs positive x, y and z lengths".to_owned(),
        ));
    }
    // A dump reports the bounding box, not the cell edges: include the tilt
    // contributions exactly the way the reader removes them.
    let xlo = 0.0_f64.min(bx).min(cx).min(bx + cx);
    let xhi = ax + 0.0_f64.max(bx).max(cx).max(bx + cx);
    let ylo = 0.0_f64.min(cy);
    let yhi = by + 0.0_f64.max(cy);
    let step = frame.step.unwrap_or(0);
    let count = checked.symbols.len();
    writer.write_all(b"ITEM: TIMESTEP\n").map_err(io_error)?;
    writeln!(writer, "{step}").map_err(io_error)?;
    writer
        .write_all(b"ITEM: NUMBER OF ATOMS\n")
        .map_err(io_error)?;
    writeln!(writer, "{count}").map_err(io_error)?;
    writer
        .write_all(b"ITEM: BOX BOUNDS xy xz yz pp pp pp\n")
        .map_err(io_error)?;
    writeln!(writer, "{xlo:.10} {xhi:.10} {bx:.10}").map_err(io_error)?;
    writeln!(writer, "{ylo:.10} {yhi:.10} {cx:.10}").map_err(io_error)?;
    writeln!(writer, "0.0 {cz:.10} {cy:.10}").map_err(io_error)?;
    writer
        .write_all(b"ITEM: ATOMS id element x y z\n")
        .map_err(io_error)?;
    for (index, symbol) in checked.symbols.iter().enumerate() {
        let frac = frame.frac.get(index).copied().ok_or_else(|| {
            Error::InvalidInput("the frame has fewer coordinates than elements".to_owned())
        })?;
        let cart = cell.to_cart(frac);
        if !cart.iter().all(|value| value.is_finite()) {
            return Err(Error::InvalidInput(
                "the frame has a non-finite Cartesian coordinate".to_owned(),
            ));
        }
        writeln!(
            writer,
            "{} {symbol} {:.10} {:.10} {:.10}",
            index + 1,
            cart[0],
            cart[1],
            cart[2]
        )
        .map_err(io_error)?;
    }
    Ok(())
}

/// Streaming LAMMPS dump writer.
pub struct LammpsWriter<W: Write> {
    writer: W,
}

impl<W: Write> LammpsWriter<W> {
    /// Wraps a writer.
    pub fn new(writer: W) -> Self {
        Self { writer }
    }

    /// Writes one frame.
    ///
    /// # Errors
    ///
    /// Returns the errors of [`write_lammps_frame`].
    pub fn write_frame(&mut self, frame: &StreamFrame) -> Result<(), Error> {
        write_lammps_frame(&mut self.writer, frame)
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
