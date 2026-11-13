//! Trajectory recording for optimizers, scans and dynamics (Task 14.2.2).
//!
//! A [`Trajectory`] is a list of [`TrajectoryFrame`] records.  Each frame keeps
//! the potential energy, the optimization objective (`E + P V`), the largest
//! force, the cell and, when requested, the fractional coordinates.  Frames are
//! plain serde data: the same object is written as one JSON document or as
//! newline-delimited JSON (`.jsonl`), so a long run can be streamed and read
//! back by any tool that understands JSON.
//!
//! Recording is off by default; [`TrajectoryOptions::default`] has
//! `every = 0`.  This keeps the memory of existing callers unchanged: a
//! trajectory is opt-in.
//!
//! **Documented precision.**  A JSON round trip preserves every number that
//! `serde_json`'s parser reconstructs exactly; the parser is not required to
//! be bit-exact, so a frame read back can differ from the in-memory one in the
//! last unit in the last place.  Use the in-memory frames when exact bits
//! matter.

use std::fs;
use std::path::Path;

use serde::{Deserialize, Serialize};
use xtalforge_core::Lattice;

use crate::Error;

/// How often a run records a trajectory frame.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct TrajectoryOptions {
    /// Record every `every` steps; `0` disables recording.
    pub every: usize,
    /// Keep fractional coordinates in every frame.
    ///
    /// Energies alone cost one record per frame; coordinates make a frame
    /// grow with the atom count.
    pub positions: bool,
}

impl TrajectoryOptions {
    /// Records every step and keeps the coordinates.
    #[must_use]
    pub fn full() -> Self {
        Self {
            every: 1,
            positions: true,
        }
    }

    /// Records energies every `every` steps, without coordinates.
    ///
    /// # Errors
    ///
    /// Returns [`Error::InvalidInput`] when `every` is zero.
    pub fn every(every: usize) -> Result<Self, Error> {
        if every == 0 {
            return Err(Error::InvalidInput(
                "a trajectory cadence must be at least one step".to_owned(),
            ));
        }
        Ok(Self {
            every,
            positions: false,
        })
    }

    /// Whether recording is enabled.
    #[must_use]
    pub fn enabled(&self) -> bool {
        self.every != 0
    }

    /// Whether a frame is due at `step`.
    #[must_use]
    pub fn wants(&self, step: usize) -> bool {
        self.enabled() && step.is_multiple_of(self.every)
    }

    /// Same cadence, with or without coordinates.
    #[must_use]
    pub fn with_positions(mut self, positions: bool) -> Self {
        self.positions = positions;
        self
    }
}

/// One recorded frame.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct TrajectoryFrame {
    /// Step number within the run.
    pub step: usize,
    /// Scan point index when the frame comes from a relaxed-scan point.
    #[serde(default)]
    pub point: Option<usize>,
    /// Simulated time in femtoseconds; zero for optimizers and scans.
    pub time_fs: f64,
    /// Potential energy in eV.
    pub energy: f64,
    /// Optimization objective in eV, `energy + pressure * volume`.
    pub enthalpy: f64,
    /// Largest force on a free site in eV/Angstrom; zero when not computed.
    pub max_force: f64,
    /// Cell volume in Angstrom^3.
    pub volume: f64,
    /// Cell at this frame.
    pub lattice: Lattice,
    /// Fractional coordinates, present when the cadence asked for them.
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub positions: Option<Vec<[f64; 3]>>,
}

/// Frames recorded by one run, in order.
#[derive(Debug, Clone, Default, PartialEq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct Trajectory {
    /// Recorded frames.
    pub frames: Vec<TrajectoryFrame>,
}

impl Trajectory {
    /// An empty trajectory.
    #[must_use]
    pub fn new() -> Self {
        Self::default()
    }

    /// Number of frames.
    #[must_use]
    pub fn len(&self) -> usize {
        self.frames.len()
    }

    /// Whether no frame was recorded.
    #[must_use]
    pub fn is_empty(&self) -> bool {
        self.frames.is_empty()
    }

    /// Appends one frame.
    pub fn push(&mut self, frame: TrajectoryFrame) {
        self.frames.push(frame);
    }

    /// Energies in frame order, in eV.
    #[must_use]
    pub fn energies(&self) -> Vec<f64> {
        self.frames.iter().map(|frame| frame.energy).collect()
    }

    /// Serializes the whole trajectory as one JSON document.
    ///
    /// # Errors
    ///
    /// Returns [`Error::Json`] when serialization fails.
    pub fn to_json(&self) -> Result<String, Error> {
        serde_json::to_string(self).map_err(|error| Error::Json(error.to_string()))
    }

    /// Parses a trajectory written by [`Trajectory::to_json`].
    ///
    /// # Errors
    ///
    /// Returns [`Error::Json`] for invalid JSON or unknown fields.
    pub fn from_json(text: &str) -> Result<Self, Error> {
        serde_json::from_str(text).map_err(|error| Error::Json(error.to_string()))
    }

    /// Serializes one frame per line, newline terminated.
    ///
    /// # Errors
    ///
    /// Returns [`Error::Json`] when a frame cannot be serialized.
    pub fn to_jsonl(&self) -> Result<String, Error> {
        let mut text = String::new();
        for frame in &self.frames {
            let line =
                serde_json::to_string(frame).map_err(|error| Error::Json(error.to_string()))?;
            text.push_str(&line);
            text.push('\n');
        }
        Ok(text)
    }

    /// Parses newline-delimited frames; blank lines are skipped.
    ///
    /// # Errors
    ///
    /// Returns [`Error::Json`] for a malformed line.
    pub fn from_jsonl(text: &str) -> Result<Self, Error> {
        let mut frames = Vec::new();
        for (index, line) in text.lines().enumerate() {
            if line.trim().is_empty() {
                continue;
            }
            let frame = serde_json::from_str(line)
                .map_err(|error| Error::Json(format!("trajectory line {}: {error}", index + 1)))?;
            frames.push(frame);
        }
        Ok(Self { frames })
    }

    /// Writes one JSON document, or JSONL when the path ends in `.jsonl`.
    ///
    /// # Errors
    ///
    /// Returns [`Error::InvalidInput`] for a filesystem failure and
    /// [`Error::Json`] for a serialization failure.
    pub fn write(&self, path: impl AsRef<Path>) -> Result<(), Error> {
        let path = path.as_ref();
        let text = if is_jsonl(path) {
            self.to_jsonl()?
        } else {
            self.to_json()?
        };
        write_text(path, &text)
    }

    /// Reads a file written by [`Trajectory::write`].
    ///
    /// # Errors
    ///
    /// Returns [`Error::InvalidInput`] for a filesystem failure and
    /// [`Error::Json`] for malformed content.
    pub fn read(path: impl AsRef<Path>) -> Result<Self, Error> {
        let path = path.as_ref();
        let text = fs::read_to_string(path).map_err(|error| {
            Error::InvalidInput(format!("trajectory {}: {error}", path.display()))
        })?;
        if is_jsonl(path) {
            Self::from_jsonl(&text)
        } else {
            Self::from_json(&text)
        }
    }
}

/// Whether a path asks for newline-delimited JSON.
pub(crate) fn is_jsonl(path: &Path) -> bool {
    path.extension()
        .and_then(std::ffi::OsStr::to_str)
        .is_some_and(|extension| extension.eq_ignore_ascii_case("jsonl"))
}

/// Writes text through a temporary file so a crash never leaves half a file.
pub(crate) fn write_text(path: &Path, text: &str) -> Result<(), Error> {
    let temporary = path.with_extension(format!(
        "{}.tmp",
        path.extension()
            .and_then(std::ffi::OsStr::to_str)
            .unwrap_or("tmp")
    ));
    fs::write(&temporary, text)
        .map_err(|error| Error::InvalidInput(format!("{}: {error}", temporary.display())))?;
    fs::rename(&temporary, path)
        .map_err(|error| Error::InvalidInput(format!("{}: {error}", path.display())))
}
