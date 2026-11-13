//! Resumable scan checkpoints (Task 14.2.2).
//!
//! A relaxed scan can write a [`ScanCheckpoint`] after every finished point.
//! Re-running the same scan with the same checkpoint path resumes at the first
//! point that was not recorded, seeded from the last relaxed structure, so a
//! cancelled or interrupted scan never repeats its expensive points.
//!
//! The checkpoint is a plain JSON document; a path ending in `.jsonl` uses a
//! header line followed by one JSON object per point.  Both forms keep the
//! coordinate definition, the requested values and the atom count, and a
//! resumed scan refuses a checkpoint that describes a different scan instead of
//! silently mixing results.
//!
//! Files are written through a temporary file and renamed, so an interrupted
//! write leaves the previous checkpoint intact.  The whole checkpoint is
//! rewritten after every finished point; that is simple and safe, and a
//! `.jsonl` file can also be read incrementally by another tool while the scan
//! runs.
//!
//! **Documented precision.**  A resumed point starts from the structure read
//! back from disk.  `serde_json`'s parser is not required to reconstruct every
//! decimal spelling bit-exactly, so a resumed scan can differ from an
//! uninterrupted one in the last unit in the last place; energies and
//! geometries agree to the tolerances a force field itself works at.

use std::path::Path;

use serde::{Deserialize, Serialize};

use crate::Error;
use crate::scan::coordinates::ScanCoordinate;
use crate::scan::driver::ScanPoint;
use crate::trajectory::{is_jsonl, write_text};

/// Marker of the JSONL header.
const FORMAT: &str = "xtalforge-scan-checkpoint";
/// Header version of the JSONL form.
const VERSION: u32 = 1;

/// Everything a scan needs to resume.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct ScanCheckpoint {
    /// Coordinate that was scanned.
    pub coordinate: ScanCoordinate,
    /// Requested values, in scan order.
    pub values: Vec<f64>,
    /// Number of sites the scan started from.
    pub n_sites: usize,
    /// Finished points, in scan order.
    pub points: Vec<ScanPoint>,
}

impl ScanCheckpoint {
    /// An empty checkpoint for one scan.
    #[must_use]
    pub fn new(coordinate: ScanCoordinate, values: Vec<f64>, n_sites: usize) -> Self {
        Self {
            coordinate,
            values,
            n_sites,
            points: Vec::new(),
        }
    }

    /// Whether every requested value has a finished point.
    #[must_use]
    pub fn is_complete(&self) -> bool {
        !self.values.is_empty() && self.points.len() >= self.values.len()
    }

    /// How many points are already finished.
    #[must_use]
    pub fn resumed_len(&self) -> usize {
        self.points.len().min(self.values.len())
    }

    /// Checks that this checkpoint belongs to the scan being started.
    ///
    /// # Errors
    ///
    /// Returns [`Error::InvalidInput`] when the coordinate, the values, the
    /// atom count or the number of points differ.
    pub fn matches(
        &self,
        coordinate: &ScanCoordinate,
        values: &[f64],
        n_sites: usize,
    ) -> Result<(), Error> {
        if &self.coordinate != coordinate {
            return Err(Error::InvalidInput(format!(
                "the checkpoint scans {} but this run scans {}",
                self.coordinate.label(),
                coordinate.label()
            )));
        }
        if self.values != values {
            return Err(Error::InvalidInput(
                "the checkpoint was written for different scan values".to_owned(),
            ));
        }
        if self.n_sites != n_sites {
            return Err(Error::InvalidInput(format!(
                "the checkpoint has {} sites but the structure has {n_sites}",
                self.n_sites
            )));
        }
        if self.points.len() > self.values.len() {
            return Err(Error::InvalidInput(format!(
                "the checkpoint has {} points for {} values",
                self.points.len(),
                self.values.len()
            )));
        }
        Ok(())
    }

    /// Serializes the checkpoint as one JSON document.
    ///
    /// # Errors
    ///
    /// Returns [`Error::Json`] when serialization fails.
    pub fn to_json(&self) -> Result<String, Error> {
        serde_json::to_string(self).map_err(|error| Error::Json(error.to_string()))
    }

    /// Parses a checkpoint written by [`ScanCheckpoint::to_json`].
    ///
    /// # Errors
    ///
    /// Returns [`Error::Json`] for invalid JSON or unknown fields.
    pub fn from_json(text: &str) -> Result<Self, Error> {
        serde_json::from_str(text).map_err(|error| Error::Json(error.to_string()))
    }

    /// Serializes a header line plus one line per point.
    ///
    /// # Errors
    ///
    /// Returns [`Error::Json`] when a record cannot be serialized.
    pub fn to_jsonl(&self) -> Result<String, Error> {
        let header = Header {
            format: FORMAT.to_owned(),
            version: VERSION,
            coordinate: self.coordinate.clone(),
            values: self.values.clone(),
            n_sites: self.n_sites,
        };
        let mut text =
            serde_json::to_string(&header).map_err(|error| Error::Json(error.to_string()))?;
        text.push('\n');
        for point in &self.points {
            let line =
                serde_json::to_string(point).map_err(|error| Error::Json(error.to_string()))?;
            text.push_str(&line);
            text.push('\n');
        }
        Ok(text)
    }

    /// Parses a checkpoint written by [`ScanCheckpoint::to_jsonl`].
    ///
    /// # Errors
    ///
    /// Returns [`Error::Json`] for a malformed header, an unsupported version
    /// or a malformed point.
    pub fn from_jsonl(text: &str) -> Result<Self, Error> {
        let mut lines = text.lines().filter(|line| !line.trim().is_empty());
        let Some(first) = lines.next() else {
            return Err(Error::Json("scan checkpoint is empty".to_owned()));
        };
        let header: Header = serde_json::from_str(first)
            .map_err(|error| Error::Json(format!("checkpoint header: {error}")))?;
        if header.format != FORMAT {
            return Err(Error::Json(format!(
                "checkpoint header names {:?}, expected {FORMAT:?}",
                header.format
            )));
        }
        if header.version != VERSION {
            return Err(Error::Json(format!(
                "checkpoint version {} is not supported",
                header.version
            )));
        }
        let mut points = Vec::new();
        for (index, line) in lines.enumerate() {
            let point = serde_json::from_str(line)
                .map_err(|error| Error::Json(format!("checkpoint point {}: {error}", index + 1)))?;
            points.push(point);
        }
        Ok(Self {
            coordinate: header.coordinate,
            values: header.values,
            n_sites: header.n_sites,
            points,
        })
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

    /// Reads a checkpoint written by [`ScanCheckpoint::write`].
    ///
    /// # Errors
    ///
    /// Returns [`Error::InvalidInput`] for a filesystem failure and
    /// [`Error::Json`] for malformed content.
    pub fn read(path: impl AsRef<Path>) -> Result<Self, Error> {
        let path = path.as_ref();
        let text = std::fs::read_to_string(path).map_err(|error| {
            Error::InvalidInput(format!("scan checkpoint {}: {error}", path.display()))
        })?;
        if is_jsonl(path) {
            Self::from_jsonl(&text)
        } else {
            Self::from_json(&text)
        }
    }
}

/// First line of the JSONL form.
#[derive(Debug, Clone, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
struct Header {
    format: String,
    version: u32,
    coordinate: ScanCoordinate,
    values: Vec<f64>,
    n_sites: usize,
}
