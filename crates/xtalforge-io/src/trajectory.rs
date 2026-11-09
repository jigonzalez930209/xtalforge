//! Trajectory model and multi-frame extended XYZ (Task 3.6).

#![allow(clippy::indexing_slicing)] // fixed-format parser, lengths checked above
use std::collections::BTreeMap;

use serde::{Deserialize, Serialize};
use xtalforge_core::{Error, Lattice, Site, Structure, StructureBuilder};

/// One frame of a trajectory.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct Frame {
    /// Fractional coordinates in atom order.
    pub frac: Vec<[f64; 3]>,
    /// Cell of this frame, when it varies.
    pub lattice: Option<[[f64; 3]; 3]>,
    /// Potential energy in eV, when present.
    pub energy: Option<f64>,
    /// Step number, when present.
    pub step: Option<u64>,
    /// Arbitrary comment-line key/value pairs.
    pub info: BTreeMap<String, String>,
}

/// Many frames sharing atom order and elements.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct Trajectory {
    /// Element symbols in atom order.
    pub elements: Vec<String>,
    /// Fixed cell when every frame uses the same one.
    pub lattice: Option<[[f64; 3]; 3]>,
    /// Frames in file order.
    pub frames: Vec<Frame>,
}

impl Trajectory {
    /// Builds a structure from one frame.
    ///
    /// # Errors
    ///
    /// Returns [`Error::InvalidInput`] for an out-of-range frame.
    pub fn structure(&self, frame: usize) -> Result<Structure, Error> {
        let frame = self.frames.get(frame).ok_or_else(|| {
            Error::InvalidInput(format!("trajectory frame {frame} is out of range"))
        })?;
        let matrix = frame
            .lattice
            .or(self.lattice)
            .ok_or_else(|| Error::InvalidInput("the frame has no lattice".to_owned()))?;
        let lattice = Lattice::new(matrix)?;
        let mut builder = StructureBuilder::new("trajectory", lattice);
        for (index, (element, frac)) in self.elements.iter().zip(frame.frac.iter()).enumerate() {
            builder = builder.site(Site::new(
                format!("{element}{}", index + 1),
                element,
                *frac,
            )?);
        }
        builder.build()
    }
}

/// Reads multi-frame extended XYZ.
///
/// # Errors
///
/// Returns [`Error::InvalidInput`] for a truncated frame or a bad lattice.
pub fn read_trajectory_xyz(text: &str) -> Result<Trajectory, Error> {
    let lines: Vec<&str> = text.lines().collect();
    let mut index = 0;
    let mut elements: Vec<String> = Vec::new();
    let mut frames = Vec::new();
    let mut fixed_lattice: Option<[[f64; 3]; 3]> = None;
    while index < lines.len() {
        if lines[index].trim().is_empty() {
            index += 1;
            continue;
        }
        let count: usize = lines[index]
            .trim()
            .parse()
            .map_err(|_| Error::InvalidInput("bad XYZ atom count".to_owned()))?;
        index += 1;
        let comment = lines.get(index).copied().unwrap_or("");
        index += 1;
        let info = parse_info(comment);
        let lattice = info.get("Lattice").and_then(|value| parse_lattice(value));
        if fixed_lattice.is_none() {
            fixed_lattice = lattice;
        }
        let mut frac = Vec::with_capacity(count);
        let mut frame_elements = Vec::with_capacity(count);
        for _ in 0..count {
            let line = lines.get(index).ok_or_else(|| {
                Error::InvalidInput("trajectory ended before all atoms".to_owned())
            })?;
            index += 1;
            let parts: Vec<&str> = line.split_whitespace().collect();
            if parts.len() < 4 {
                return Err(Error::InvalidInput(
                    "an XYZ atom line needs element and three coordinates".to_owned(),
                ));
            }
            let element = parts[0].to_owned();
            let cart = [
                parts[1]
                    .parse::<f64>()
                    .map_err(|_| Error::InvalidInput("bad XYZ x".to_owned()))?,
                parts[2]
                    .parse::<f64>()
                    .map_err(|_| Error::InvalidInput("bad XYZ y".to_owned()))?,
                parts[3]
                    .parse::<f64>()
                    .map_err(|_| Error::InvalidInput("bad XYZ z".to_owned()))?,
            ];
            frame_elements.push(element);
            frac.push(if let Some(matrix) = lattice {
                Lattice::new(matrix)?.to_frac(cart)
            } else {
                cart
            });
        }
        if elements.is_empty() {
            elements = frame_elements;
        } else if elements != frame_elements {
            return Err(Error::InvalidInput(
                "trajectory frames disagree on atom order or elements".to_owned(),
            ));
        }
        frames.push(Frame {
            frac,
            lattice,
            energy: info.get("energy").and_then(|value| value.parse().ok()),
            step: info.get("step").and_then(|value| value.parse().ok()),
            info,
        });
    }
    Ok(Trajectory {
        elements,
        lattice: fixed_lattice,
        frames,
    })
}

/// Writes a trajectory as extended XYZ.
///
/// # Errors
///
/// Returns [`Error::InvalidInput`] when a frame has no lattice.
pub fn write_trajectory_xyz(trajectory: &Trajectory) -> Result<String, Error> {
    let mut out = String::new();
    for frame in &trajectory.frames {
        let matrix = frame.lattice.or(trajectory.lattice).ok_or_else(|| {
            Error::InvalidInput("a trajectory frame needs a lattice to write XYZ".to_owned())
        })?;
        let lattice = Lattice::new(matrix)?;
        out.push_str(&format!(
            "{}\nLattice=\"{}\"",
            trajectory.elements.len(),
            matrix
                .iter()
                .flatten()
                .map(|value| format!("{value:.10}"))
                .collect::<Vec<_>>()
                .join(" ")
        ));
        if let Some(energy) = frame.energy {
            out.push_str(&format!(" energy={energy:.10}"));
        }
        if let Some(step) = frame.step {
            out.push_str(&format!(" step={step}"));
        }
        out.push('\n');
        for (element, frac) in trajectory.elements.iter().zip(frame.frac.iter()) {
            let cart = lattice.to_cart(*frac);
            out.push_str(&format!(
                "{element} {:.10} {:.10} {:.10}\n",
                cart[0], cart[1], cart[2]
            ));
        }
    }
    Ok(out)
}

pub(super) fn parse_info(comment: &str) -> BTreeMap<String, String> {
    let mut info = BTreeMap::new();
    let mut tokens = Vec::new();
    let mut current = String::new();
    let mut quoted = false;
    for character in comment.chars() {
        match character {
            '"' => {
                quoted = !quoted;
                current.push(character);
            }
            ' ' | '\t' if !quoted => {
                if !current.is_empty() {
                    tokens.push(std::mem::take(&mut current));
                }
            }
            _ => current.push(character),
        }
    }
    if !current.is_empty() {
        tokens.push(current);
    }
    for token in tokens {
        if let Some((key, value)) = token.split_once('=') {
            info.insert(key.to_owned(), value.trim_matches('"').to_owned());
        }
    }
    info
}

pub(super) fn parse_lattice(text: &str) -> Option<[[f64; 3]; 3]> {
    let values: Vec<f64> = text
        .split_whitespace()
        .filter_map(|value| value.parse().ok())
        .collect();
    if values.len() != 9 {
        return None;
    }
    Some([
        [values[0], values[1], values[2]],
        [values[3], values[4], values[5]],
        [values[6], values[7], values[8]],
    ])
}

mod foreign;
mod indexed;

pub use foreign::{read_geo_end, read_lammps_dump, read_lammps_dump_with_types};
pub use indexed::{XyzIndex, read_trajectory, read_xdatcar};
