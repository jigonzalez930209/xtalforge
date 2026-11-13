//! Lazy XYZ indexing and the VASP XDATCAR reader (Task 3.6).

use std::collections::BTreeMap;

use xtalforge_core::{Error, Lattice, Structure};

use super::foreign::{read_geo_end, read_lammps_dump};
use super::{Frame, Trajectory, parse_info, parse_lattice, read_trajectory_xyz};

/// Lazily indexed multi-frame extended XYZ.
///
/// The index records where every frame starts and parses one frame on demand,
/// so a 10,000-frame file does not materialise 10,000 `Frame` values. The
/// text itself is borrowed; file-backed random access belongs to the caller.
#[derive(Debug, Clone)]
pub struct XyzIndex<'a> {
    text: &'a str,
    frames: Vec<FrameIndex>,
    elements: Vec<String>,
    lattice: Option<[[f64; 3]; 3]>,
}

#[derive(Debug, Clone)]
struct FrameIndex {
    start: usize,
    energy: Option<f64>,
    step: Option<u64>,
}

impl<'a> XyzIndex<'a> {
    /// Scans the frame headers and atom elements without parsing coordinates.
    ///
    /// # Errors
    ///
    /// Returns [`Error::InvalidInput`] for a truncated frame, a bad count or
    /// frames that disagree on atom order.
    pub fn new(text: &'a str) -> Result<Self, Error> {
        let mut offsets = Vec::new();
        let mut byte = 0;
        for line in text.split_inclusive('\n') {
            offsets.push((byte, line.trim_end_matches(['\n', '\r'])));
            byte += line.len();
        }
        let mut index = 0;
        let mut frames = Vec::new();
        let mut elements: Vec<String> = Vec::new();
        let mut lattice = None;
        while index < offsets.len() {
            if offsets[index].1.trim().is_empty() {
                index += 1;
                continue;
            }
            let start = offsets[index].0;
            let count: usize = offsets[index]
                .1
                .trim()
                .parse()
                .map_err(|_| Error::InvalidInput("bad XYZ atom count".to_owned()))?;
            index += 1;
            let comment = offsets.get(index).map_or("", |(_, line)| *line);
            index += 1;
            let info = parse_info(comment);
            let frame_lattice = info.get("Lattice").and_then(|value| parse_lattice(value));
            if lattice.is_none() {
                lattice = frame_lattice;
            }
            let mut frame_elements = Vec::with_capacity(count);
            for _ in 0..count {
                let line = offsets.get(index).ok_or_else(|| {
                    Error::InvalidInput("trajectory ended before all atoms".to_owned())
                })?;
                index += 1;
                let element = line
                    .1
                    .split_whitespace()
                    .next()
                    .ok_or_else(|| Error::InvalidInput("empty XYZ atom line".to_owned()))?;
                frame_elements.push(element.to_owned());
            }
            if elements.is_empty() {
                elements = frame_elements;
            } else if elements != frame_elements {
                return Err(Error::InvalidInput(
                    "trajectory frames disagree on atom order or elements".to_owned(),
                ));
            }
            frames.push(FrameIndex {
                start,
                energy: info.get("energy").and_then(|value| value.parse().ok()),
                step: info.get("step").and_then(|value| value.parse().ok()),
            });
        }
        Ok(Self {
            text,
            frames,
            elements,
            lattice,
        })
    }

    /// Number of frames.
    #[must_use]
    pub fn len(&self) -> usize {
        self.frames.len()
    }

    /// Whether there are no frames.
    #[must_use]
    pub fn is_empty(&self) -> bool {
        self.frames.is_empty()
    }

    /// Element symbols in atom order.
    #[must_use]
    pub fn elements(&self) -> &[String] {
        &self.elements
    }

    /// Fixed lattice when every indexed frame carries the same one.
    #[must_use]
    pub fn lattice(&self) -> Option<[[f64; 3]; 3]> {
        self.lattice
    }

    /// Energy of every frame, when the comment line carried one.
    #[must_use]
    pub fn energies(&self) -> Vec<Option<f64>> {
        self.frames.iter().map(|frame| frame.energy).collect()
    }

    /// Step numbers of every frame, when present.
    #[must_use]
    pub fn steps(&self) -> Vec<Option<u64>> {
        self.frames.iter().map(|frame| frame.step).collect()
    }

    /// Parses one frame on demand.
    ///
    /// # Errors
    ///
    /// Returns [`Error::InvalidInput`] for an out-of-range frame or bad
    /// coordinates.
    pub fn frame(&self, index: usize) -> Result<Frame, Error> {
        let start = self
            .frames
            .get(index)
            .map(|frame| frame.start)
            .ok_or_else(|| {
                Error::InvalidInput(format!("trajectory frame {index} is out of range"))
            })?;
        let end = self
            .frames
            .get(index + 1)
            .map_or(self.text.len(), |frame| frame.start);
        let one = read_trajectory_xyz(&self.text[start..end])?;
        one.frames
            .into_iter()
            .next()
            .ok_or_else(|| Error::InvalidInput(format!("trajectory frame {index} has no atoms")))
    }

    /// Builds a structure from one frame on demand.
    ///
    /// # Errors
    ///
    /// Returns the errors of [`XyzIndex::frame`] or [`Trajectory::structure`].
    pub fn structure(&self, index: usize) -> Result<Structure, Error> {
        let trajectory = Trajectory {
            elements: self.elements.clone(),
            lattice: self.lattice,
            frames: vec![self.frame(index)?],
        };
        trajectory.structure(0)
    }
}

/// Reads a VASP XDATCAR trajectory.
///
/// Supports the fixed-cell `Direct` and `Cartesian` forms written by VASP and
/// ASE. The first frame may omit the `Direct configuration=` separator.
///
/// # Errors
///
/// Returns [`Error::InvalidInput`] for a truncated header, an unsupported
/// negative scale or a frame with too few coordinates.
pub fn read_xdatcar(text: &str) -> Result<Trajectory, Error> {
    let lines: Vec<&str> = text.lines().collect();
    let mut index = 0;
    while index < lines.len() && lines[index].trim().is_empty() {
        index += 1;
    }
    if index + 6 >= lines.len() {
        return Err(Error::InvalidInput(
            "XDATCAR header is truncated".to_owned(),
        ));
    }
    index += 1; // comment
    let scale: f64 = lines[index]
        .trim()
        .parse()
        .map_err(|_| Error::InvalidInput("XDATCAR has a bad scale".to_owned()))?;
    if !scale.is_finite() || scale <= 0.0 {
        return Err(Error::InvalidInput(
            "XDATCAR needs a positive finite scale".to_owned(),
        ));
    }
    index += 1;
    let mut matrix = [[0.0; 3]; 3];
    for row in &mut matrix {
        let values: Vec<f64> = lines[index]
            .split_whitespace()
            .map(|value| value.parse::<f64>())
            .collect::<Result<_, _>>()
            .map_err(|_| Error::InvalidInput("XDATCAR has a bad lattice row".to_owned()))?;
        if values.len() < 3 {
            return Err(Error::InvalidInput(
                "XDATCAR lattice row needs three values".to_owned(),
            ));
        }
        *row = [values[0] * scale, values[1] * scale, values[2] * scale];
        index += 1;
    }
    let species: Vec<&str> = lines[index].split_whitespace().collect();
    index += 1;
    let counts: Vec<usize> = lines[index]
        .split_whitespace()
        .map(|value| value.parse::<usize>())
        .collect::<Result<_, _>>()
        .map_err(|_| Error::InvalidInput("XDATCAR has bad species counts".to_owned()))?;
    index += 1;
    if species.len() != counts.len() || species.is_empty() {
        return Err(Error::InvalidInput(
            "XDATCAR species and counts disagree".to_owned(),
        ));
    }
    let mut elements = Vec::new();
    for (element, count) in species.iter().zip(counts.iter()) {
        for _ in 0..*count {
            elements.push((*element).to_owned());
        }
    }
    let lattice = Lattice::new(matrix)?;
    let mut frames = Vec::new();
    while index < lines.len() {
        let line = lines[index].trim();
        if line.is_empty() {
            index += 1;
            continue;
        }
        let cartesian = if line.starts_with("Direct configuration")
            || line.starts_with("Cartesian configuration")
        {
            let cartesian = line.starts_with("Cartesian");
            index += 1;
            cartesian
        } else {
            false
        };
        let mut frac = Vec::with_capacity(elements.len());
        for _ in 0..elements.len() {
            let line = lines.get(index).ok_or_else(|| {
                Error::InvalidInput("XDATCAR ended before all coordinates".to_owned())
            })?;
            index += 1;
            let values: Vec<f64> = line
                .split_whitespace()
                .take(3)
                .map(|value| value.parse::<f64>())
                .collect::<Result<_, _>>()
                .map_err(|_| Error::InvalidInput("XDATCAR has a bad coordinate".to_owned()))?;
            if values.len() < 3 {
                return Err(Error::InvalidInput(
                    "XDATCAR coordinate needs three values".to_owned(),
                ));
            }
            let point = [values[0], values[1], values[2]];
            frac.push(if cartesian {
                lattice.to_frac(point)
            } else {
                point
            });
        }
        frames.push(Frame {
            frac,
            lattice: None,
            energy: None,
            step: None,
            info: BTreeMap::new(),
            forces: None,
            stress: None,
        });
    }
    Ok(Trajectory {
        elements,
        lattice: Some(matrix),
        frames,
    })
}

/// Detects LAMMPS dump, VASP `OUTCAR`/`vasprun.xml`, extended XYZ, DFTB+
/// `geo_end` or XDATCAR.
///
/// # Errors
///
/// Returns the error of the selected reader.
pub fn read_trajectory(text: &str) -> Result<Trajectory, Error> {
    if crate::vasp::looks_like_vasprun(text) {
        return crate::vasp::read_vasprun(text);
    }
    if crate::vasp::looks_like_outcar(text) {
        return crate::vasp::read_outcar(text);
    }
    let first = text
        .lines()
        .find(|line| !line.trim().is_empty())
        .unwrap_or("");
    if first.trim_start().starts_with("ITEM: TIMESTEP") {
        return read_lammps_dump(text);
    }
    let tokens: Vec<&str> = first.split_whitespace().collect();
    if tokens
        .first()
        .is_some_and(|token| token.parse::<usize>().is_ok())
    {
        return if tokens.len() == 1 {
            read_trajectory_xyz(text)
        } else {
            read_geo_end(text)
        };
    }
    read_xdatcar(text)
}
