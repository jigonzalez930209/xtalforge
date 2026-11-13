//! LAMMPS dump and DFTB+ `geo_end` trajectories (Task 3.6.3).

use std::collections::BTreeMap;

use xtalforge_core::{Error, Lattice};

use super::{Frame, Trajectory, read_trajectory_xyz};
use crate::read_gen;

const PAD: f64 = 5.0;

/// Reads a LAMMPS text dump with an `element` column.
///
/// # Errors
///
/// Returns [`Error::InvalidInput`] for a missing element column, a truncated
/// frame, a bad box or an inconsistent element sequence. Use
/// [`read_lammps_dump_with_types`] when the dump only has numeric types.
pub fn read_lammps_dump(text: &str) -> Result<Trajectory, Error> {
    read_lammps_dump_with_types(text, &[])
}

/// Reads a LAMMPS text dump, mapping numeric types through `types`.
///
/// `types` is indexed by LAMMPS type number minus one. An `element` column
/// takes precedence when present.
///
/// # Errors
///
/// Returns [`Error::InvalidInput`] for a truncated frame, a bad box, an
/// unknown type or an inconsistent element sequence.
pub fn read_lammps_dump_with_types(text: &str, types: &[&str]) -> Result<Trajectory, Error> {
    let lines: Vec<&str> = text.lines().collect();
    let mut index = 0;
    let mut elements: Vec<String> = Vec::new();
    let mut frames = Vec::new();
    let mut matrix = None;
    while index < lines.len() {
        if lines[index].trim() != "ITEM: TIMESTEP" {
            index += 1;
            continue;
        }
        index += 1;
        let step = parse_u64(lines.get(index), "LAMMPS timestep")?;
        index += 1;
        expect(lines.get(index), "ITEM: NUMBER OF ATOMS")?;
        index += 1;
        let count = parse_usize(lines.get(index), "LAMMPS atom count")?;
        index += 1;
        let header = lines.get(index).ok_or_else(|| truncated("box header"))?;
        if !header.starts_with("ITEM: BOX BOUNDS") {
            return Err(Error::InvalidInput(
                "LAMMPS dump is missing ITEM: BOX BOUNDS".to_owned(),
            ));
        }
        index += 1;
        let mut bounds = [[0.0; 3]; 3];
        for (row, bound) in bounds.iter_mut().enumerate() {
            let line = lines.get(index).ok_or_else(|| truncated("box bounds"))?;
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
            *bound = [values[0], values[1], values.get(2).copied().unwrap_or(0.0)];
            index += 1;
            let _ = row;
        }
        let cell = lammps_lattice(bounds)?;
        if matrix.is_none() {
            matrix = Some(cell.matrix());
        }
        let atom_header = lines.get(index).ok_or_else(|| truncated("atoms header"))?;
        if !atom_header.starts_with("ITEM: ATOMS") {
            return Err(Error::InvalidInput(
                "LAMMPS dump is missing ITEM: ATOMS".to_owned(),
            ));
        }
        let columns: Vec<&str> = atom_header.split_whitespace().skip(2).collect();
        index += 1;
        let element_column = columns.iter().position(|column| *column == "element");
        let type_column = columns.iter().position(|column| *column == "type");
        let (mode, x, y, z) = coordinate_columns(&columns)?;
        let mut frame_elements = Vec::with_capacity(count);
        let mut frac = Vec::with_capacity(count);
        for _ in 0..count {
            let line = lines.get(index).ok_or_else(|| truncated("atom row"))?;
            index += 1;
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
                types
                    .get(kind.saturating_sub(1))
                    .map(|value| (*value).to_owned())
                    .ok_or_else(|| {
                        Error::InvalidInput(format!("LAMMPS type {kind} has no element mapping"))
                    })?
            } else {
                return Err(Error::InvalidInput(
                    "LAMMPS dump needs an element column or a type mapping".to_owned(),
                ));
            };
            let point = [value(&values, x)?, value(&values, y)?, value(&values, z)?];
            frame_elements.push(element);
            frac.push(match mode {
                CoordinateMode::Fractional => point,
                CoordinateMode::Cartesian => cell.to_frac(point),
            });
        }
        if elements.is_empty() {
            elements = frame_elements;
        } else if elements != frame_elements {
            return Err(Error::InvalidInput(
                "LAMMPS frames disagree on atom order or elements".to_owned(),
            ));
        }
        frames.push(Frame {
            frac,
            lattice: None,
            energy: None,
            step: Some(step),
            info: BTreeMap::new(),
            forces: None,
            stress: None,
        });
    }
    if frames.is_empty() {
        return Err(Error::InvalidInput("LAMMPS dump has no frames".to_owned()));
    }
    Ok(Trajectory {
        elements,
        lattice: matrix,
        frames,
    })
}

/// Reads a DFTB+ `geo_end.xyz` or `geo_end.gen`.
///
/// A plain XYZ without a cell is placed in one shared cubic box with 5 A of
/// vacuum, the same convention as a molecular structure. A `.gen` final
/// geometry becomes a one-frame trajectory.
///
/// # Errors
///
/// Returns the errors of [`read_trajectory_xyz`] or [`read_gen`].
pub fn read_geo_end(text: &str) -> Result<Trajectory, Error> {
    let first = text
        .lines()
        .find(|line| !line.trim().is_empty())
        .unwrap_or("");
    let mut tokens = first.split_whitespace();
    let count = tokens.next().and_then(|value| value.parse::<usize>().ok());
    if count.is_some() && tokens.next().is_some() {
        let structure = read_gen(text)?;
        let p1 = structure.p1()?;
        let frame = Frame {
            frac: p1.atoms.iter().map(|atom| atom.frac).collect(),
            lattice: None,
            energy: None,
            step: None,
            info: BTreeMap::new(),
            forces: None,
            stress: None,
        };
        return Ok(Trajectory {
            elements: p1.atoms.iter().map(|atom| atom.element.clone()).collect(),
            lattice: Some(structure.lattice().matrix()),
            frames: vec![frame],
        });
    }
    Ok(boxed(read_trajectory_xyz(text)?))
}

fn boxed(mut trajectory: Trajectory) -> Trajectory {
    if trajectory.lattice.is_some() || trajectory.frames.is_empty() {
        return trajectory;
    }
    let mut low = [f64::INFINITY; 3];
    let mut high = [f64::NEG_INFINITY; 3];
    for frame in &trajectory.frames {
        for point in &frame.frac {
            for axis in 0..3 {
                low[axis] = low[axis].min(point[axis]);
                high[axis] = high[axis].max(point[axis]);
            }
        }
    }
    if !low.iter().chain(high.iter()).all(|value| value.is_finite()) {
        return trajectory;
    }
    let side = (0..3)
        .map(|axis| high[axis] - low[axis] + 2.0 * PAD)
        .fold(10.0, f64::max);
    let Ok(lattice) = Lattice::cubic(side) else {
        return trajectory;
    };
    let shift = [PAD - low[0], PAD - low[1], PAD - low[2]];
    for frame in &mut trajectory.frames {
        for point in &mut frame.frac {
            for axis in 0..3 {
                point[axis] += shift[axis];
            }
        }
    }
    trajectory.lattice = Some(lattice.matrix());
    trajectory
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

fn expect(line: Option<&&str>, expected: &str) -> Result<(), Error> {
    if line.is_some_and(|line| line.trim() == expected) {
        Ok(())
    } else {
        Err(Error::InvalidInput(format!(
            "LAMMPS dump is missing {expected}"
        )))
    }
}

fn parse_u64(line: Option<&&str>, label: &str) -> Result<u64, Error> {
    line.ok_or_else(|| truncated(label))?
        .trim()
        .parse()
        .map_err(|_| Error::InvalidInput(format!("{label} is not an integer")))
}

fn parse_usize(line: Option<&&str>, label: &str) -> Result<usize, Error> {
    line.ok_or_else(|| truncated(label))?
        .trim()
        .parse()
        .map_err(|_| Error::InvalidInput(format!("{label} is not an integer")))
}

fn value(values: &[&str], index: usize) -> Result<f64, Error> {
    values
        .get(index)
        .ok_or_else(|| truncated("coordinate column"))?
        .parse()
        .map_err(|_| Error::InvalidInput("LAMMPS coordinate is not a number".to_owned()))
}

fn truncated(part: &str) -> Error {
    Error::InvalidInput(format!("LAMMPS dump ended before {part}"))
}
