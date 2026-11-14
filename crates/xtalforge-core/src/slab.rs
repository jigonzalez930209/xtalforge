//! Slabs cut along a lattice plane (Task 1.8).
//!
//! The surface basis is ported from Crystal-Builder `core/slab.py`. The bonds
//! of the parent cell are carried down to the copies they join: every bond is
//! taken to each copy of its first atom and kept when its partner is in the
//! slab too. A bond whose partner is outside the slab is severed by the cut and
//! counted in [`Slab::cut`]. Bonds are never perceived afresh from the slab
//! geometry, so a bond the user removed stays removed and a bond they drew
//! keeps its order.

#![allow(clippy::indexing_slicing)] // fixed 3x3 integer matrices

use std::f64::consts::PI;

use crate::bonding::{BondRules, perceive_bonds};
use crate::{Error, Lattice, Site, Structure, StructureBuilder};

mod carry;
mod terminations;

pub(crate) use terminations::validate as validate_termination_options;
pub use terminations::{
    SlabTermination, SlabTerminationReport, TerminationOptions, slab_cut_shifts, slab_terminations,
    slab_terminations_with_options,
};

/// How far above the bottom of the cell a slab starts, in Angstrom.
pub const SLAB_FLOOR: f64 = 0.5;

/// A slab and what it took to make it.
#[derive(Debug, Clone, PartialEq)]
pub struct Slab {
    /// The slab structure, always P1.
    pub structure: Structure,
    /// Reduced Miller index.
    pub hkl: [i32; 3],
    /// Number of layers.
    pub layers: u32,
    /// Slab thickness in Angstrom.
    pub thickness: f64,
    /// Vacuum above the slab in Angstrom.
    pub vacuum: f64,
    /// Chemical bonds severed by the cut.
    pub cut: usize,
}

/// Integer basis with v1 and v2 in the plane `(hkl)` and `h . v3 = 1`.
///
/// # Errors
///
/// Returns [`Error::InvalidInput`] for `(0, 0, 0)`.
pub fn surface_basis(hkl: [i32; 3]) -> Result<[[i32; 3]; 3], Error> {
    let [h0, h1, h2] = hkl;
    if h0 == 0 && h1 == 0 && h2 == 0 {
        return Err(Error::InvalidInput(
            "a plane needs a Miller index that is not (0 0 0)".to_owned(),
        ));
    }
    let g = gcd(gcd(h0.abs(), h1.abs()), h2.abs()).max(1);
    let mut r = [h0 / g, h1 / g, h2 / g];
    let mut q = [[1i32, 0, 0], [0, 1, 0], [0, 0, 1]];
    while r.iter().filter(|value| **value != 0).count() > 1 {
        let pivot = (0..3)
            .filter(|index| r[*index] != 0)
            .min_by_key(|index| r[*index].abs())
            .unwrap_or(0);
        for k in 0..3 {
            if k != pivot && r[k] != 0 {
                let m = r[k].div_euclid(r[pivot]);
                r[k] -= m * r[pivot];
                for row in &mut q {
                    row[k] -= m * row[pivot];
                }
            }
        }
    }
    let last = (0..3).find(|index| r[*index] != 0).unwrap_or(0);
    if r[last] < 0 {
        for row in &mut q {
            row[last] = -row[last];
        }
    }
    let order = {
        let mut order = Vec::new();
        for index in 0..3 {
            if index != last {
                order.push(index);
            }
        }
        order.push(last);
        order
    };
    let mut basis = [[0i32; 3]; 3];
    for (row, column) in order.iter().enumerate() {
        for index in 0..3 {
            basis[row][index] = q[index][*column];
        }
    }
    if determinant_i32(basis) < 0 {
        basis[0] = [-basis[0][0], -basis[0][1], -basis[0][2]];
    }
    Ok(basis)
}

/// Builds `layers` spacings of `(hkl)` with `vacuum` Angstrom above.
///
/// # Errors
///
/// Returns [`Error::InvalidInput`] for an empty structure, a zero or negative
/// layer count, a negative vacuum or a degenerate plane.
pub fn make_slab(
    structure: &Structure,
    hkl: [i32; 3],
    layers: u32,
    vacuum: f64,
    shift: f64,
) -> Result<Slab, Error> {
    if layers == 0 {
        return Err(Error::InvalidInput(
            "a slab needs at least one layer".to_owned(),
        ));
    }
    if !vacuum.is_finite() || vacuum < 0.0 {
        return Err(Error::InvalidInput(
            "the vacuum cannot be negative".to_owned(),
        ));
    }
    if !shift.is_finite() {
        return Err(Error::InvalidInput(
            "the slab shift must be finite".to_owned(),
        ));
    }
    let p1 = structure.p1()?;
    if p1.atoms.is_empty() {
        return Err(Error::InvalidInput(
            "there are no atoms to cut a slab from".to_owned(),
        ));
    }
    let graph = perceive_bonds(structure, BondRules::from_structure(structure)?)?;
    let basis = surface_basis(hkl)?;
    let basis_float = basis.map(|row| row.map(f64::from));
    let inverse = Lattice::new(basis_float)?.inverse();
    let matrix = p1.lattice.matrix();
    let shift_vector = basis_float[2].map(|value| value * shift);

    let mut y = Vec::with_capacity(p1.atoms.len());
    for atom in &p1.atoms {
        let shifted = [
            atom.frac[0] - shift_vector[0],
            atom.frac[1] - shift_vector[1],
            atom.frac[2] - shift_vector[2],
        ];
        y.push(multiply_row(shifted, inverse));
    }
    let whole: Vec<[i32; 3]> = y
        .iter()
        .map(|point| point.map(|value| (value + 1e-9).floor() as i32))
        .collect();
    let first: Vec<[f64; 3]> = y
        .iter()
        .zip(&whole)
        .map(|(point, whole)| {
            [
                point[0] - f64::from(whole[0]),
                point[1] - f64::from(whole[1]),
                point[2] - f64::from(whole[2]),
            ]
        })
        .collect();

    let a = multiply_row(basis_float[0], matrix);
    let b = multiply_row(basis_float[1], matrix);
    let c = multiply_row(basis_float[2], matrix);
    let normal = normalize(cross(a, b))?;
    let rise = dot(c, normal);
    if rise.abs() < 1e-12 {
        return Err(Error::InvalidInput(
            "the plane is degenerate in this cell".to_owned(),
        ));
    }
    let in_plane = [
        c[0] - rise * normal[0],
        c[1] - rise * normal[1],
        c[2] - rise * normal[2],
    ];
    let along = solve_2x2(a, b, in_plane)?;
    let thickness = rise.abs() * f64::from(layers);
    let height = thickness + vacuum;
    let floor = SLAB_FLOOR.min(vacuum / 2.0) / height;
    let sign = if rise < 0.0 { -1.0 } else { 1.0 };
    let lattice = Lattice::new([
        a,
        b,
        [
            normal[0] * sign * height,
            normal[1] * sign * height,
            normal[2] * sign * height,
        ],
    ])?;

    let mut builder = StructureBuilder::new(structure.name(), lattice);
    let mut cells = Vec::with_capacity(p1.atoms.len() * layers as usize);
    for layer in 0..layers {
        for (index, atom) in p1.atoms.iter().enumerate() {
            let z = first[index][2] + f64::from(layer);
            let plane = [
                first[index][0] + along[0] * z,
                first[index][1] + along[1] * z,
            ];
            let cell = [(plane[0] + 1e-9).floor(), (plane[1] + 1e-9).floor()];
            let frac = [
                plane[0] - cell[0],
                plane[1] - cell[1],
                floor + z * rise.abs() / height,
            ];
            let site = Site::new(
                format!("{}{}", atom.element, index + 1),
                &atom.element,
                frac,
            )?;
            builder = builder.site(site.with_occupancy(atom.occupancy)?);
            cells.push([cell[0] as i32, cell[1] as i32]);
        }
    }
    let n_atoms = p1.atoms.len();
    let (carried, cut) = carry::SlabCarry {
        structure,
        p1: &p1,
        graph: &graph,
        layers,
        whole: &whole,
        cells: &cells,
        inverse,
        n_atoms,
    }
    .bonds();
    for bond in carried {
        builder = builder.explicit_bond(bond);
    }
    for (key, value) in structure.metadata() {
        builder = builder.metadata(key.clone(), value.clone());
    }
    let gcd_value = gcd(gcd(hkl[0].abs(), hkl[1].abs()), hkl[2].abs()).max(1);
    Ok(Slab {
        structure: builder.build()?,
        hkl: [hkl[0] / gcd_value, hkl[1] / gcd_value, hkl[2] / gcd_value],
        layers,
        thickness,
        vacuum,
        cut,
    })
}

/// Angle between two planes in degrees, for a quick surface check.
#[must_use]
pub fn plane_angle(first: [i32; 3], second: [i32; 3]) -> f64 {
    let dot = f64::from(first[0] * second[0] + first[1] * second[1] + first[2] * second[2]);
    let first_norm =
        f64::from(first[0] * first[0] + first[1] * first[1] + first[2] * first[2]).sqrt();
    let second_norm =
        f64::from(second[0] * second[0] + second[1] * second[1] + second[2] * second[2]).sqrt();
    (dot / (first_norm * second_norm)).clamp(-1.0, 1.0).acos() * 180.0 / PI
}

fn gcd(first: i32, second: i32) -> i32 {
    let mut a = first;
    let mut b = second;
    while b != 0 {
        let remainder = a % b;
        a = b;
        b = remainder;
    }
    a.abs()
}

fn determinant_i32(matrix: [[i32; 3]; 3]) -> i32 {
    let [[a, b, c], [d, e, f], [g, h, i]] = matrix;
    a * (e * i - f * h) - b * (d * i - f * g) + c * (d * h - e * g)
}

fn multiply_row(vector: [f64; 3], matrix: [[f64; 3]; 3]) -> [f64; 3] {
    let [x, y, z] = vector;
    let [[a, b, c], [d, e, f], [g, h, i]] = matrix;
    [
        x * a + y * d + z * g,
        x * b + y * e + z * h,
        x * c + y * f + z * i,
    ]
}

fn cross(first: [f64; 3], second: [f64; 3]) -> [f64; 3] {
    let [a, b, c] = first;
    let [d, e, f] = second;
    [b * f - c * e, c * d - a * f, a * e - b * d]
}

fn dot(first: [f64; 3], second: [f64; 3]) -> f64 {
    first[0] * second[0] + first[1] * second[1] + first[2] * second[2]
}

fn normalize(vector: [f64; 3]) -> Result<[f64; 3], Error> {
    let length = dot(vector, vector).sqrt();
    if length < 1e-12 {
        return Err(Error::InvalidInput("the plane has no normal".to_owned()));
    }
    Ok([vector[0] / length, vector[1] / length, vector[2] / length])
}

fn solve_2x2(a: [f64; 3], b: [f64; 3], target: [f64; 3]) -> Result<[f64; 3], Error> {
    let aa = dot(a, a);
    let ab = dot(a, b);
    let bb = dot(b, b);
    let ta = dot(target, a);
    let tb = dot(target, b);
    let det = aa * bb - ab * ab;
    if det.abs() < 1e-12 {
        return Err(Error::InvalidInput(
            "the surface vectors are degenerate".to_owned(),
        ));
    }
    let x = (ta * bb - tb * ab) / det;
    let y = (tb * aa - ta * ab) / det;
    Ok([x, y, 0.0])
}
