//! Periodic neighbor search (Task 1.10).
//!
//! The index is a linked-cell grid in fractional coordinates. Its bin count
//! per axis comes from the perpendicular width of the cell, so any pair
//! within the cutoff lies in the same or an adjacent bin, also in triclinic
//! cells. Candidate pairs then enumerate the integer images that can reach the
//! cutoff and are filtered by the exact Cartesian distance.

use std::collections::HashMap;

use crate::spatial::{bin_counts, bin_index, offsets, perpendicular_widths, wrap};
use crate::{Error, Structure};

/// A pair of atoms within a cutoff, with the lattice image of `j`.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Neighbor {
    /// Index of the first atom in the P1 cell.
    pub i: usize,
    /// Index of the second atom in the P1 cell.
    pub j: usize,
    /// Lattice image of `j` relative to `i`.
    pub image: [i32; 3],
    /// Cartesian distance in Angstrom.
    pub distance: f64,
}

/// All periodic neighbors within `cutoff` Angstrom.
///
/// The search uses a linked-cell index in fractional space and is exact for
/// triclinic cells: bins are sized from the perpendicular widths of the cell,
/// and every candidate pair is checked against the cutoff with its Cartesian
/// distance. Results are sorted deterministically by `(i, j, image)`.
///
/// # Errors
///
/// Returns [`Error::InvalidInput`] for a non-positive or non-finite cutoff.
pub fn neighbors(structure: &Structure, cutoff: f64) -> Result<Vec<Neighbor>, Error> {
    validate_cutoff(cutoff)?;
    let p1 = structure.p1()?;
    let fracs: Vec<[f64; 3]> = p1.atoms.iter().map(|atom| atom.frac).collect();
    search(&p1.lattice, &fracs, cutoff)
}

fn validate_cutoff(cutoff: f64) -> Result<(), Error> {
    if !cutoff.is_finite() || cutoff <= 0.0 {
        return Err(Error::InvalidInput(format!(
            "neighbor cutoff must be positive and finite, got {cutoff}"
        )));
    }
    Ok(())
}

/// Neighbor search over fractional positions in one cell.
///
/// [`neighbors`] and the structure-of-arrays view share this path, so both
/// return the same sorted list.
pub(crate) fn search(
    lattice: &crate::Lattice,
    fracs: &[[f64; 3]],
    cutoff: f64,
) -> Result<Vec<Neighbor>, Error> {
    search_view(lattice, FractionalRows(fracs), cutoff)
}

/// Neighbor search directly over SoA coordinate columns (Task 17.2.3).
///
/// No `[f64; 3]` row array is materialised; the linked-cell grid reads the
/// three columns in place. The result is identical to [`search`].
pub(crate) fn search_columns(
    lattice: &crate::Lattice,
    x: &[f64],
    y: &[f64],
    z: &[f64],
    cutoff: f64,
) -> Result<Vec<Neighbor>, Error> {
    if x.len() != y.len() || x.len() != z.len() {
        return Err(Error::InvalidInput(
            "neighbor columns must have the same length".to_owned(),
        ));
    }
    search_view(lattice, Columns { x, y, z }, cutoff)
}

trait Coordinates {
    fn len(&self) -> usize;
    fn get(&self, index: usize) -> [f64; 3];
}

struct FractionalRows<'a>(&'a [[f64; 3]]);

impl Coordinates for FractionalRows<'_> {
    fn len(&self) -> usize {
        self.0.len()
    }

    fn get(&self, index: usize) -> [f64; 3] {
        self.0.get(index).copied().unwrap_or([0.0; 3])
    }
}

struct Columns<'a> {
    x: &'a [f64],
    y: &'a [f64],
    z: &'a [f64],
}

impl Coordinates for Columns<'_> {
    fn len(&self) -> usize {
        self.x.len()
    }

    fn get(&self, index: usize) -> [f64; 3] {
        [
            self.x.get(index).copied().unwrap_or(0.0),
            self.y.get(index).copied().unwrap_or(0.0),
            self.z.get(index).copied().unwrap_or(0.0),
        ]
    }
}

fn search_view<C: Coordinates>(
    lattice: &crate::Lattice,
    coordinates: C,
    cutoff: f64,
) -> Result<Vec<Neighbor>, Error> {
    validate_cutoff(cutoff)?;
    if coordinates.len() == 0 {
        return Ok(Vec::new());
    }
    let spacing = perpendicular_widths(lattice)?;
    let counts = bin_counts(spacing, cutoff);
    let matrix = lattice.matrix();

    let mut bins: HashMap<[usize; 3], Vec<usize>> = HashMap::new();
    let mut atom_bins: Vec<[usize; 3]> = Vec::with_capacity(coordinates.len());
    let [count_x, count_y, count_z] = counts;
    for index in 0..coordinates.len() {
        let [x, y, z] = coordinates.get(index);
        let frac = [x.rem_euclid(1.0), y.rem_euclid(1.0), z.rem_euclid(1.0)];
        let [frac_x, frac_y, frac_z] = frac;
        let bin = [
            bin_index(frac_x, count_x),
            bin_index(frac_y, count_y),
            bin_index(frac_z, count_z),
        ];
        bins.entry(bin).or_default().push(index);
        atom_bins.push(bin);
    }

    let search = Search {
        cutoff,
        spacing,
        matrix,
        coordinates: &coordinates,
    };
    let mut found = Vec::new();
    let offsets_x = offsets(count_x);
    let offsets_y = offsets(count_y);
    let offsets_z = offsets(count_z);
    for i in 0..coordinates.len() {
        let Some(base) = atom_bins.get(i) else {
            continue;
        };
        let [base_x, base_y, base_z] = *base;
        for &offset_x in offsets_x {
            let bin_x = wrap(base_x, offset_x, count_x);
            for &offset_y in offsets_y {
                let bin_y = wrap(base_y, offset_y, count_y);
                for &offset_z in offsets_z {
                    let bin_z = wrap(base_z, offset_z, count_z);
                    let Some(bucket) = bins.get(&[bin_x, bin_y, bin_z]) else {
                        continue;
                    };
                    for &j in bucket {
                        if j < i {
                            continue;
                        }
                        search.pair(i, j, &mut found);
                    }
                }
            }
        }
    }
    found.sort_by(|first, second| {
        (first.i, first.j, first.image, first.distance.to_bits()).cmp(&(
            second.i,
            second.j,
            second.image,
            second.distance.to_bits(),
        ))
    });
    Ok(found)
}

/// Neighbors of one atom within `cutoff`.
///
/// # Errors
///
/// Returns the errors of [`neighbors`].
pub fn neighbors_of(
    structure: &Structure,
    index: usize,
    cutoff: f64,
) -> Result<Vec<Neighbor>, Error> {
    Ok(neighbors(structure, cutoff)?
        .into_iter()
        .filter(|neighbor| neighbor.i == index || neighbor.j == index)
        .collect())
}

struct Search<'a, C: Coordinates> {
    cutoff: f64,
    spacing: [f64; 3],
    matrix: [[f64; 3]; 3],
    coordinates: &'a C,
}

impl<C: Coordinates> Search<'_, C> {
    fn pair(&self, i: usize, j: usize, found: &mut Vec<Neighbor>) {
        let first = self.coordinates.get(i);
        let second = self.coordinates.get(j);
        let [x_i, y_i, z_i] = first;
        let [x_j, y_j, z_j] = second;
        let delta = [x_j - x_i, y_j - y_i, z_j - z_i];
        let reach = self.spacing.map(|width| self.cutoff / width + 1e-9);
        let [delta_x, delta_y, delta_z] = delta;
        let [reach_x, reach_y, reach_z] = reach;
        let low_x = (-delta_x - reach_x).ceil() as i64;
        let high_x = (-delta_x + reach_x).floor() as i64;
        let low_y = (-delta_y - reach_y).ceil() as i64;
        let high_y = (-delta_y + reach_y).floor() as i64;
        let low_z = (-delta_z - reach_z).ceil() as i64;
        let high_z = (-delta_z + reach_z).floor() as i64;
        for image_x in low_x..=high_x {
            for image_y in low_y..=high_y {
                for image_z in low_z..=high_z {
                    if i == j && image_x == 0 && image_y == 0 && image_z == 0 {
                        continue;
                    }
                    let shift = [
                        delta_x + image_x as f64,
                        delta_y + image_y as f64,
                        delta_z + image_z as f64,
                    ];
                    let distance = norm(multiply_row(shift, self.matrix));
                    if distance > self.cutoff + 1e-12 {
                        continue;
                    }
                    if let (Ok(ix), Ok(iy), Ok(iz)) = (
                        i32::try_from(image_x),
                        i32::try_from(image_y),
                        i32::try_from(image_z),
                    ) {
                        found.push(Neighbor {
                            i,
                            j,
                            image: [ix, iy, iz],
                            distance,
                        });
                    }
                }
            }
        }
    }
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

fn norm(vector: [f64; 3]) -> f64 {
    vector.iter().map(|value| value * value).sum::<f64>().sqrt()
}
