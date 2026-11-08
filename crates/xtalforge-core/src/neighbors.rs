//! Periodic neighbor search (Task 1.10).

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
/// The search is direct and periodic: every pair is checked over the integer
/// images that can reach the cutoff, with the range derived from the smallest
/// interplanar spacing. Results are sorted deterministically by `(i, j, image)`.
///
/// # Errors
///
/// Returns [`Error::InvalidInput`] for a non-positive or non-finite cutoff.
pub fn neighbors(structure: &Structure, cutoff: f64) -> Result<Vec<Neighbor>, Error> {
    if !cutoff.is_finite() || cutoff <= 0.0 {
        return Err(Error::InvalidInput(format!(
            "neighbor cutoff must be positive and finite, got {cutoff}"
        )));
    }
    let p1 = structure.p1()?;
    if p1.atoms.is_empty() {
        return Ok(Vec::new());
    }
    let reciprocal = p1.lattice.reciprocal(false)?;
    let minimum_spacing = reciprocal
        .matrix()
        .iter()
        .map(|row| 1.0 / row.iter().map(|value| value * value).sum::<f64>().sqrt())
        .fold(f64::INFINITY, f64::min);
    let range = (cutoff / minimum_spacing).ceil() as i32 + 1;
    let matrix = p1.lattice.matrix();
    let mut found = Vec::new();
    for (i, atom_i) in p1.atoms.iter().enumerate() {
        for (j, atom_j) in p1.atoms.iter().enumerate().skip(i) {
            for x in -range..=range {
                for y in -range..=range {
                    for z in -range..=range {
                        let image = [x, y, z];
                        if i == j && image == [0, 0, 0] {
                            continue;
                        }
                        let frac_i = atom_i.frac;
                        let frac_j = atom_j.frac;
                        let delta = [
                            frac_j[0] + f64::from(x) - frac_i[0],
                            frac_j[1] + f64::from(y) - frac_i[1],
                            frac_j[2] + f64::from(z) - frac_i[2],
                        ];
                        let distance = norm(multiply_row(delta, matrix));
                        if distance <= cutoff + 1e-12 {
                            found.push(Neighbor {
                                i,
                                j,
                                image,
                                distance,
                            });
                        }
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
