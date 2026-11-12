//! The supercell closest to a cube for a target cell count.

use crate::{Error, Structure};

/// The supercell that is closest to a cube for a target cell count.
#[derive(Debug, Clone, PartialEq)]
pub struct SupercellChoice {
    /// Integer transformation matrix; rows are new vectors in old units.
    pub matrix: [[f64; 3]; 3],
    /// The transformed structure.
    pub structure: Structure,
    /// Determinant of the matrix (number of original cells).
    pub cells: u32,
    /// Shape score; lower is closer to a cube.
    pub score: f64,
}

/// Finds the integer supercell closest to a cube for `cells` original cells.
///
/// The search enumerates short integer vectors and keeps the triple with the
/// determinant `cells` that minimizes length variance and angle deviation
/// (ASE's `find_optimal_cell_shape` objective).
///
/// # Errors
///
/// Returns [`Error::InvalidInput`] for `cells == 0` or when no triple exists,
/// or a supercell error.
pub fn optimal_supercell(structure: &Structure, cells: u32) -> Result<SupercellChoice, Error> {
    if cells == 0 {
        return Err(Error::InvalidInput(
            "the target cell count must be positive".to_owned(),
        ));
    }
    let limit = (f64::from(cells).cbrt() * 2.0 + 3.0).ceil() as i32;
    let mut vectors: Vec<[i32; 3]> = Vec::new();
    for x in -limit..=limit {
        for y in -limit..=limit {
            for z in -limit..=limit {
                if x == 0 && y == 0 && z == 0 {
                    continue;
                }
                vectors.push([x, y, z]);
            }
        }
    }
    vectors.sort_by_key(|vector| {
        (
            vector[0] * vector[0] + vector[1] * vector[1] + vector[2] * vector[2],
            *vector,
        )
    });
    vectors.truncate(140);
    let target = i64::from(cells);
    let matrix = structure.lattice().matrix();
    let mut best: Option<([[f64; 3]; 3], f64)> = None;
    for first in 0..vectors.len() {
        for second in (first + 1)..vectors.len() {
            for third in (second + 1)..vectors.len() {
                let (Some(&a), Some(&b), Some(&c)) =
                    (vectors.get(first), vectors.get(second), vectors.get(third))
                else {
                    continue;
                };
                let determinant = determinant_i32([a, b, c]);
                if determinant != target {
                    continue;
                }
                let rows = [a, b, c].map(|row| row.map(f64::from));
                let transformed = multiply(rows, matrix);
                let lengths = transformed.map(norm);
                let mean = (lengths[0] + lengths[1] + lengths[2]) / 3.0;
                if mean <= 0.0 {
                    continue;
                }
                let mut score = 0.0;
                for length in lengths {
                    score += ((length - mean) / mean).powi(2);
                }
                for (i, j) in [(0, 1), (0, 2), (1, 2)] {
                    let cosine = dot(transformed[i], transformed[j]) / (lengths[i] * lengths[j]);
                    score += cosine * cosine;
                }
                let candidate = (rows, score);
                let replace = match &best {
                    Some((current, current_score)) => {
                        if score < *current_score - 1e-12 {
                            true
                        } else if (score - *current_score).abs() <= 1e-12 {
                            rows < *current
                        } else {
                            false
                        }
                    }
                    None => true,
                };
                if replace {
                    best = Some(candidate);
                }
            }
        }
    }
    let Some((matrix, score)) = best else {
        return Err(Error::InvalidInput(format!(
            "no supercell with {cells} cells was found"
        )));
    };
    let structure = structure.supercell_matrix(matrix)?;
    Ok(SupercellChoice {
        matrix,
        structure,
        cells,
        score,
    })
}

fn determinant_i32(matrix: [[i32; 3]; 3]) -> i64 {
    let [[a, b, c], [d, e, f], [g, h, i]] = matrix;
    i64::from(a) * (i64::from(e) * i64::from(i) - i64::from(f) * i64::from(h))
        - i64::from(b) * (i64::from(d) * i64::from(i) - i64::from(f) * i64::from(g))
        + i64::from(c) * (i64::from(d) * i64::from(h) - i64::from(e) * i64::from(g))
}

fn multiply(left: [[f64; 3]; 3], right: [[f64; 3]; 3]) -> [[f64; 3]; 3] {
    let mut result = [[0.0; 3]; 3];
    for row in 0..3 {
        for column in 0..3 {
            for k in 0..3 {
                result[row][column] += left[row][k] * right[k][column];
            }
        }
    }
    result
}

fn dot(first: [f64; 3], second: [f64; 3]) -> f64 {
    first[0] * second[0] + first[1] * second[1] + first[2] * second[2]
}

fn norm(vector: [f64; 3]) -> f64 {
    dot(vector, vector).sqrt()
}
