#![allow(clippy::indexing_slicing, clippy::needless_range_loop)] // dense matrix kernel

//! Coulomb-matrix descriptor (Task 16.3.3).

use serde::{Deserialize, Serialize};
use xtalforge_core::{Error, Structure};

use super::{canonical_cell, to_json};

/// Options for [`coulomb_matrix`].
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct CoulombMatrixOptions {
    /// Largest accepted P1 atom count.
    ///
    /// The matrix has `n * n` entries, so the default caps it at 512 atoms
    /// (about 2 MiB of `f64`). The function returns an error instead of
    /// allocating an unbounded matrix.
    pub max_atoms: usize,
}

impl Default for CoulombMatrixOptions {
    fn default() -> Self {
        Self { max_atoms: 512 }
    }
}

/// The Coulomb matrix of a structure.
///
/// The matrix follows the standard molecular definition
/// (`M_ij = 0.5 Z_i^2.4` on the diagonal and `M_ij = Z_i Z_j / r_ij`
/// elsewhere), with the minimum-image distance for the periodic structure.
/// Atoms are in canonical order: atomic number, fractional position, label.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct CoulombMatrix {
    /// Element symbols in canonical atom order.
    pub elements: Vec<String>,
    /// Atomic numbers in canonical atom order.
    pub atomic_numbers: Vec<u8>,
    /// Fractional positions in canonical atom order.
    pub positions: Vec<[f64; 3]>,
    /// Symmetric matrix, row-major.
    pub matrix: Vec<Vec<f64>>,
}

impl CoulombMatrix {
    /// Number of atoms and the side of the square matrix.
    #[must_use]
    pub fn n_atoms(&self) -> usize {
        self.atomic_numbers.len()
    }

    /// Feature vector with the matrix flattened row by row.
    #[must_use]
    pub fn to_vector(&self) -> Vec<f64> {
        self.matrix
            .iter()
            .flat_map(|row| row.iter().copied())
            .collect()
    }

    /// Serializes the matrix as pretty JSON.
    ///
    /// # Errors
    ///
    /// Returns [`Error::InvalidInput`] when the serializer fails.
    pub fn to_json(&self) -> Result<String, Error> {
        to_json(self, "Coulomb matrix")
    }
}

/// Computes the Coulomb matrix of a structure.
///
/// # Errors
///
/// Returns [`Error::InvalidInput`] for an empty structure, an unknown element,
/// a zero maximum atom count, a structure above the maximum, or two sites at
/// the same position (the matrix would be infinite).
pub fn coulomb_matrix(
    structure: &Structure,
    options: &CoulombMatrixOptions,
) -> Result<CoulombMatrix, Error> {
    if options.max_atoms == 0 {
        return Err(Error::InvalidInput(
            "Coulomb matrix max_atoms must be positive".to_owned(),
        ));
    }
    let (lattice, sites) = canonical_cell(structure)?;
    if sites.len() > options.max_atoms {
        return Err(Error::InvalidInput(format!(
            "Coulomb matrix needs {} atoms but max_atoms is {}",
            sites.len(),
            options.max_atoms
        )));
    }
    let n_atoms = sites.len();
    let mut matrix = vec![vec![0.0; n_atoms]; n_atoms];
    for first in 0..n_atoms {
        let z_first = f64::from(sites[first].z);
        matrix[first][first] = 0.5 * z_first.powf(2.4);
        for second in (first + 1)..n_atoms {
            let z_second = f64::from(sites[second].z);
            let delta = [
                sites[second].frac[0] - sites[first].frac[0],
                sites[second].frac[1] - sites[first].frac[1],
                sites[second].frac[2] - sites[first].frac[2],
            ];
            let (_, distance) = lattice.minimum_image(delta);
            if !distance.is_finite() || distance <= 1e-8 {
                return Err(Error::InvalidInput(format!(
                    "Coulomb matrix needs distinct sites; atoms {first} and {second} coincide"
                )));
            }
            let value = z_first * z_second / distance;
            matrix[first][second] = value;
            matrix[second][first] = value;
        }
    }
    Ok(CoulombMatrix {
        elements: sites.iter().map(|site| site.element.clone()).collect(),
        atomic_numbers: sites.iter().map(|site| site.z).collect(),
        positions: sites.iter().map(|site| site.frac).collect(),
        matrix,
    })
}
