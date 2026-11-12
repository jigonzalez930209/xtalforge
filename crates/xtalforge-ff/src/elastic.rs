#![allow(clippy::indexing_slicing, clippy::needless_range_loop)] // fixed-size numeric kernels

//! Elastic tensor from stress finite differences (Task 9.6.1).
//!
//! Six Voigt strain components are applied with a central difference, the
//! analytic stress of any [`Calculator`] is read, and the 6x6 stiffness is
//! inverted to Voigt-Reuss-Hill moduli. Stability is the positive
//! definiteness of the stiffness matrix.

use serde::{Deserialize, Serialize};

use xtalforge_core::{Error, Structure, voigt_to_tensor};

use crate::{Calculator, Error as FfError};

/// Options for [`elastic_tensor`].
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct ElasticOptions {
    /// Magnitude of the applied engineering strain.
    pub strain: f64,
}

impl Default for ElasticOptions {
    fn default() -> Self {
        Self { strain: 0.005 }
    }
}

/// Elastic tensor and derived moduli.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct ElasticTensor {
    /// Stiffness matrix in GPa, Voigt notation.
    pub stiffness: [[f64; 6]; 6],
    /// Compliance matrix in 1/GPa.
    pub compliance: [[f64; 6]; 6],
    /// Voigt bulk modulus in GPa.
    pub bulk_modulus_voigt: f64,
    /// Reuss bulk modulus in GPa.
    pub bulk_modulus_reuss: f64,
    /// Hill bulk modulus in GPa.
    pub bulk_modulus_hill: f64,
    /// Voigt shear modulus in GPa.
    pub shear_modulus_voigt: f64,
    /// Reuss shear modulus in GPa.
    pub shear_modulus_reuss: f64,
    /// Hill shear modulus in GPa.
    pub shear_modulus_hill: f64,
    /// Hill Young modulus in GPa.
    pub young_modulus: f64,
    /// Hill Poisson ratio.
    pub poisson_ratio: f64,
    /// Whether the stiffness matrix is positive definite.
    pub stable: bool,
    /// Number of stress evaluations.
    pub evaluations: usize,
}

/// Computes the elastic tensor by central differences of the stress.
///
/// # Errors
///
/// Returns [`Error::InvalidInput`] for a non-positive strain, a calculator
/// without stress, or a singular stiffness matrix, and the calculator errors.
pub fn elastic_tensor(
    structure: &Structure,
    calculator: &dyn Calculator,
    options: &ElasticOptions,
) -> Result<ElasticTensor, Error> {
    if !options.strain.is_finite() || options.strain <= 0.0 {
        return Err(Error::InvalidInput(
            "the elastic strain must be positive and finite".to_owned(),
        ));
    }
    if !calculator.provides_stress() {
        return Err(Error::InvalidInput(format!(
            "calculator {} does not provide stress",
            calculator.name()
        )));
    }
    let mut stiffness = [[0.0; 6]; 6];
    let mut evaluations = 0;
    for component in 0..6 {
        let mut plus = [0.0; 6];
        let mut minus = [0.0; 6];
        plus[component] = options.strain;
        minus[component] = -options.strain;
        let plus_structure = xtalforge_core::strain(structure, voigt_to_tensor(plus)?)?;
        let minus_structure = xtalforge_core::strain(structure, voigt_to_tensor(minus)?)?;
        let plus_stress = calculator
            .stress(&plus_structure)
            .map_err(map_error)?
            .ok_or_else(|| {
                Error::InvalidInput(format!(
                    "calculator {} returned no stress",
                    calculator.name()
                ))
            })?;
        let minus_stress = calculator
            .stress(&minus_structure)
            .map_err(map_error)?
            .ok_or_else(|| {
                Error::InvalidInput(format!(
                    "calculator {} returned no stress",
                    calculator.name()
                ))
            })?;
        evaluations += 2;
        let plus_voigt = stress_voigt(plus_stress);
        let minus_voigt = stress_voigt(minus_stress);
        for row in 0..6 {
            stiffness[row][component] =
                (plus_voigt[row] - minus_voigt[row]) / (2.0 * options.strain);
        }
    }
    for row in 0..6 {
        for column in (row + 1)..6 {
            let value = 0.5 * (stiffness[row][column] + stiffness[column][row]);
            stiffness[row][column] = value;
            stiffness[column][row] = value;
        }
    }
    let Some(compliance) = inverse6(stiffness) else {
        return Err(Error::InvalidInput(
            "the stiffness matrix is singular".to_owned(),
        ));
    };
    let bulk_voigt = (stiffness[0][0]
        + stiffness[1][1]
        + stiffness[2][2]
        + 2.0 * (stiffness[0][1] + stiffness[0][2] + stiffness[1][2]))
        / 9.0;
    let shear_voigt = (stiffness[0][0] + stiffness[1][1] + stiffness[2][2]
        - (stiffness[0][1] + stiffness[0][2] + stiffness[1][2])
        + 3.0 * (stiffness[3][3] + stiffness[4][4] + stiffness[5][5]))
        / 15.0;
    let bulk_reuss = 1.0
        / (compliance[0][0]
            + compliance[1][1]
            + compliance[2][2]
            + 2.0 * (compliance[0][1] + compliance[0][2] + compliance[1][2]));
    let shear_reuss = 15.0
        / (4.0 * (compliance[0][0] + compliance[1][1] + compliance[2][2])
            - 4.0 * (compliance[0][1] + compliance[0][2] + compliance[1][2])
            + 3.0 * (compliance[3][3] + compliance[4][4] + compliance[5][5]));
    let bulk_hill = 0.5 * (bulk_voigt + bulk_reuss);
    let shear_hill = 0.5 * (shear_voigt + shear_reuss);
    let young = if shear_hill.abs() > 1e-12 {
        9.0 * bulk_hill * shear_hill / (3.0 * bulk_hill + shear_hill)
    } else {
        0.0
    };
    let poisson = if bulk_hill.abs() > 1e-12 {
        (3.0 * bulk_hill - 2.0 * shear_hill) / (2.0 * (3.0 * bulk_hill + shear_hill))
    } else {
        0.0
    };
    Ok(ElasticTensor {
        stiffness,
        compliance,
        bulk_modulus_voigt: bulk_voigt,
        bulk_modulus_reuss: bulk_reuss,
        bulk_modulus_hill: bulk_hill,
        shear_modulus_voigt: shear_voigt,
        shear_modulus_reuss: shear_reuss,
        shear_modulus_hill: shear_hill,
        young_modulus: young,
        poisson_ratio: poisson,
        stable: cholesky(stiffness).is_some(),
        evaluations,
    })
}

/// Engineering strain tensors of the six Voigt deformations.
#[must_use]
pub fn deformation_set(strain: f64) -> Vec<[[f64; 3]; 3]> {
    let mut set = Vec::with_capacity(6);
    for component in 0..6 {
        let mut voigt = [0.0; 6];
        voigt[component] = strain;
        if let Ok(tensor) = voigt_to_tensor(voigt) {
            set.push(tensor);
        }
    }
    set
}

fn stress_voigt(stress: [[f64; 3]; 3]) -> [f64; 6] {
    [
        stress[0][0],
        stress[1][1],
        stress[2][2],
        stress[1][2],
        stress[0][2],
        stress[0][1],
    ]
}

fn map_error(error: FfError) -> Error {
    match error {
        FfError::InvalidInput(message) => Error::InvalidInput(message),
        other => Error::InvalidInput(other.to_string()),
    }
}

fn inverse6(matrix: [[f64; 6]; 6]) -> Option<[[f64; 6]; 6]> {
    let mut augmented = [[0.0; 12]; 6];
    for row in 0..6 {
        for column in 0..6 {
            augmented[row][column] = matrix[row][column];
        }
        augmented[row][6 + row] = 1.0;
    }
    for pivot in 0..6 {
        let mut best = pivot;
        for row in (pivot + 1)..6 {
            if augmented[row][pivot].abs() > augmented[best][pivot].abs() {
                best = row;
            }
        }
        if augmented[best][pivot].abs() < 1e-14 {
            return None;
        }
        augmented.swap(pivot, best);
        let divisor = augmented[pivot][pivot];
        for column in 0..12 {
            augmented[pivot][column] /= divisor;
        }
        for row in 0..6 {
            if row == pivot {
                continue;
            }
            let factor = augmented[row][pivot];
            for column in 0..12 {
                augmented[row][column] -= factor * augmented[pivot][column];
            }
        }
    }
    let mut inverse = [[0.0; 6]; 6];
    for row in 0..6 {
        for column in 0..6 {
            inverse[row][column] = augmented[row][6 + column];
        }
    }
    Some(inverse)
}

fn cholesky(matrix: [[f64; 6]; 6]) -> Option<()> {
    let mut lower = [[0.0; 6]; 6];
    for row in 0..6 {
        for column in 0..=row {
            let mut sum = matrix[row][column];
            for k in 0..column {
                sum -= lower[row][k] * lower[column][k];
            }
            if row == column {
                if sum <= 1e-12 {
                    return None;
                }
                lower[row][column] = sum.sqrt();
            } else {
                let divisor = lower[column][column];
                if divisor.abs() < 1e-12 {
                    return None;
                }
                lower[row][column] = sum / divisor;
            }
        }
    }
    Some(())
}
