//! Free-function wrappers for the lattice API.

use super::Lattice;
use crate::Error;

/// Convenience wrapper for `Lattice::from_parameters`.
///
/// # Errors
///
/// Returns the errors of [`Lattice::from_parameters`].
pub fn lattice_info(
    a: f64,
    b: f64,
    c: f64,
    alpha: f64,
    beta: f64,
    gamma: f64,
) -> Result<Lattice, Error> {
    Lattice::from_parameters(a, b, c, alpha, beta, gamma)
}

/// Minimum-image vector and distance for a fractional displacement.
///
/// # Errors
///
/// Returns the errors of [`Lattice::from_parameters`].
pub fn lattice_minimum_image(
    a: f64,
    b: f64,
    c: f64,
    alpha: f64,
    beta: f64,
    gamma: f64,
    delta: [f64; 3],
) -> Result<([f64; 3], f64), Error> {
    Ok(Lattice::from_parameters(a, b, c, alpha, beta, gamma)?.minimum_image(delta))
}
