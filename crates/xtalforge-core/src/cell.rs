//! Cell geometry: volumes from the six lattice parameters.
//!
//! Lengths are in Angstrom and angles in degrees at the API boundary, as
//! everywhere in xtalforge.

use crate::Error;

/// Volume of the cell described by `a, b, c, alpha, beta, gamma`.
///
/// Angles are in degrees. The result is in cubic Angstrom.
///
/// # Errors
///
/// Returns [`Error::InvalidInput`] when a length is not positive, an angle is
/// not strictly between 0 and 180 degrees, a value is not finite, or the
/// parameters describe a degenerate (zero-volume) cell.
pub fn cell_volume(
    a: f64,
    b: f64,
    c: f64,
    alpha: f64,
    beta: f64,
    gamma: f64,
) -> Result<f64, Error> {
    let values = [a, b, c, alpha, beta, gamma];
    if !values.iter().all(|value| value.is_finite()) {
        return Err(Error::InvalidInput(
            "cell parameters must be finite numbers".to_owned(),
        ));
    }
    if a <= 0.0 || b <= 0.0 || c <= 0.0 {
        return Err(Error::InvalidInput(
            "cell lengths must be strictly positive".to_owned(),
        ));
    }
    if !(valid_angle(alpha) && valid_angle(beta) && valid_angle(gamma)) {
        return Err(Error::InvalidInput(
            "cell angles must lie strictly between 0 and 180 degrees".to_owned(),
        ));
    }

    let cos_alpha = alpha.to_radians().cos();
    let cos_beta = beta.to_radians().cos();
    let cos_gamma = gamma.to_radians().cos();
    let factor = 1.0 - cos_alpha * cos_alpha - cos_beta * cos_beta - cos_gamma * cos_gamma
        + 2.0 * cos_alpha * cos_beta * cos_gamma;

    if factor <= 0.0 || !factor.is_finite() {
        return Err(Error::InvalidInput(
            "cell parameters describe a degenerate cell".to_owned(),
        ));
    }

    Ok(a * b * c * factor.sqrt())
}

fn valid_angle(angle: f64) -> bool {
    angle > 0.0 && angle < 180.0
}
