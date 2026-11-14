//! Gaussian smearing and moments of sampled densities of states.

use xtalforge_core::Error;

use super::super::validate::invalid;

/// Square root of two pi, used by the Gaussian kernel.
const SQRT_TWO_PI: f64 = 2.506_628_274_631_000_7;

/// Deposits Gaussian-broadened levels on an energy grid.
///
/// Each level adds a unit-area Gaussian, so the integral of the result over a
/// wide grid is the number of levels. The grid must have at least two strictly
/// increasing, finite points and `sigma` must be positive.
pub fn smear_levels(levels: &[f64], grid: &[f64], sigma: f64) -> Result<Vec<f64>, Error> {
    validate_grid(grid)?;
    validate_sigma(sigma)?;
    if levels.iter().any(|level| !level.is_finite()) {
        return Err(invalid("a level energy is not finite"));
    }
    let mut values = vec![0.0; grid.len()];
    for level in levels {
        for (slot, energy) in values.iter_mut().zip(grid) {
            *slot += gaussian_kernel(energy - level, sigma);
        }
    }
    Ok(values)
}

/// Smooths a sampled DOS with a Gaussian of width `sigma`.
///
/// The convolution uses the trapezoid measure of the source grid, so a
/// normalized input stays normalized except for the tails that leave the grid.
/// Both arrays must have the same length, the grid must increase strictly and
/// `sigma` must be positive.
pub fn smear_dos(energies: &[f64], values: &[f64], sigma: f64) -> Result<Vec<f64>, Error> {
    validate_grid(energies)?;
    validate_sigma(sigma)?;
    if values.len() != energies.len() {
        return Err(invalid(format!(
            "the DOS has {} energies but {} values",
            energies.len(),
            values.len()
        )));
    }
    if values.iter().any(|value| !value.is_finite()) {
        return Err(invalid("a DOS value is not finite"));
    }
    let measures = trapezoid_measures(energies);
    let mut smoothed = vec![0.0; energies.len()];
    for ((center, value), measure) in energies.iter().zip(values).zip(&measures) {
        if *value == 0.0 {
            continue;
        }
        let weighted = measure * value;
        for (slot, energy) in smoothed.iter_mut().zip(energies) {
            *slot += weighted * gaussian_kernel(energy - center, sigma);
        }
    }
    Ok(smoothed)
}

/// Centre of a sampled DOS, `integral(E * DOS) / integral(DOS)`, in eV.
///
/// Uses the trapezoid rule. Energies must increase strictly, values must be
/// finite and non-negative, and the total weight must be positive.
pub fn d_band_centre(energies: &[f64], values: &[f64]) -> Result<f64, Error> {
    validate_grid(energies)?;
    if values.len() != energies.len() {
        return Err(invalid(format!(
            "the DOS has {} energies but {} values",
            energies.len(),
            values.len()
        )));
    }
    if values
        .iter()
        .any(|value| !value.is_finite() || *value < 0.0)
    {
        return Err(invalid("a DOS value must be finite and non-negative"));
    }
    let mut numerator = 0.0;
    let mut denominator = 0.0;
    for (energy_pair, value_pair) in energies.windows(2).zip(values.windows(2)) {
        if let ([left, right], [low, high]) = (energy_pair, value_pair) {
            let width = right - left;
            numerator += 0.5 * (left * low + right * high) * width;
            denominator += 0.5 * (low + high) * width;
        }
    }
    if denominator <= 0.0 {
        return Err(invalid("the DOS has zero total weight"));
    }
    Ok(numerator / denominator)
}

/// Normalized Gaussian kernel at a distance.
fn gaussian_kernel(distance: f64, sigma: f64) -> f64 {
    let scaled = distance / sigma;
    (-0.5 * scaled * scaled).exp() / (sigma * SQRT_TWO_PI)
}

/// Checks a grid for length, finiteness and strict growth.
fn validate_grid(grid: &[f64]) -> Result<(), Error> {
    if grid.len() < 2 {
        return Err(invalid("a DOS grid needs at least two points"));
    }
    if grid.iter().any(|value| !value.is_finite()) {
        return Err(invalid("a DOS grid value is not finite"));
    }
    for pair in grid.windows(2) {
        if let [left, right] = pair
            && right <= left
        {
            return Err(invalid("a DOS grid must increase strictly"));
        }
    }
    Ok(())
}

/// Checks the Gaussian width.
pub(super) fn validate_sigma(sigma: f64) -> Result<(), Error> {
    if !sigma.is_finite() || sigma <= 0.0 {
        return Err(invalid(format!(
            "the Gaussian width must be finite and positive, got {sigma}"
        )));
    }
    Ok(())
}

/// Trapezoid integration weights of a grid.
fn trapezoid_measures(grid: &[f64]) -> Vec<f64> {
    let mut measures = Vec::with_capacity(grid.len());
    for pair in grid.windows(2) {
        if let [left, right] = pair {
            let half = 0.5 * (right - left);
            if let Some(last) = measures.last_mut() {
                *last += half;
            } else {
                measures.push(half);
            }
            measures.push(half);
        }
    }
    measures
}

/// Trapezoid integral of sampled values.
pub(super) fn integrate(energies: &[f64], values: &[f64]) -> Result<f64, Error> {
    validate_grid(energies)?;
    if values.len() != energies.len() {
        return Err(invalid(format!(
            "the DOS has {} energies but {} values",
            energies.len(),
            values.len()
        )));
    }
    if values.iter().any(|value| !value.is_finite()) {
        return Err(invalid("a DOS value is not finite"));
    }
    let mut total = 0.0;
    for (energy_pair, value_pair) in energies.windows(2).zip(values.windows(2)) {
        if let ([left, right], [low, high]) = (energy_pair, value_pair) {
            total += 0.5 * (low + high) * (right - left);
        }
    }
    Ok(total)
}
