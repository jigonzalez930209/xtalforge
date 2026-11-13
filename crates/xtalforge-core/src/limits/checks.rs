//! Arithmetic estimates and the named size checks (Task 1.3, Task 15.2.5).

use crate::Error;
use crate::lattice::Lattice;

use super::{Budget, HESSIAN_MIB_PER_ATOM, Limits, Verdict, count};

/// Atoms drawn for a display range: cell atoms times the range volume.
#[must_use]
pub fn drawn_atoms(n_atoms: u64, ranges: [(f64, f64); 3]) -> u64 {
    let volume = ranges
        .iter()
        .map(|(low, high)| (high - low).max(0.0))
        .product::<f64>();
    (n_atoms as f64 * volume).round() as u64
}

/// `|det P|` for a 3x3 transformation matrix.
///
/// # Errors
///
/// Returns [`Error::InvalidInput`] for a singular matrix or non-finite values.
pub fn multiplier(p: [[f64; 3]; 3]) -> Result<u64, Error> {
    if !p.iter().flatten().all(|value| value.is_finite()) {
        return Err(Error::InvalidInput(
            "transformation matrix contains non-finite values".to_owned(),
        ));
    }
    let [[a, b, c], [d, e, f], [g, h, i]] = p;
    let determinant = a * (e * i - f * h) - b * (d * i - f * g) + c * (d * h - e * g);
    if determinant.abs() < 1e-9 {
        return Err(Error::InvalidInput(
            "transformation matrix is singular".to_owned(),
        ));
    }
    Ok(determinant.abs().round() as u64)
}

/// Atoms in the cell a transformation builds from `n_atoms`.
///
/// # Errors
///
/// Returns the errors of [`multiplier`].
pub fn supercell_atoms(n_atoms: u64, p: [[f64; 3]; 3]) -> Result<u64, Error> {
    Ok(n_atoms.saturating_mul(multiplier(p)?))
}

/// Grid shape from cell lengths and a target spacing.
///
/// # Errors
///
/// Returns [`Error::InvalidInput`] for a non-positive or non-finite spacing.
pub fn grid_shape(lattice: &Lattice, spacing: f64) -> Result<[u32; 3], Error> {
    if !spacing.is_finite() || spacing <= 0.0 {
        return Err(Error::InvalidInput(format!(
            "grid spacing must be positive and finite, got {spacing}"
        )));
    }
    let spacing = spacing.max(0.05);
    Ok(lattice
        .lengths()
        .map(|length| (length / spacing).round().max(8.0) as u32))
}

/// Points in a porosity grid.
///
/// # Errors
///
/// Returns the errors of [`grid_shape`].
pub fn grid_points(lattice: &Lattice, spacing: f64) -> Result<u64, Error> {
    let [a, b, c] = grid_shape(lattice, spacing)?;
    Ok(u64::from(a) * u64::from(b) * u64::from(c))
}

fn cube_that_fits(per_cell: f64, limit: u64, cells: &str) -> String {
    let k = ((limit as f64 / per_cell.max(1.0)).powf(1.0 / 3.0) + 1e-9).floor() as u64;
    if k < 1 {
        format!("one cell alone is {} atoms", count(per_cell.round() as u64))
    } else {
        format!("{k} x {k} x {k}{cells} is the largest that fits")
    }
}

/// Checks the atoms a display range would draw.
#[must_use]
pub fn check_drawn(limits: &Limits, n_atoms: u64, ranges: [(f64, f64); 3]) -> Verdict {
    let estimate = drawn_atoms(n_atoms, ranges);
    limits.judge(
        Budget::Drawn,
        estimate,
        "the picture would draw {} atoms",
        |limit| cube_that_fits(n_atoms as f64, limit, " cells"),
    )
}

/// Checks a supercell transformation.
///
/// # Errors
///
/// Returns the errors of [`supercell_atoms`].
pub fn check_supercell(limits: &Limits, n_atoms: u64, p: [[f64; 3]; 3]) -> Result<Verdict, Error> {
    let estimate = supercell_atoms(n_atoms, p)?;
    Ok(limits.judge(
        Budget::Supercell,
        estimate,
        "the new cell would have {} atoms",
        |limit| cube_that_fits(n_atoms as f64, limit, ""),
    ))
}

/// Checks a porosity grid.
///
/// # Errors
///
/// Returns the errors of [`grid_points`].
pub fn check_grid(limits: &Limits, lattice: &Lattice, spacing: f64) -> Result<Verdict, Error> {
    let estimate = grid_points(lattice, spacing)?;
    Ok(limits.judge(
        Budget::Grid,
        estimate,
        "the grid would have {} points",
        |limit| {
            let mut step = (spacing / 0.05).ceil() * 0.05;
            step = step.max(0.05);
            while grid_points(lattice, step).unwrap_or(u64::MAX) > limit && step < 5.0 {
                step += 0.05;
            }
            format!(
                "at {step:.2} A it would have {}",
                grid_points(lattice, step)
                    .map(count)
                    .unwrap_or_else(|_| "an invalid count".to_owned())
            )
        },
    ))
}

/// Checks a phonon supercell.
///
/// # Errors
///
/// Returns the errors of [`supercell_atoms`].
pub fn check_phonons(limits: &Limits, n_atoms: u64, p: [[f64; 3]; 3]) -> Result<Verdict, Error> {
    let estimate = supercell_atoms(n_atoms, p)?;
    Ok(limits.judge(
        Budget::Phonons,
        estimate,
        "the phonon supercell would have {} atoms",
        |limit| {
            format!(
                "{}, or Gamma only needs none",
                cube_that_fits(n_atoms as f64, limit, "")
            )
        },
    ))
}

/// Checks an analytic Hessian over `n_atoms` free atoms.
#[must_use]
pub fn check_hessian(limits: &Limits, n_atoms: u64) -> Verdict {
    limits.judge(
        Budget::Hessian,
        n_atoms,
        &format!(
            "an analytic Hessian over {{}} atoms would hold about {:.1} GB",
            n_atoms as f64 * HESSIAN_MIB_PER_ATOM / 1024.0
        ),
        |limit| {
            format!(
                "up to {} atoms fit, and finite differences need no more than one evaluation does",
                count(limit)
            )
        },
    )
}
