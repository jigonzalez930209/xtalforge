//! Equation of state through a set of scaled volumes (Task 4.4.5).
//!
//! Ported from Crystal-Builder `xtal/ff/eos.py`, reduced to a deterministic
//! in-crate fit.  Two forms are offered:
//!
//! * third-order Birch-Murnaghan, what papers usually report;
//! * Murnaghan, the older integrated form and a useful check.
//!
//! The fit refuses fewer than four finite points, a minimum outside the
//! scanned volumes, or a non-positive bulk modulus rather than returning a
//! plausible number from an extrapolation.  Volumes are Angstrom cubed,
//! energies eV for the whole cell, and [`EosFit::b0`] is GPa.

use serde::{Deserialize, Serialize};

use crate::Error;
use crate::math::solve_linear;
use crate::optimize::GPA_TO_EV_PER_A3;

/// Minimum number of points the four-parameter fit accepts.
pub const MINIMUM_EOS_POINTS: usize = 4;

/// Which energy-volume form is fitted.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default, Serialize, Deserialize)]
pub enum EosForm {
    /// Third-order Birch-Murnaghan.
    #[default]
    BirchMurnaghan,
    /// Integrated Murnaghan form.
    Murnaghan,
}

/// One fitted equation of state.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct EosFit {
    /// Form that was fitted.
    pub form: EosForm,
    /// Equilibrium volume in Angstrom cubed.
    pub v0: f64,
    /// Energy at the equilibrium volume in eV.
    pub e0: f64,
    /// Bulk modulus in GPa.
    pub b0: f64,
    /// Pressure derivative of the bulk modulus, dimensionless.
    pub b0_prime: f64,
    /// Root-mean-square residual in eV.
    pub rms: f64,
    /// Number of finite points used.
    pub n_points: usize,
}

impl EosFit {
    /// Fitted energy at a volume, in eV.
    #[must_use]
    pub fn energy(&self, volume: f64) -> f64 {
        model(
            self.form,
            volume,
            self.e0,
            self.v0,
            self.b0 * GPA_TO_EV_PER_A3,
            self.b0_prime,
        )
    }
}

/// Fits an equation of state to volume-energy points.
///
/// Non-finite pairs are left out, as an unconverged scan point is not a
/// number.  The points must bracket the minimum: if the lowest energy is at
/// the first or last volume, `V0` would be an extrapolation and the fit is
/// refused.
///
/// # Errors
///
/// Returns [`Error::InvalidInput`] for fewer than
/// [`MINIMUM_EOS_POINTS`] finite points, non-positive volumes, a minimum
/// outside the scan, a curve that does not turn upwards, or a fit that does
/// not converge.
pub fn fit_eos(volumes: &[f64], energies: &[f64], form: EosForm) -> Result<EosFit, Error> {
    if volumes.len() != energies.len() {
        return Err(Error::InvalidInput(format!(
            "got {} volumes and {} energies",
            volumes.len(),
            energies.len()
        )));
    }
    let mut points: Vec<(f64, f64)> = volumes
        .iter()
        .zip(energies)
        .filter(|(volume, energy)| volume.is_finite() && energy.is_finite() && **volume > 0.0)
        .map(|(volume, energy)| (*volume, *energy))
        .collect();
    if points.len() < MINIMUM_EOS_POINTS {
        return Err(Error::InvalidInput(format!(
            "{} finite points; an equation of state needs at least {MINIMUM_EOS_POINTS}",
            points.len()
        )));
    }
    points.sort_by(|first, second| first.0.total_cmp(&second.0));
    let lowest = points
        .iter()
        .enumerate()
        .min_by(|(_, first), (_, second)| first.1.total_cmp(&second.1))
        .map(|(index, _)| index)
        .unwrap_or(0);
    if lowest == 0 || lowest + 1 == points.len() {
        return Err(Error::InvalidInput(
            "the lowest energy is at an end of the volume scan; the minimum is not bracketed"
                .to_owned(),
        ));
    }

    let initial = quadratic_start(&points)?;
    let mut parameters = initial;
    let mut sse = sum_squares(&points, form, parameters);
    let mut damping = 1e-3_f64;
    for _iteration in 0..200 {
        let jacobian = jacobian(&points, form, parameters);
        let residuals = residuals(&points, form, parameters);
        let mut normal = vec![vec![0.0; 4]; 4];
        let mut right = vec![0.0; 4];
        for (row, residual) in jacobian.iter().zip(&residuals) {
            for i in 0..4 {
                let first = row.get(i).copied().unwrap_or(0.0);
                if let Some(slot) = right.get_mut(i) {
                    *slot -= first * residual;
                }
                for j in 0..4 {
                    let second = row.get(j).copied().unwrap_or(0.0);
                    if let Some(slot) = normal.get_mut(i).and_then(|values| values.get_mut(j)) {
                        *slot += first * second;
                    }
                }
            }
        }
        for i in 0..4 {
            if let Some(slot) = normal.get_mut(i).and_then(|row| row.get_mut(i)) {
                *slot += damping * (*slot).abs().max(1e-12);
            }
        }
        let Some(delta) = solve_linear(&mut normal, &mut right) else {
            damping *= 10.0;
            if damping > 1e12 {
                break;
            }
            continue;
        };
        let mut trial = parameters;
        let mut valid = true;
        for (index, value) in trial.iter_mut().enumerate() {
            *value += delta.get(index).copied().unwrap_or(0.0);
        }
        if trial[1] <= 0.0
            || trial[2] <= 0.0
            || (form == EosForm::Murnaghan && (trial[3] - 1.0).abs() < 1e-6)
        {
            valid = false;
        }
        if valid {
            let trial_sse = sum_squares(&points, form, trial);
            if trial_sse.is_finite() && trial_sse <= sse {
                let change = (sse - trial_sse).abs();
                parameters = trial;
                sse = trial_sse;
                damping = (damping * 0.3).max(1e-12);
                if change <= 1e-14 * sse.max(1.0) {
                    break;
                }
                continue;
            }
        }
        damping *= 10.0;
        if damping > 1e12 {
            break;
        }
    }

    let [e0, v0, b0, b0_prime] = parameters;
    let lowest_volume = points.first().map(|point| point.0).unwrap_or(0.0);
    let highest_volume = points.last().map(|point| point.0).unwrap_or(0.0);
    if !(lowest_volume..=highest_volume).contains(&v0) || b0 <= 0.0 {
        return Err(Error::InvalidInput(
            "the equation of state puts the minimum outside the scan; widen it".to_owned(),
        ));
    }
    let rms = (sse / points.len() as f64).sqrt();
    Ok(EosFit {
        form,
        v0,
        e0,
        b0: b0 / GPA_TO_EV_PER_A3,
        b0_prime,
        rms,
        n_points: points.len(),
    })
}

/// Quadratic start: vertex and curvature give `V0`, `E0` and `B0`.
fn quadratic_start(points: &[(f64, f64)]) -> Result<[f64; 4], Error> {
    let n = points.len() as f64;
    let mut s = [0.0; 5];
    let mut t = [0.0; 3];
    for (volume, energy) in points {
        let mut power = 1.0;
        for value in &mut s {
            *value += power;
            power *= volume;
        }
        t[0] += energy;
        t[1] += energy * volume;
        t[2] += energy * volume * volume;
    }
    let mut matrix = vec![
        vec![s[4], s[3], s[2]],
        vec![s[3], s[2], s[1]],
        vec![s[2], s[1], n],
    ];
    let mut right = vec![t[2], t[1], t[0]];
    let solution = solve_linear(&mut matrix, &mut right).ok_or_else(|| {
        Error::InvalidInput(
            "the volumes are degenerate; the quadratic start is singular".to_owned(),
        )
    })?;
    let a = solution.first().copied().unwrap_or(0.0);
    let b = solution.get(1).copied().unwrap_or(0.0);
    let c = solution.get(2).copied().unwrap_or(0.0);
    if a <= 0.0 {
        return Err(Error::InvalidInput(
            "the energy does not curve upwards, so there is no minimum to fit".to_owned(),
        ));
    }
    let v0 = -b / (2.0 * a);
    let e0 = a * v0 * v0 + b * v0 + c;
    let b0 = 2.0 * a * v0;
    if v0 <= 0.0 || b0 <= 0.0 {
        return Err(Error::InvalidInput(
            "the quadratic start has no positive minimum".to_owned(),
        ));
    }
    Ok([e0, v0, b0, 4.0])
}

fn residuals(points: &[(f64, f64)], form: EosForm, parameters: [f64; 4]) -> Vec<f64> {
    let [e0, v0, b0, b0_prime] = parameters;
    points
        .iter()
        .map(|(volume, energy)| model(form, *volume, e0, v0, b0, b0_prime) - energy)
        .collect()
}

fn jacobian(points: &[(f64, f64)], form: EosForm, parameters: [f64; 4]) -> Vec<[f64; 4]> {
    let mut result = vec![[0.0; 4]; points.len()];
    for (index, (volume, _)) in points.iter().enumerate() {
        for (column, value) in parameters.iter().enumerate() {
            let step = value.abs().max(1.0) * 1e-6;
            let mut high = parameters;
            let mut low = parameters;
            if let Some(value) = high.get_mut(column) {
                *value += step;
            }
            if let Some(value) = low.get_mut(column) {
                *value -= step;
            }
            let up = model(form, *volume, high[0], high[1], high[2], high[3]);
            let down = model(form, *volume, low[0], low[1], low[2], low[3]);
            if let Some(row) = result.get_mut(index)
                && let Some(slot) = row.get_mut(column)
            {
                *slot = (up - down) / (2.0 * step);
            }
        }
    }
    result
}

fn sum_squares(points: &[(f64, f64)], form: EosForm, parameters: [f64; 4]) -> f64 {
    residuals(points, form, parameters)
        .iter()
        .map(|value| value * value)
        .sum()
}

/// Third-order Birch-Murnaghan, `b0` in energy per volume.
#[must_use]
pub fn birch_murnaghan(volume: f64, e0: f64, v0: f64, b0: f64, b0_prime: f64) -> f64 {
    let eta = (v0 / volume).powf(2.0 / 3.0) - 1.0;
    e0 + 9.0 * v0 * b0 / 16.0 * (eta.powi(3) * b0_prime + eta.powi(2) * (6.0 - 4.0 * (eta + 1.0)))
}

/// Murnaghan's integrated form, `b0` in energy per volume.
#[must_use]
pub fn murnaghan(volume: f64, e0: f64, v0: f64, b0: f64, b0_prime: f64) -> f64 {
    let exponent = b0_prime - 1.0;
    e0 + b0 * volume / b0_prime * ((v0 / volume).powf(b0_prime) / exponent + 1.0)
        - b0 * v0 / exponent
}

fn model(form: EosForm, volume: f64, e0: f64, v0: f64, b0: f64, b0_prime: f64) -> f64 {
    match form {
        EosForm::BirchMurnaghan => birch_murnaghan(volume, e0, v0, b0, b0_prime),
        EosForm::Murnaghan => murnaghan(volume, e0, v0, b0, b0_prime),
    }
}
