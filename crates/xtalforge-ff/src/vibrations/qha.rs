//! Quasi-harmonic free energy from sampled volumes (Task 4.8.5, core subset).
//!
//! The full task scans volumes, relaxes each one, computes vibrations at every
//! converged volume and minimizes `F(V, T) = E(V) + F_vib(V, T)` at every
//! temperature.  This module implements the free-energy evaluation and a
//! three-point parabolic refinement of the minimum for volumes the caller has
//! already prepared.  It does **not** relax volumes, fit an equation of state,
//! write checkpoints or compute thermal expansion and bulk modulus; those are
//! recorded as pending rather than approximated.
//!
//! The frequencies at each volume are used as supplied.  Linear interpolation
//! between volumes is therefore the caller's choice in how densely it samples
//! the volume axis; this code interpolates only the free-energy minimum, not
//! each mode frequency.  That is stated here so a smooth curve is not mistaken
//! for a mode-Grueneisen model.

use serde::{Deserialize, Serialize};

use crate::Error;

use super::thermo::{ThermoOptions, thermochemistry};

/// One volume point with its static energy and frequencies.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct QhaPoint {
    /// Volume in Angstrom^3.
    pub volume: f64,
    /// Static (electronic or force-field) energy at this volume, in eV.
    pub static_energy: f64,
    /// Harmonic frequencies at this volume, in cm^-1.
    pub frequencies: Vec<f64>,
}

/// Free energy at one sampled volume.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct QhaFreeEnergy {
    /// Volume in Angstrom^3.
    pub volume: f64,
    /// Helmholtz free energy `E + ZPE + U_vib - T S` in eV.
    pub free_energy: f64,
}

/// Result of a quasi-harmonic free-energy evaluation.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct QhaResult {
    /// One entry per input volume, in input order.
    pub points: Vec<QhaFreeEnergy>,
    /// Volume of the lowest free energy, in Angstrom^3.
    pub minimum_volume: f64,
    /// Lowest free energy, in eV.
    pub minimum_free_energy: f64,
}

/// Evaluates the quasi-harmonic free energy at each supplied volume.
///
/// Volumes must be positive, finite and strictly increasing; at least two
/// points are required.  With three or more points the lowest sampled value is
/// refined by a parabola through its neighbours when the best point is
/// interior; otherwise the best sample is returned unchanged.  The result is
/// Helmholtz free energy (no pressure-volume term); pass the static energy with
/// the desired reference and treat [`ThermoOptions::pressure`] as the
/// documented standard state for the vibrational entropy.
///
/// # Errors
///
/// Returns [`Error::InvalidInput`] for bad options, fewer than two points,
/// non-positive or non-increasing volumes, a non-finite energy or an invalid
/// frequency list.
pub fn quasi_harmonic_free_energy(
    points: &[QhaPoint],
    options: &ThermoOptions,
) -> Result<QhaResult, Error> {
    options.validate()?;
    if points.len() < 2 {
        return Err(Error::InvalidInput(
            "quasi-harmonic free energy needs at least two volume points".to_owned(),
        ));
    }
    let mut evaluated = Vec::with_capacity(points.len());
    let mut previous: Option<f64> = None;
    for (index, point) in points.iter().enumerate() {
        if !point.volume.is_finite() || point.volume <= 0.0 {
            return Err(Error::InvalidInput(format!(
                "volume point {index} must be positive and finite, got {}",
                point.volume
            )));
        }
        if let Some(last) = previous
            && point.volume <= last
        {
            return Err(Error::InvalidInput(format!(
                "volume points must be strictly increasing; point {index} is {} after {last}",
                point.volume
            )));
        }
        previous = Some(point.volume);
        if !point.static_energy.is_finite() {
            return Err(Error::InvalidInput(format!(
                "static energy at volume point {index} is not finite"
            )));
        }
        let thermo = thermochemistry(&point.frequencies, options)?;
        let free_energy = point.static_energy + thermo.zpe + thermo.thermal_u
            - options.temperature * thermo.entropy;
        if !free_energy.is_finite() {
            return Err(Error::InvalidInput(format!(
                "free energy at volume point {index} is not finite"
            )));
        }
        evaluated.push(QhaFreeEnergy {
            volume: point.volume,
            free_energy,
        });
    }

    let best = evaluated
        .iter()
        .enumerate()
        .min_by(|(first_index, first), (second_index, second)| {
            first
                .free_energy
                .partial_cmp(&second.free_energy)
                .unwrap_or(std::cmp::Ordering::Equal)
                .then(first_index.cmp(second_index))
        })
        .map(|(index, point)| (index, point.clone()))
        .ok_or_else(|| Error::InvalidInput("no volume points were evaluated".to_owned()))?;
    let (index, sampled) = best;
    let mut minimum_volume = sampled.volume;
    let mut minimum_free_energy = sampled.free_energy;
    if index > 0
        && index + 1 < evaluated.len()
        && let (Some(left), Some(right)) = (evaluated.get(index - 1), evaluated.get(index + 1))
    {
        let x1 = left.volume;
        let x2 = sampled.volume;
        let x3 = right.volume;
        let y1 = left.free_energy;
        let y2 = sampled.free_energy;
        let y3 = right.free_energy;
        let denominator = x1 * (y2 - y3) + x2 * (y3 - y1) + x3 * (y1 - y2);
        if denominator.abs() > f64::EPSILON {
            let vertex = 0.5 * (x1 * x1 * (y2 - y3) + x2 * x2 * (y3 - y1) + x3 * x3 * (y1 - y2))
                / denominator;
            if vertex > x1 && vertex < x3 {
                let first = (vertex - x2) * (vertex - x3) / ((x1 - x2) * (x1 - x3));
                let second = (vertex - x1) * (vertex - x3) / ((x2 - x1) * (x2 - x3));
                let third = (vertex - x1) * (vertex - x2) / ((x3 - x1) * (x3 - x2));
                let value = y1 * first + y2 * second + y3 * third;
                if value.is_finite() && value <= sampled.free_energy {
                    minimum_volume = vertex;
                    minimum_free_energy = value;
                }
            }
        }
    }
    Ok(QhaResult {
        points: evaluated,
        minimum_volume,
        minimum_free_energy,
    })
}
