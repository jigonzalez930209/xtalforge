//! Model assembly and linear solve for the Pawley fit (Task 7.6.2).

#![allow(clippy::indexing_slicing)] // bounded local matrices

use xtalforge_io::xy::PatternPoint;

use crate::profiles::Profile;

use super::PawleyOptions;
use super::cell::Family;
use super::design::{Design, Linear};
use super::esd::invert;

/// Layout of the non-linear parameter vector.
#[derive(Clone)]
pub(super) struct Layout {
    family: Family,
    base: [f64; 6],
    fixed_fwhm: f64,
    free_cell: usize,
    refine_zero: bool,
    refine_size: bool,
    refine_strain: bool,
}

impl Layout {
    pub(super) fn new(family: Family, base: [f64; 6], options: &PawleyOptions) -> Self {
        Self {
            family,
            base,
            fixed_fwhm: options.fwhm,
            free_cell: if options.refine_cell {
                family.free_indices().len()
            } else {
                0
            },
            refine_zero: options.refine_zero,
            refine_size: options.refine_size,
            refine_strain: options.refine_strain,
        }
    }

    pub(super) fn len(&self) -> usize {
        self.free_cell
            + usize::from(self.refine_zero)
            + usize::from(self.refine_size)
            + usize::from(self.refine_strain)
    }

    pub(super) fn initial(&self, options: &PawleyOptions) -> Vec<f64> {
        let mut values = if self.free_cell > 0 {
            self.family.pack(self.base)
        } else {
            Vec::new()
        };
        if self.refine_zero {
            values.push(0.0);
        }
        if self.refine_size {
            values.push(options.fwhm);
        }
        if self.refine_strain {
            values.push(0.0);
        }
        values
    }

    pub(super) fn cell(&self, parameters: &[f64]) -> [f64; 6] {
        if self.free_cell == 0 {
            // With every cell number fixed the base cell is the model; without
            // this branch the family unpacking would read an empty vector and
            // collapse the cell to zero.
            return self.base;
        }
        self.family
            .unpack(self.base, parameters.get(..self.free_cell).unwrap_or(&[]))
    }

    pub(super) fn zero(&self, parameters: &[f64]) -> f64 {
        if self.refine_zero {
            parameters.get(self.free_cell).copied().unwrap_or(0.0)
        } else {
            0.0
        }
    }

    pub(super) fn size_value(&self, parameters: &[f64]) -> f64 {
        let index = self.free_cell + usize::from(self.refine_zero);
        if self.refine_size {
            parameters.get(index).copied().unwrap_or(self.fixed_fwhm)
        } else {
            self.fixed_fwhm
        }
    }

    pub(super) fn strain_value(&self, parameters: &[f64]) -> f64 {
        let index = self.free_cell + usize::from(self.refine_zero) + usize::from(self.refine_size);
        if self.refine_strain {
            parameters.get(index).copied().unwrap_or(0.0)
        } else {
            0.0
        }
    }

    pub(super) fn steps(&self) -> Vec<f64> {
        let mut steps = vec![0.05; self.free_cell];
        if self.refine_zero {
            steps.push(0.05);
        }
        if self.refine_size {
            steps.push(0.02);
        }
        if self.refine_strain {
            steps.push(0.02);
        }
        steps
    }

    pub(super) fn errors(
        &self,
        parameters: &[f64],
        design: &Design,
        linear: &Linear,
        variance: f64,
    ) -> ([f64; 6], f64) {
        let mut cell_esds = [0.0; 6];
        let mut zero_esd = 0.0;
        let count = self.len();
        if count > 0 {
            let jacobian = design.nonlinear_jacobian(parameters, linear);
            let mut normal = vec![vec![0.0; count]; count];
            for row in &jacobian {
                for first in 0..count {
                    for second in 0..count {
                        normal[first][second] += row[first] * row[second];
                    }
                }
            }
            if let Some(covariance) = invert(&normal) {
                for (slot, index) in self.family.free_indices().iter().enumerate() {
                    if let Some(value) = cell_esds.get_mut(*index) {
                        *value = (covariance[slot][slot] * variance).max(0.0).sqrt();
                    }
                }
                if self.refine_zero {
                    zero_esd = (covariance[self.free_cell][self.free_cell] * variance)
                        .max(0.0)
                        .sqrt();
                }
            }
        }
        (cell_esds, zero_esd)
    }
}

pub(super) fn chebyshev(order: usize, x: f64) -> Vec<f64> {
    let mut terms = Vec::with_capacity(order + 1);
    terms.push(1.0);
    if order >= 1 {
        terms.push(x);
    }
    for index in 2..=order {
        let previous = terms.get(index - 1).copied().unwrap_or(0.0);
        let before = terms.get(index - 2).copied().unwrap_or(0.0);
        terms.push(2.0 * x * previous - before);
    }
    terms
}

pub(super) fn profile_value(options: &PawleyOptions, center: f64, fwhm: f64, x: f64) -> f64 {
    let shape = Profile {
        kind: options.kind,
        center,
        fwhm,
        eta: options.eta,
        height: 1.0,
    }
    .value(x);
    // Area-normalised Gaussian, so a fitted intensity is the integrated
    // intensity of the line and matches the PXRD reflection scale.
    let width = fwhm.max(1e-9);
    shape * 2.0 / width * (std::f64::consts::LN_2 / std::f64::consts::PI).sqrt()
}

pub(super) fn weight(point: &PatternPoint) -> f64 {
    match point.uncertainty {
        Some(value) if value.is_finite() && value > 0.0 => 1.0 / (value * value),
        _ => 1.0 / point.intensity.abs().max(1.0),
    }
}

pub(super) fn factors(
    observed: &[PatternPoint],
    fit: &[f64],
    parameters: usize,
) -> (f64, f64, f64, f64) {
    let mut numerator = 0.0;
    let mut denominator = 0.0;
    let mut weighted = 0.0;
    let mut weighted_observed = 0.0;
    for (point, value) in observed.iter().zip(fit.iter()) {
        let weight = weight(point);
        numerator += (point.intensity - value).abs();
        denominator += point.intensity.abs();
        weighted += weight * (point.intensity - value).powi(2);
        weighted_observed += weight * point.intensity.powi(2);
    }
    let degrees = observed.len().saturating_sub(parameters).max(1) as f64;
    let rp = if denominator > 0.0 {
        numerator / denominator
    } else {
        0.0
    };
    let rwp = if weighted_observed > 0.0 {
        (weighted / weighted_observed).sqrt()
    } else {
        0.0
    };
    let rexp = if weighted_observed > 0.0 {
        (degrees / weighted_observed).sqrt()
    } else {
        0.0
    };
    (rp, rwp, rexp, weighted / degrees)
}
