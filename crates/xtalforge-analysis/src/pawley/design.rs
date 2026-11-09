//! Design matrix, matching and linear solve for the Pawley fit (Task 7.6.2).

#![allow(clippy::indexing_slicing)] // bounded local matrices

use xtalforge_core::Error;
use xtalforge_io::xy::PatternPoint;

use crate::indexing::CellMetric;

use super::PawleyOptions;
use super::esd::{invert, solve};
use super::model::{Layout, chebyshev, profile_value, weight};

/// One allowed reflection position.
pub(super) struct ReflectionData {
    pub(super) hkl: [i32; 3],
    pub(super) d: f64,
    pub(super) two_theta: f64,
}

/// Linear parameters at one non-linear point.
#[derive(Clone)]
pub(super) struct Linear {
    pub(super) intensities: Vec<f64>,
    pub(super) background: Vec<f64>,
    pub(super) covariance: Vec<Vec<f64>>,
}

impl Linear {
    pub(super) fn value(&self, reflection: usize) -> f64 {
        self.intensities.get(reflection).copied().unwrap_or(0.0)
    }

    pub(super) fn esd(&self, index: usize, variance: f64) -> f64 {
        self.covariance
            .get(index)
            .and_then(|row| row.get(index))
            .map_or(0.0, |value| (value * variance).max(0.0).sqrt())
    }
}

pub(super) struct Design {
    observed: Vec<PatternPoint>,
    pub(super) reflections: Vec<ReflectionData>,
    widths: Vec<f64>,
    background: Vec<Vec<f64>>,
    active: Vec<Vec<usize>>,
    background_order: usize,
    weights: Vec<f64>,
    options: PawleyOptions,
    layout: Layout,
    hkls: Vec<[i32; 3]>,
    wavelength: f64,
}

impl Design {
    pub(super) fn new(
        observed: &[PatternPoint],
        hkls: &[[i32; 3]],
        options: &PawleyOptions,
        layout: &Layout,
        parameters: &[f64],
        wavelength: f64,
    ) -> Result<Self, Error> {
        let cell = layout.cell(parameters);
        let metric = CellMetric::new(cell[0], cell[1], cell[2], cell[3], cell[4], cell[5]);
        let zero = layout.zero(parameters);
        let size = layout.size_value(parameters);
        let strain = layout.strain_value(parameters);
        let first = observed.first().map_or(0.0, |point| point.two_theta);
        let last = observed.last().map_or(0.0, |point| point.two_theta);
        let mut reflections = Vec::new();
        for hkl in hkls {
            let Some(q) = metric.reciprocal_squared(*hkl) else {
                continue;
            };
            let d = 1.0 / q.sqrt();
            let sine = wavelength / (2.0 * d);
            if sine >= 1.0 {
                continue;
            }
            let two_theta = 2.0 * sine.asin().to_degrees() + zero;
            if two_theta < first - 5.0 || two_theta > last + 5.0 {
                continue;
            }
            reflections.push(ReflectionData {
                hkl: *hkl,
                d,
                two_theta,
            });
        }
        if reflections.is_empty() {
            return Err(Error::InvalidInput(
                "the cell predicts no reflection in the observed range".to_owned(),
            ));
        }
        let widths = reflections
            .iter()
            .map(|reflection| {
                let theta = (reflection.two_theta - zero).to_radians() / 2.0;
                let term = strain * theta.tan();
                (size * size + term * term).sqrt().max(0.01)
            })
            .collect::<Vec<_>>();
        let background_order = options.background_order;
        let width = last - first;
        let background = observed
            .iter()
            .map(|point| {
                let x = if width > 0.0 {
                    2.0 * (point.two_theta - first) / width - 1.0
                } else {
                    0.0
                };
                chebyshev(background_order, x)
            })
            .collect::<Vec<_>>();
        let active = observed
            .iter()
            .map(|point| {
                reflections
                    .iter()
                    .enumerate()
                    .filter(|(index, reflection)| {
                        widths.get(*index).is_some_and(|width| {
                            (point.two_theta - reflection.two_theta).abs() <= 6.0 * width
                        })
                    })
                    .map(|(index, _)| index)
                    .collect::<Vec<_>>()
            })
            .collect::<Vec<_>>();
        let values = observed.iter().map(weight).collect::<Vec<_>>();
        Ok(Self {
            observed: observed.to_vec(),
            reflections,
            widths,
            background,
            active,
            background_order,
            weights: values,
            options: *options,
            layout: layout.clone(),
            hkls: hkls.to_vec(),
            wavelength,
        })
    }

    pub(super) fn refit(&self, parameters: &[f64]) -> Result<Self, Error> {
        Self::new(
            &self.observed,
            &self.hkls,
            &self.options,
            &self.layout,
            parameters,
            self.wavelength,
        )
    }

    pub(super) fn parameter_count(&self) -> usize {
        self.reflections.len() + self.background_order + 1
    }

    pub(super) fn solve(&self) -> Result<Linear, Error> {
        let count = self.reflections.len() + self.background_order + 1;
        let reflection_count = self.reflections.len();
        let mut at_a = vec![vec![0.0; count]; count];
        let mut at_y = vec![0.0; count];
        for (row, point) in self.observed.iter().enumerate() {
            let mut entries: Vec<(usize, f64)> = Vec::new();
            for &index in self.active.get(row).into_iter().flatten() {
                let Some(reflection) = self.reflections.get(index) else {
                    continue;
                };
                let Some(width) = self.widths.get(index).copied() else {
                    continue;
                };
                entries.push((
                    index,
                    profile_value(&self.options, reflection.two_theta, width, point.two_theta),
                ));
            }
            if let Some(terms) = self.background.get(row) {
                for (offset, term) in terms.iter().enumerate() {
                    entries.push((reflection_count + offset, *term));
                }
            }
            let weight = self.weights.get(row).copied().unwrap_or(1.0);
            for &(first, first_value) in &entries {
                if let Some(slot) = at_y.get_mut(first) {
                    *slot += first_value * weight * point.intensity;
                }
                for &(second, second_value) in &entries {
                    if let Some(cell) = at_a.get_mut(first).and_then(|line| line.get_mut(second)) {
                        *cell += first_value * weight * second_value;
                    }
                }
            }
        }
        let scale = at_a
            .iter()
            .flat_map(|row| row.iter())
            .fold(0.0_f64, |value, entry| value.max(entry.abs()))
            .max(1.0);
        for (index, row) in at_a.iter_mut().enumerate() {
            if let Some(value) = row.get_mut(index) {
                *value += 1e-9 * scale;
            }
        }
        let Some(solution) = solve(&at_a, &at_y) else {
            return Err(Error::InvalidInput(
                "the Pawley normal matrix is singular".to_owned(),
            ));
        };
        let covariance = invert(&at_a).unwrap_or_else(|| vec![vec![0.0; count]; count]);
        Ok(Linear {
            intensities: solution.iter().take(reflection_count).copied().collect(),
            background: solution.iter().skip(reflection_count).copied().collect(),
            covariance,
        })
    }

    pub(super) fn model(&self, linear: &Linear) -> Vec<f64> {
        let mut values = vec![0.0; self.observed.len()];
        for (row, point) in self.observed.iter().enumerate() {
            let mut value = 0.0;
            for &index in self.active.get(row).into_iter().flatten() {
                let Some(reflection) = self.reflections.get(index) else {
                    continue;
                };
                let Some(width) = self.widths.get(index).copied() else {
                    continue;
                };
                value += linear.value(index)
                    * profile_value(&self.options, reflection.two_theta, width, point.two_theta);
            }
            if let Some(terms) = self.background.get(row) {
                for (offset, coefficient) in linear.background.iter().enumerate() {
                    if let Some(term) = terms.get(offset) {
                        value += coefficient * term;
                    }
                }
            }
            if let Some(slot) = values.get_mut(row) {
                *slot = value;
            }
        }
        values
    }

    pub(super) fn residuals(&self, linear: &Linear) -> Vec<f64> {
        let model = self.model(linear);
        self.observed
            .iter()
            .zip(model.iter())
            .enumerate()
            .map(|(row, (point, value))| {
                let weight = self.weights.get(row).copied().unwrap_or(1.0);
                weight.sqrt() * (point.intensity - value)
            })
            .collect()
    }

    pub(super) fn chi_squared(&self, parameters: &[f64]) -> Result<f64, Error> {
        let design = self.refit(parameters)?;
        let linear = design.solve()?;
        Ok(design
            .residuals(&linear)
            .iter()
            .map(|value| value * value)
            .sum())
    }

    pub(super) fn reduced_variance(&self, _parameters: &[f64]) -> f64 {
        let linear = self.solve().unwrap_or_else(|_| Linear {
            intensities: vec![0.0; self.reflections.len()],
            background: vec![0.0; self.background_order + 1],
            covariance: Vec::new(),
        });
        let residuals = self.residuals(&linear);
        let chi_squared: f64 = residuals.iter().map(|value| value * value).sum();
        chi_squared
            / self
                .observed
                .len()
                .saturating_sub(self.parameter_count())
                .max(1) as f64
    }

    pub(super) fn nonlinear_jacobian(&self, parameters: &[f64], linear: &Linear) -> Vec<Vec<f64>> {
        let count = parameters.len();
        let rows = self.observed.len();
        let mut jacobian = vec![vec![0.0; count]; rows];
        for column in 0..count {
            let value = parameters.get(column).copied().unwrap_or(0.0);
            let step = 1e-5 * value.abs().max(1.0);
            let mut plus = parameters.to_vec();
            let mut minus = parameters.to_vec();
            if let Some(slot) = plus.get_mut(column) {
                *slot += step;
            }
            if let Some(slot) = minus.get_mut(column) {
                *slot -= step;
            }
            let (Ok(plus_design), Ok(minus_design)) = (self.refit(&plus), self.refit(&minus))
            else {
                continue;
            };
            let plus_residual = plus_design.residuals(linear);
            let minus_residual = minus_design.residuals(linear);
            for row in 0..rows {
                if let Some(slot) = jacobian.get_mut(row).and_then(|line| line.get_mut(column)) {
                    *slot = (plus_residual.get(row).copied().unwrap_or(0.0)
                        - minus_residual.get(row).copied().unwrap_or(0.0))
                        / (2.0 * step);
                }
            }
        }
        jacobian
    }
}
