//! Internal model evaluation for [`super::refine_phase`].

use xtalforge_core::{Error, Lattice, Site, Structure, StructureBuilder};
use xtalforge_io::xy::PatternPoint;

use crate::pxrd::PxrdOptions;
use crate::rietveld::RietveldOptions;

use super::curve::{CurveOptions, evaluate_curve};
use super::fit::{fit_scale_background, weighted_rwp};
use super::wyckoff::allowed_directions;

/// Baseline `U_iso` in Angstrom^2 for a site without a displacement parameter.
pub(super) const DEFAULT_U_ISO: f64 = 0.005;

/// One refinement parameter.
#[derive(Debug, Clone)]
pub(super) struct Parameter {
    pub(super) name: String,
    pub(super) kind: ParameterKind,
}

/// Kind of a refinement parameter.
#[derive(Debug, Clone, Copy)]
pub(super) enum ParameterKind {
    /// Displacement of one asymmetric site along a Wyckoff-allowed direction.
    Position { site: usize, direction: [f64; 3] },
    /// Global multiplier on every site's isotropic ADP.
    AdpScale,
    /// Occupancy of one asymmetric site.
    Occupancy { site: usize },
    /// Zero shift in degrees two-theta.
    Zero,
    /// One free cell number.
    Cell { index: usize },
    /// Multiplier on the profile FWHM.
    WidthScale,
    /// March-Dollase preferred-orientation coefficient.
    Texture,
}

/// One fitted model at fixed parameters.
pub(super) struct Model {
    pub(super) structure: Structure,
    pub(super) scale: f64,
    pub(super) background: Vec<f64>,
    pub(super) model: Vec<f64>,
    pub(super) rwp: f64,
}

/// Everything the model evaluator needs, borrowed once per refinement.
pub(super) struct ModelContext<'a> {
    pub(super) base: &'a Structure,
    pub(super) base_u: &'a [f64],
    pub(super) parameters: &'a [Parameter],
    pub(super) observed: &'a [PatternPoint],
    pub(super) pxrd: &'a PxrdOptions,
    pub(super) texture_axis: [f64; 3],
    pub(super) doublet: bool,
}

impl ModelContext<'_> {
    fn curve_options(&self, values: &[f64]) -> CurveOptions {
        let mut options = CurveOptions {
            texture_axis: self.texture_axis,
            doublet: self.doublet,
            ..CurveOptions::default()
        };
        for (parameter, value) in self.parameters.iter().zip(values.iter()) {
            match parameter.kind {
                ParameterKind::Zero => options.zero = *value,
                ParameterKind::WidthScale => options.width_scale = *value,
                ParameterKind::Texture => options.texture = *value,
                _ => {}
            }
        }
        options
    }

    /// Builds the structure for one parameter vector.
    pub(super) fn structure(&self, values: &[f64]) -> Result<Structure, Error> {
        make_structure(self.base, self.base_u, self.parameters, values)
    }

    /// Builds the full model, fitting scale and background analytically.
    pub(super) fn model(&self, values: &[f64], terms: usize) -> Result<Model, Error> {
        let structure = self.structure(values)?;
        let options = self.curve_options(values);
        let curve = evaluate_curve(&structure, self.observed, self.pxrd, &options)?;
        let (scale, background, model) = fit_scale_background(&curve, self.observed, terms)?;
        let rwp = weighted_rwp(self.observed, &model)?;
        Ok(Model {
            structure,
            scale,
            background,
            model,
            rwp,
        })
    }

    /// Central-difference Jacobian of the model curve, weighted by the inverse
    /// square root of the observed intensity.
    pub(super) fn jacobian(
        &self,
        values: &[f64],
        indices: &[usize],
        terms: usize,
    ) -> Result<Vec<Vec<f64>>, Error> {
        let mut jacobian = vec![vec![0.0; indices.len()]; self.observed.len()];
        for (column, index) in indices.iter().enumerate() {
            let Some(parameter) = self.parameters.get(*index) else {
                continue;
            };
            let (low, high) = bounds(parameter.kind);
            let h = finite_difference(parameter.kind);
            let current = values.get(*index).copied().unwrap_or(0.0);
            let plus_value = (current + h).min(high);
            let minus_value = (current - h).max(low);
            let denominator = plus_value - minus_value;
            if denominator.abs() < 1e-12 {
                continue;
            }
            let mut plus = values.to_vec();
            if let Some(slot) = plus.get_mut(*index) {
                *slot = plus_value;
            }
            let mut minus = values.to_vec();
            if let Some(slot) = minus.get_mut(*index) {
                *slot = minus_value;
            }
            // The full model, not the raw curve: the analytic scale and
            // background refit after every perturbation, so the derivative
            // already includes their projection onto the structural change.
            let plus_model = self.model(&plus, terms)?;
            let minus_model = self.model(&minus, terms)?;
            for (row, point) in self.observed.iter().enumerate() {
                let Some(plus_value) = plus_model.model.get(row).copied() else {
                    continue;
                };
                let Some(minus_value) = minus_model.model.get(row).copied() else {
                    continue;
                };
                let root = 1.0 / point.intensity.abs().max(1.0).sqrt();
                let derivative = (plus_value - minus_value) / denominator;
                if let Some(slot) = jacobian.get_mut(row).and_then(|line| line.get_mut(column)) {
                    *slot = root * derivative;
                }
            }
        }
        Ok(jacobian)
    }
}

pub(super) fn pxrd_options(options: &RietveldOptions, observed: &[PatternPoint]) -> PxrdOptions {
    let first = observed
        .first()
        .map_or(options.two_theta_range[0], |point| point.two_theta);
    let last = observed
        .last()
        .map_or(options.two_theta_range[1], |point| point.two_theta);
    PxrdOptions {
        radiation: options.radiation,
        two_theta_range: [
            options.two_theta_range[0].min(first).max(0.0),
            options.two_theta_range[1].max(last).min(180.0),
        ],
        step: options.step,
        fwhm: options.fwhm,
    }
}

pub(super) fn build_parameters(
    structure: &Structure,
    options: &super::RefineOptions,
) -> Vec<Parameter> {
    let mut parameters = Vec::new();
    if options.refine_zero {
        parameters.push(Parameter {
            name: "zero".to_owned(),
            kind: ParameterKind::Zero,
        });
    }
    if options.refine_cell {
        for index in free_cell_indices(structure.lattice().parameters()) {
            parameters.push(Parameter {
                name: format!("cell.{index}"),
                kind: ParameterKind::Cell { index },
            });
        }
    }
    if options.refine_width {
        parameters.push(Parameter {
            name: "width_scale".to_owned(),
            kind: ParameterKind::WidthScale,
        });
    }
    if options.refine_texture {
        parameters.push(Parameter {
            name: "texture".to_owned(),
            kind: ParameterKind::Texture,
        });
    }
    if options.refine_positions {
        for (site, entry) in structure.sites().iter().enumerate() {
            let directions = allowed_directions(structure.operations(), entry.frac());
            for (axis, direction) in directions.iter().enumerate() {
                parameters.push(Parameter {
                    name: format!("{}.d{axis}", entry.label()),
                    kind: ParameterKind::Position {
                        site,
                        direction: *direction,
                    },
                });
            }
        }
    }
    if options.refine_adps {
        parameters.push(Parameter {
            name: "adp_scale".to_owned(),
            kind: ParameterKind::AdpScale,
        });
    }
    if options.refine_occupancies {
        for (site, entry) in structure.sites().iter().enumerate() {
            parameters.push(Parameter {
                name: format!("{}.occupancy", entry.label()),
                kind: ParameterKind::Occupancy { site },
            });
        }
    }
    parameters
}

/// Cell numbers the metric family leaves free.
/// Restores the equal-length ties of the metric family after one free number
/// changed.
fn enforce_cell_ties(cell: &mut [f64; 6], base: [f64; 6]) {
    let close = |first: f64, second: f64| (first - second).abs() <= 0.01;
    let right = |angle: f64| (angle - 90.0).abs() <= 0.01;
    let [a, b, c, alpha, beta, gamma] = base;
    if close(a, b) && close(b, c) && right(alpha) && right(beta) && right(gamma) {
        cell[1] = cell[0];
        cell[2] = cell[0];
    } else if close(a, b) && right(alpha) && right(beta) && (right(gamma) || close(gamma, 120.0)) {
        cell[1] = cell[0];
    }
    cell[3] = base[3];
    cell[4] = base[4];
    cell[5] = base[5];
}

pub(super) fn free_cell_indices(parameters: [f64; 6]) -> Vec<usize> {
    let [a, b, c, alpha, beta, gamma] = parameters;
    let close = |first: f64, second: f64| (first - second).abs() <= 0.01;
    let right = |angle: f64| (angle - 90.0).abs() <= 0.01;
    if close(a, b) && close(b, c) && right(alpha) && right(beta) && right(gamma) {
        vec![0]
    } else if close(a, b) && right(alpha) && right(beta) && (right(gamma) || close(gamma, 120.0)) {
        vec![0, 2]
    } else if right(alpha) && right(beta) && right(gamma) {
        vec![0, 1, 2]
    } else if right(alpha) && right(gamma) {
        vec![0, 1, 2, 4]
    } else {
        vec![0, 1, 2, 3, 4, 5]
    }
}

pub(super) fn initial_value(kind: ParameterKind, structure: &Structure) -> f64 {
    match kind {
        ParameterKind::Position { .. } | ParameterKind::Zero => 0.0,
        ParameterKind::AdpScale | ParameterKind::WidthScale | ParameterKind::Texture => 1.0,
        ParameterKind::Occupancy { site } => structure
            .sites()
            .get(site)
            .map_or(1.0, |entry| entry.occupancy()),
        ParameterKind::Cell { index } => structure
            .lattice()
            .parameters()
            .get(index)
            .copied()
            .unwrap_or(1.0),
    }
}

pub(super) fn bounds(kind: ParameterKind) -> (f64, f64) {
    match kind {
        ParameterKind::Position { .. } => (-0.5, 0.5),
        ParameterKind::AdpScale => (0.05, 20.0),
        ParameterKind::Occupancy { .. } => (0.0, 1.0),
        ParameterKind::Zero => (-2.0, 2.0),
        ParameterKind::Cell { index } if index < 3 => (0.5, 200.0),
        ParameterKind::Cell { .. } => (30.0, 150.0),
        ParameterKind::WidthScale => (0.1, 10.0),
        ParameterKind::Texture => (0.1, 10.0),
    }
}

fn finite_difference(kind: ParameterKind) -> f64 {
    match kind {
        ParameterKind::Position { .. } => 0.002,
        ParameterKind::AdpScale => 0.02,
        ParameterKind::Occupancy { .. } => 0.005,
        ParameterKind::Zero => 0.005,
        ParameterKind::Cell { index } if index < 3 => 0.01,
        ParameterKind::Cell { .. } => 0.02,
        ParameterKind::WidthScale | ParameterKind::Texture => 0.01,
    }
}

pub(super) fn clamp(value: f64, low: f64, high: f64) -> f64 {
    if value.is_finite() {
        value.max(low).min(high)
    } else {
        low
    }
}

fn make_structure(
    base: &Structure,
    base_u: &[f64],
    parameters: &[Parameter],
    values: &[f64],
) -> Result<Structure, Error> {
    let mut sites = base.sites().to_vec();
    let mut adp_scale: Option<f64> = None;
    let mut cell = base.lattice().parameters();
    let mut cell_changed = false;
    for (parameter, value) in parameters.iter().zip(values.iter()) {
        match parameter.kind {
            ParameterKind::Position { site, direction } => {
                let Some(current) = sites.get(site).cloned() else {
                    continue;
                };
                let mut frac = current.frac();
                for (axis, component) in direction.iter().enumerate() {
                    if let Some(slot) = frac.get_mut(axis) {
                        *slot = (*slot + component * value).rem_euclid(1.0);
                    }
                }
                if let Some(target) = sites.get_mut(site) {
                    *target = current.with_frac(frac)?;
                }
            }
            ParameterKind::AdpScale => adp_scale = Some(*value),
            ParameterKind::Occupancy { site } => {
                let Some(current) = sites.get(site).cloned() else {
                    continue;
                };
                if let Some(target) = sites.get_mut(site) {
                    *target = current.with_occupancy(value.clamp(0.0, 1.0))?;
                }
            }
            ParameterKind::Cell { index } => {
                if let Some(slot) = cell.get_mut(index) {
                    *slot = *value;
                    cell_changed = true;
                }
            }
            ParameterKind::Zero | ParameterKind::WidthScale | ParameterKind::Texture => {}
        }
    }
    if let Some(scale) = adp_scale {
        for (index, site) in sites.iter_mut().enumerate() {
            let base_u = base_u.get(index).copied().unwrap_or(DEFAULT_U_ISO);
            if let Some(aniso) = site.u_aniso() {
                // Anisotropic parameters take precedence in the scattering
                // engine, so the global multiplier scales all six components
                // and leaves `U_iso` alone.
                *site = site
                    .clone()
                    .with_u_aniso(aniso.map(|value| value * scale))?;
            } else {
                *site = site.clone().with_u_iso((base_u * scale).max(0.0))?;
            }
        }
    }
    if cell_changed {
        // The free cell numbers move one at a time; the metric family ties
        // must be restored after every move or a cubic a would no longer
        // equal b.
        enforce_cell_ties(&mut cell, base.lattice().parameters());
    }
    let lattice = if cell_changed {
        Lattice::from_parameters(cell[0], cell[1], cell[2], cell[3], cell[4], cell[5])?
    } else {
        *base.lattice()
    };
    rebuild(base, sites, lattice)
}

fn rebuild(base: &Structure, sites: Vec<Site>, lattice: Lattice) -> Result<Structure, Error> {
    let mut builder = StructureBuilder::new(base.name(), lattice).hall(base.hall());
    for site in sites {
        builder = builder.site(site);
    }
    for operation in base.operations() {
        builder = builder.operation(*operation);
    }
    for bond in base.explicit_bonds() {
        builder = builder.explicit_bond(*bond);
    }
    for bond in base.topology_bonds() {
        builder = builder.topology_bond(*bond);
    }
    for (key, value) in base.metadata() {
        builder = builder.metadata(key.clone(), value.clone());
    }
    builder.build()
}
