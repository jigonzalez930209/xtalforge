//! Internal model evaluation for [`super::refine_phase`].

use std::f64::consts::{LN_2, PI};

use xtalforge_core::{Error, Site, Structure, StructureBuilder};
use xtalforge_io::xy::PatternPoint;

use crate::pxrd::{PxrdOptions, calculate_pxrd};
use crate::rietveld::RietveldOptions;

use super::fit::{fit_scale_background, weighted_rwp};

/// Baseline `U_iso` in Angstrom^2 for a site without a displacement parameter.
pub(super) const DEFAULT_U_ISO: f64 = 0.005;

/// One refinement parameter.
#[derive(Debug, Clone)]
pub(super) struct Parameter {
    pub(super) name: String,
    pub(super) kind: ParameterKind,
}

/// Kind of a structural parameter.
#[derive(Debug, Clone, Copy)]
pub(super) enum ParameterKind {
    /// One fractional coordinate of one asymmetric site.
    Position { site: usize, axis: usize },
    /// Global multiplier on every site's isotropic ADP.
    AdpScale,
    /// Occupancy of one asymmetric site.
    Occupancy { site: usize },
}

/// One fitted model at fixed structural parameters.
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
}

impl ModelContext<'_> {
    /// Builds the curve for one parameter vector without fitting anything else.
    pub(super) fn curve(&self, values: &[f64]) -> Result<Vec<f64>, Error> {
        let structure = make_structure(self.base, self.base_u, self.parameters, values)?;
        evaluate_curve(&structure, self.observed, self.pxrd)
    }

    /// Builds the full model, fitting scale and background analytically.
    pub(super) fn model(&self, values: &[f64], terms: usize) -> Result<Model, Error> {
        let structure = make_structure(self.base, self.base_u, self.parameters, values)?;
        let curve = evaluate_curve(&structure, self.observed, self.pxrd)?;
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

    /// Central-difference Jacobian of the model curve with respect to the
    /// structural parameters, weighted by the inverse square root of the
    /// observed intensity.
    pub(super) fn jacobian(&self, values: &[f64]) -> Result<Vec<Vec<f64>>, Error> {
        let mut jacobian = vec![vec![0.0; self.parameters.len()]; self.observed.len()];
        for (column, parameter) in self.parameters.iter().enumerate() {
            let (low, high) = bounds(parameter.kind);
            let h = finite_difference(parameter.kind);
            let current = values.get(column).copied().unwrap_or(0.0);
            let plus_value = (current + h).min(high);
            let minus_value = (current - h).max(low);
            let denominator = plus_value - minus_value;
            if denominator.abs() < 1e-12 {
                continue;
            }
            let mut plus = values.to_vec();
            if let Some(slot) = plus.get_mut(column) {
                *slot = plus_value;
            }
            let mut minus = values.to_vec();
            if let Some(slot) = minus.get_mut(column) {
                *slot = minus_value;
            }
            let plus_curve = self.curve(&plus)?;
            let minus_curve = self.curve(&minus)?;
            for (row, point) in self.observed.iter().enumerate() {
                let Some(plus_value) = plus_curve.get(row).copied() else {
                    continue;
                };
                let Some(minus_value) = minus_curve.get(row).copied() else {
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
    if options.refine_positions {
        for (site, entry) in structure.sites().iter().enumerate() {
            for (axis, suffix) in ["x", "y", "z"].into_iter().enumerate() {
                parameters.push(Parameter {
                    name: format!("{}.{suffix}", entry.label()),
                    kind: ParameterKind::Position { site, axis },
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

pub(super) fn initial_value(kind: ParameterKind, sites: &[Site]) -> f64 {
    match kind {
        ParameterKind::Position { site, axis } => sites
            .get(site)
            .and_then(|entry| entry.frac().get(axis).copied())
            .unwrap_or(0.0),
        ParameterKind::AdpScale => 1.0,
        ParameterKind::Occupancy { site } => sites.get(site).map_or(1.0, |entry| entry.occupancy()),
    }
}

pub(super) fn bounds(kind: ParameterKind) -> (f64, f64) {
    match kind {
        ParameterKind::Position { .. } => (-1.0, 2.0),
        ParameterKind::AdpScale => (0.05, 20.0),
        ParameterKind::Occupancy { .. } => (0.0, 1.0),
    }
}

pub(super) fn max_step(kind: ParameterKind) -> f64 {
    match kind {
        ParameterKind::Position { .. } => 0.05,
        ParameterKind::AdpScale => 0.25,
        ParameterKind::Occupancy { .. } => 0.10,
    }
}

fn finite_difference(kind: ParameterKind) -> f64 {
    match kind {
        ParameterKind::Position { .. } => 0.002,
        ParameterKind::AdpScale => 0.02,
        ParameterKind::Occupancy { .. } => 0.005,
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
    for (parameter, value) in parameters.iter().zip(values.iter()) {
        match parameter.kind {
            ParameterKind::Position { site, axis } => {
                let Some(current) = sites.get(site).cloned() else {
                    continue;
                };
                let mut frac = current.frac();
                if let Some(slot) = frac.get_mut(axis) {
                    *slot = *value;
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
        }
    }
    if let Some(scale) = adp_scale {
        for (index, site) in sites.iter_mut().enumerate() {
            let base_u = base_u.get(index).copied().unwrap_or(DEFAULT_U_ISO);
            *site = site.clone().with_u_iso((base_u * scale).max(0.0))?;
        }
    }
    rebuild(base, sites)
}

fn rebuild(base: &Structure, sites: Vec<Site>) -> Result<Structure, Error> {
    let mut builder = StructureBuilder::new(base.name(), *base.lattice()).hall(base.hall());
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

fn evaluate_curve(
    structure: &Structure,
    observed: &[PatternPoint],
    pxrd: &PxrdOptions,
) -> Result<Vec<f64>, Error> {
    let pattern = calculate_pxrd(structure, pxrd)?;
    if pattern.reflections.is_empty() {
        return Err(Error::InvalidInput(
            "the structure has no reflections in the observed range".to_owned(),
        ));
    }
    // Same area-normalised Gaussian convention as `pxrd::profile::sample_pattern`,
    // so a synthetic pattern built from `calculate_pxrd` fits with a scale of 1.
    let normalisation = 2.0 / pxrd.fwhm * (LN_2 / PI).sqrt();
    let mut curve = vec![0.0; observed.len()];
    for (slot, point) in curve.iter_mut().zip(observed.iter()) {
        let mut value = 0.0;
        for reflection in &pattern.reflections {
            let delta = (point.two_theta - reflection.two_theta) / pxrd.fwhm;
            value += reflection.intensity * normalisation * (-4.0 * LN_2 * delta * delta).exp();
        }
        *slot = value;
    }
    Ok(curve)
}
