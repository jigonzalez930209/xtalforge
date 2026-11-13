//! A calculator wrapper that adds scan penalties (Task 4.5.2).

use xtalforge_core::Structure;

use crate::Error;
use crate::calculator::Calculator;
use crate::math::sub;
use crate::scan::constraints::Constraint;

/// A calculator with a quadratic penalty added for each constraint.
pub(crate) struct PenalizedCalculator<'a> {
    inner: &'a dyn Calculator,
    constraints: Vec<Constraint>,
    constant: f64,
}

impl<'a> PenalizedCalculator<'a> {
    /// Wraps `inner` with `constraints` at penalty constant `constant`.
    pub(crate) fn new(
        inner: &'a dyn Calculator,
        constraints: Vec<Constraint>,
        constant: f64,
    ) -> Self {
        Self {
            inner,
            constraints,
            constant,
        }
    }
}

impl Calculator for PenalizedCalculator<'_> {
    fn name(&self) -> &'static str {
        "penalized"
    }

    fn energy(&self, structure: &Structure) -> Result<f64, Error> {
        let mut energy = self.inner.energy(structure)?;
        let positions = structure.cart_sites();
        for constraint in &self.constraints {
            let (penalty, _) =
                constraint.penalty(&positions, structure.lattice(), self.constant)?;
            energy += penalty;
        }
        Ok(energy)
    }

    fn energy_forces(&self, structure: &Structure) -> Result<(f64, Vec<[f64; 3]>), Error> {
        let (mut energy, mut forces) = self.inner.energy_forces(structure)?;
        let positions = structure.cart_sites();
        for constraint in &self.constraints {
            let (penalty, gradient) =
                constraint.penalty(&positions, structure.lattice(), self.constant)?;
            energy += penalty;
            for (force, component) in forces.iter_mut().zip(&gradient) {
                *force = sub(*force, *component);
            }
        }
        Ok((energy, forces))
    }

    fn provides_forces(&self) -> bool {
        true
    }

    fn provides_stress(&self) -> bool {
        false
    }
}
