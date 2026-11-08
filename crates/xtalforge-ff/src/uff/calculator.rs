//! The UFF calculator registered as `uff` (Task 4.2.4).
//!
//! Ported from Crystal-Builder `xtal/ff/uff/calculator.py`. This wrapper is
//! deliberately thin: [`UffModel`] owns the typing and topology, and the
//! [`Calculator`] methods build one model per query so the trait remains
//! usable from an optimizer that hands in a new structure each step. Callers
//! that iterate should hold a [`UffModel`] and use its `*_at` methods.

use xtalforge_core::Structure;

use crate::Error;
use crate::calculator::{Calculator, CalculatorOptions};
use crate::uff::model::UffModel;

/// UFF and UFF4MOF over a structure's P1 cell.
#[derive(Debug, Clone, Default)]
pub struct UffCalculator {
    options: CalculatorOptions,
}

impl UffCalculator {
    /// Builds a calculator from options.
    ///
    /// # Errors
    ///
    /// Returns [`Error::InvalidInput`] when the options are not usable.
    pub fn new(options: CalculatorOptions) -> Result<Self, Error> {
        options.validate()?;
        Ok(Self { options })
    }

    /// A calculator with default options (UFF4MOF parameters, van der Waals
    /// on, 12 Angstrom cutoff).
    ///
    /// # Errors
    ///
    /// Returns the errors of [`UffCalculator::new`].
    pub fn with_defaults() -> Result<Self, Error> {
        Self::new(CalculatorOptions::default())
    }

    /// Options this calculator was built with.
    #[must_use]
    pub fn options(&self) -> &CalculatorOptions {
        &self.options
    }

    /// Prepares a model for one structure.
    ///
    /// # Errors
    ///
    /// Returns the errors of [`UffModel::build`].
    pub fn model(&self, structure: &Structure) -> Result<UffModel, Error> {
        UffModel::build(structure, self.options.clone())
    }
}

impl Calculator for UffCalculator {
    fn name(&self) -> &'static str {
        "uff"
    }

    fn energy(&self, structure: &Structure) -> Result<f64, Error> {
        self.model(structure)?.energy()
    }

    fn energy_forces(&self, structure: &Structure) -> Result<(f64, Vec<[f64; 3]>), Error> {
        self.model(structure)?.energy_forces()
    }

    fn provides_forces(&self) -> bool {
        true
    }

    fn provides_stress(&self) -> bool {
        true
    }

    fn stress(&self, structure: &Structure) -> Result<Option<[[f64; 3]; 3]>, Error> {
        Ok(Some(self.model(structure)?.stress()?))
    }
}
