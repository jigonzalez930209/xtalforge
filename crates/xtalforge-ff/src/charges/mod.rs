//! Atomic charges and periodic electrostatics (Tasks 4.3.1 to 4.3.5).
//!
//! Four entry points share one [`Charges`] result and one periodic-image
//! convention:
//!
//! * [`qeq_charges`] — Rappé-Goddard charge equilibration with the UFF
//!   electronegativity and hardness table (Task 4.3.1).
//! * [`eqeq_charges`] — EQeq, the same solver with each element expanded
//!   about its usual oxidation state and measured ionisation energies
//!   (Task 4.3.2).
//! * [`eeq_charges`] — EEQ-style equilibration with an Ohno-Klopman kernel
//!   and an explicit Coulomb radius per element (Task 4.3.3).
//! * [`ewald_energy_forces`] — converged periodic Coulomb energy and forces
//!   for a set of point charges (Task 4.3.5).
//!
//! Charges are in units of the elementary charge and energies are eV.  Every
//! charge routine works on the P1 cell: `per_atom` has one entry per
//! [`xtalforge_core::P1Cell`] atom, and dummy `X` atoms are excluded from the
//! solve and carry zero.  [`Charges::per_site`] and [`Charges::apply_to`]
//! turn that vector into the per-asymmetric-site values a structure stores
//! (Task 4.3.4).

use std::collections::BTreeMap;

use serde::{Deserialize, Serialize};
use xtalforge_core::Structure;

use crate::Error;
use crate::structure_util::rebuild;

pub(crate) mod data;
pub(crate) mod eeq;
pub(crate) mod eqeq;
pub mod ewald;
pub(crate) mod periodic;
pub(crate) mod qeq;

pub use ewald::{EwaldOptions, EwaldResult, ewald_energy, ewald_energy_forces};

/// `e^2 / (4 pi eps0)` in eV Angstrom.
///
/// This is the unit that makes electronegativities and hardnesses in eV and
/// distances in Angstrom meet in one Coulomb expression.
pub const COULOMB_EV_ANGSTROM: f64 = 14.399_645;

/// Largest cell the dense charge-equilibration solver accepts.
///
/// The method solves an `(n + 1) x (n + 1)` system and builds an `n x n`
/// interaction matrix, so the cost is quadratic in memory and cubic in time.
/// Refusing a cell this large is better than appearing to hang.
pub const MAX_CHARGE_ATOMS: usize = 2000;

/// Charges for the atoms of a P1 cell.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct Charges {
    /// Charge of each P1 atom, in units of the elementary charge.
    pub per_atom: Vec<f64>,
    /// Sum of `per_atom`, in the same units.
    pub total: f64,
    /// Number of solver passes; the direct dense solve reports one.
    pub iterations: usize,
}

impl Charges {
    /// An empty charge set.
    #[must_use]
    pub fn empty() -> Self {
        Self {
            per_atom: Vec::new(),
            total: 0.0,
            iterations: 0,
        }
    }

    /// Number of charges.
    #[must_use]
    pub fn len(&self) -> usize {
        self.per_atom.len()
    }

    /// Whether there are no charges.
    #[must_use]
    pub fn is_empty(&self) -> bool {
        self.per_atom.is_empty()
    }

    /// Charge of each asymmetric site: the mean over its P1 images.
    ///
    /// A site's images are one atom as far as the crystal is concerned, so a
    /// symmetric cell gives every image the same charge to rounding and this
    /// mean is the value the structure should store.  A site with no image —
    /// which a valid P1 expansion cannot produce — gets zero.
    ///
    /// # Errors
    ///
    /// Returns [`Error::InvalidInput`] when the length does not match the
    /// structure's P1 atom count, or a core error when the expansion fails.
    pub fn per_site(&self, structure: &Structure) -> Result<Vec<f64>, Error> {
        let cell = structure.p1()?;
        if self.per_atom.len() != cell.n_atoms() {
            return Err(Error::InvalidInput(format!(
                "got {} charges for {} P1 atoms",
                self.per_atom.len(),
                cell.n_atoms()
            )));
        }
        let mut sums = vec![0.0; structure.sites().len()];
        let mut counts = vec![0usize; structure.sites().len()];
        for (atom, charge) in cell.atoms.iter().zip(&self.per_atom) {
            if let (Some(sum), Some(count)) =
                (sums.get_mut(atom.site_idx), counts.get_mut(atom.site_idx))
            {
                *sum += charge;
                *count += 1;
            }
        }
        Ok(sums
            .into_iter()
            .zip(counts)
            .map(|(sum, count)| if count == 0 { 0.0 } else { sum / count as f64 })
            .collect())
    }

    /// Copy of `structure` with each site's charge set from this result.
    ///
    /// The structure keeps its lattice, space group, bonds and metadata; only
    /// `Site::with_charge` changes (Task 4.3.4).  Writing the result to CIF
    /// `_atom_site_charge` or a LAMMPS data file is then a format concern, not
    /// a force-field one.
    ///
    /// # Errors
    ///
    /// Returns the errors of [`Charges::per_site`] or of rebuilding the
    /// structure.
    pub fn apply_to(&self, structure: &Structure) -> Result<Structure, Error> {
        let means = self.per_site(structure)?;
        let sites = structure
            .sites()
            .iter()
            .zip(means)
            .map(|(site, charge)| site.clone().with_charge(charge))
            .collect::<Result<Vec<_>, _>>()?;
        rebuild(structure, *structure.lattice(), sites)
    }

    /// Total charge from the per-atom values, without trusting the stored
    /// `total`.
    #[must_use]
    pub fn recompute_total(&self) -> f64 {
        self.per_atom.iter().sum()
    }
}

/// Common options for charge equilibration.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct ChargeOptions {
    /// Total charge the cell must carry, in units of the elementary charge.
    pub total_charge: f64,
    /// Real-space cutoff in Angstrom for the periodic interaction sums.
    pub cutoff: f64,
    /// Relative accuracy of the Ewald split used by EQeq and Ewald.
    pub accuracy: f64,
    /// Dielectric screening of the Coulomb interaction.
    ///
    /// QEq ignores it (its shielding is in the Ohno-Klopman kernel); EQeq
    /// uses it and its reference default is 1.2.
    pub dielectric: f64,
    /// Per-element charge centres for EQeq, overriding the table.
    ///
    /// An empty map uses the common oxidation state in the built-in table
    /// (the reference's `CHARGE_CENTRES`).
    pub centres: BTreeMap<String, i32>,
}

impl Default for ChargeOptions {
    fn default() -> Self {
        Self {
            total_charge: 0.0,
            cutoff: 12.0,
            accuracy: 1e-6,
            dielectric: 1.0,
            centres: BTreeMap::new(),
        }
    }
}

impl ChargeOptions {
    /// Options for the published EQeq defaults: dielectric 1.2.
    #[must_use]
    pub fn eqeq() -> Self {
        Self {
            dielectric: 1.2,
            ..Self::default()
        }
    }

    /// Creates options with a non-zero total charge.
    ///
    /// # Errors
    ///
    /// Returns [`Error::InvalidInput`] when the charge is not finite.
    pub fn with_total_charge(mut self, total_charge: f64) -> Result<Self, Error> {
        if !total_charge.is_finite() {
            return Err(Error::InvalidInput(
                "total charge must be finite".to_owned(),
            ));
        }
        self.total_charge = total_charge;
        Ok(self)
    }

    /// Checks that the options are usable.
    ///
    /// # Errors
    ///
    /// Returns [`Error::InvalidInput`] for a non-finite total charge, a
    /// non-positive cutoff or dielectric, or an accuracy outside `(0, 1)`.
    pub fn validate(&self) -> Result<(), Error> {
        if !self.total_charge.is_finite() {
            return Err(Error::InvalidInput(
                "total charge must be finite".to_owned(),
            ));
        }
        if !self.cutoff.is_finite() || self.cutoff <= 0.0 {
            return Err(Error::InvalidInput(format!(
                "charge cutoff must be positive and finite, got {}",
                self.cutoff
            )));
        }
        if !(self.accuracy.is_finite() && 0.0 < self.accuracy && self.accuracy < 1.0) {
            return Err(Error::InvalidInput(format!(
                "Ewald accuracy must be in (0, 1), got {}",
                self.accuracy
            )));
        }
        if !self.dielectric.is_finite() || self.dielectric <= 0.0 {
            return Err(Error::InvalidInput(format!(
                "dielectric must be positive and finite, got {}",
                self.dielectric
            )));
        }
        Ok(())
    }
}

/// QEq charges for a structure (Task 4.3.1).
///
/// Every atom is given the same electronegativity, subject to the cell
/// carrying `options.total_charge`.  The interaction is the Ohno-Klopman
/// shielded Coulomb form summed over periodic images inside the cutoff, as in
/// Crystal-Builder `xtal/ff/uff/qeq.py`; it is not the published Slater
/// integral, and the charges are an estimate to review rather than a
/// published result.
///
/// # Errors
///
/// Returns [`Error::InvalidInput`] for bad options or a singular system,
/// [`Error::UnknownElement`] for an element without QEq parameters, and a
/// core error when the P1 expansion fails.
pub fn qeq_charges(structure: &Structure, options: &ChargeOptions) -> Result<Charges, Error> {
    options.validate()?;
    qeq::solve(structure, options)
}

/// EQeq charges for a structure (Task 4.3.2).
///
/// The solver is QEq's; only the parameters change.  Each element is expanded
/// about its usual charge centre and its `chi` and hardness come from measured
/// ionisation energies, which is the correction that keeps a framework's
/// metal from running away from the neutral-atom parabola.
///
/// The periodic Coulomb matrix is the converged Ewald one plus the reference's
/// Gaussian overlap correction; the reference sums its published charges with
/// a dielectric of 1.2, available as [`ChargeOptions::eqeq`].
///
/// # Errors
///
/// Returns [`Error::InvalidInput`] for bad options or a singular system,
/// [`Error::UnknownElement`] for an element outside the EQeq table, and a core
/// error when the P1 expansion fails.
pub fn eqeq_charges(structure: &Structure, options: &ChargeOptions) -> Result<Charges, Error> {
    options.validate()?;
    eqeq::solve(structure, options)
}

/// EEQ-style charges for a structure (Task 4.3.3).
///
/// The diagonal is each element's chemical hardness and the off-diagonal is
/// the Ohno-Klopman kernel with an explicit Coulomb radius per element; the
/// total charge enters through one Lagrange multiplier and every periodic
/// image inside `options.cutoff` is summed.
///
/// The module documents the approximation precisely: sci-form's fitted EEQ2019
/// parameter table is not vendored here, so `chi` and `J` come from the local
/// QEq table and the radius is the Ohno radius implied by that hardness.  The
/// default dielectric of one therefore reproduces [`qeq_charges`]; use
/// [`ChargeOptions::eqeq`] (dielectric 1.2) for the screened EEQ-style result.
///
/// # Errors
///
/// Returns [`Error::InvalidInput`] for bad options or a singular system,
/// [`Error::UnknownElement`] for an element outside the local table, and a
/// core error when the P1 expansion fails.
pub fn eeq_charges(structure: &Structure, options: &ChargeOptions) -> Result<Charges, Error> {
    options.validate()?;
    eeq::solve(structure, options)
}

/// Applies charges to a structure's sites.
///
/// A convenience wrapper for [`Charges::apply_to`].
///
/// # Errors
///
/// Returns the errors of [`Charges::apply_to`].
pub fn apply_charges(structure: &Structure, charges: &Charges) -> Result<Structure, Error> {
    charges.apply_to(structure)
}
