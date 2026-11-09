//! Prepare a structure for simulation (Task 4.7).
//!
//! A refined structure describes the average over all cells.  This module runs
//! the force-field part of the preparation pipeline as an explicit, ordered
//! sequence and records every decision in a [`PrepareReport`]:
//!
//! 1. relax the geometry with FIRE over the variables [`OptimizeOptions`]
//!    describes (optional);
//! 2. assign QEq charges to the P1 cell and store them per asymmetric site
//!    (optional);
//! 3. evaluate the periodic Coulomb energy of those P1 charges with Ewald
//!    summation (optional);
//! 4. build the Gamma-point Hessian and its normal modes (optional).
//!
//! The steps run in that fixed order, so the same input and options give the
//! same report on every platform.  A disabled step is still listed with
//! `done = false` and a `skipped` detail, which makes the report explain the
//! pipeline rather than only its products.  An error from any step is returned
//! unchanged; a failed preparation never reports a partial success.
//!
//! **Approximations.**  The energy in the report is the calculator's energy on
//! the prepared P1 geometry and does not include the Ewald term, which is
//! reported separately.  Charges are QEq, not a function of the calculator, so
//! applying them does not change a UFF energy.  The vibration step uses the
//! default crystal projection ([`HessianOptions::default`]); a caller that
//! wants the rigid rotations of an isolated molecule removed must run
//! [`crate::vibrations::vibrational_modes`] itself with
//! [`HessianOptions::molecule`](crate::vibrations::HessianOptions::molecule).
//! Disorder, partial occupancy, solvent removal, hydrogen addition and charge
//! balance decisions from `core/prepare.py` are not part of this core task and
//! are neither guessed nor hidden.

mod cleanup;

use serde::{Deserialize, Serialize};
use xtalforge_core::{P1Cell, Site, Structure, StructureBuilder};

use crate::Error;
use crate::calculator::Calculator;
use crate::charges::{ChargeOptions, Charges, EwaldOptions, ewald_energy, qeq_charges};
use crate::optimize::{OptimizationResult, OptimizeOptions, fire};
use crate::vibrations::{HessianOptions, hessian, vibrational_modes};

pub use cleanup::{ChargeBalance, CloseContact};

/// Options for one preparation run.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct PrepareOptions {
    /// Whether duplicate sites at the same position are merged.
    pub merge_duplicates: bool,
    /// Whether one component of every disorder group is chosen.
    pub order_disorder: bool,
    /// Whether known solvent molecules are removed from the pores.
    pub remove_solvent: bool,
    /// Whether missing hydrogens are added from the perceived valence.
    pub hydrogens: bool,
    /// Whether close contacts are listed.
    pub check_close_contacts: bool,
    /// Whether the final cell charge is checked.
    pub check_charge_balance: bool,
    /// Whether the geometry is relaxed with FIRE.
    pub minimize: bool,
    /// Whether QEq charges are computed and stored on the sites.
    pub charges: bool,
    /// Whether the periodic Coulomb energy of the P1 charges is evaluated.
    pub ewald: bool,
    /// Whether the Gamma-point Hessian and its frequencies are computed.
    pub vibrations: bool,
    /// Minimizer settings, used when `minimize` is on.
    pub optimize: OptimizeOptions,
    /// Charge settings, used when `charges` is on.
    pub charge_options: ChargeOptions,
    /// Ewald settings, used when `ewald` is on.
    pub ewald_options: EwaldOptions,
}

impl Default for PrepareOptions {
    fn default() -> Self {
        Self {
            merge_duplicates: false,
            order_disorder: false,
            remove_solvent: false,
            hydrogens: false,
            check_close_contacts: false,
            check_charge_balance: false,
            minimize: true,
            charges: true,
            ewald: false,
            vibrations: false,
            optimize: OptimizeOptions::default(),
            charge_options: ChargeOptions::default(),
            ewald_options: EwaldOptions::default(),
        }
    }
}

impl PrepareOptions {
    /// Checks the option groups the enabled steps will use.
    ///
    /// # Errors
    ///
    /// Returns [`Error::InvalidInput`] for the first invalid enabled group, in
    /// pipeline order.
    pub fn validate(&self) -> Result<(), Error> {
        if self.minimize {
            self.optimize.validate()?;
        }
        if self.charges {
            self.charge_options.validate()?;
        }
        if self.ewald {
            self.ewald_options.validate()?;
        }
        Ok(())
    }
}

/// One recorded preparation step.
#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct PrepareStep {
    /// Stable step name: `minimize`, `charges`, `ewald` or `vibrations`.
    pub name: String,
    /// Whether the step ran.
    pub done: bool,
    /// Human-readable detail, including `skipped` for a disabled step.
    pub detail: String,
}

impl PrepareStep {
    fn new(name: &str, done: bool, detail: String) -> Self {
        Self {
            name: name.to_owned(),
            done,
            detail,
        }
    }
}

/// What a preparation run produced.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct PrepareReport {
    /// Prepared structure after every selected step.
    pub structure: Structure,
    /// Every step in pipeline order, including skipped ones.
    pub steps: Vec<PrepareStep>,
    /// Final calculator energy in eV on the prepared geometry.
    pub energy: f64,
    /// Total charge of the P1 cell.
    ///
    /// When the charge step ran this is the solver's constrained total;
    /// otherwise it is the sum of the sites' stored charges.
    pub total_charge: f64,
    /// Ewald energy in eV, present only when the Ewald step ran.
    pub ewald_energy: Option<f64>,
    /// Gamma-point frequencies in cm^-1, empty unless the vibration step ran.
    pub frequencies: Vec<f64>,
    /// Full minimizer result, present only when the minimize step ran.
    pub optimized: Option<OptimizationResult>,
    /// Sites removed by the duplicate, disorder or solvent steps.
    pub removed_sites: usize,
    /// Close contacts, empty unless that step ran.
    pub close_contacts: Vec<CloseContact>,
    /// Final charge balance, present only when that step ran.
    pub charge_balance: Option<ChargeBalance>,
}

/// Runs the preparation pipeline on `structure`.
///
/// The input is not modified.  Each enabled step sees the output of the
/// previous one, and the report keeps the steps, the final energy, the period
/// charges and the optional products.
///
/// # Errors
///
/// Returns [`Error::InvalidInput`] for invalid options and any error from the
/// calculator, the charge solver, Ewald summation or the vibration analysis.
/// Errors are returned at the step that produced them.
pub fn prepare(
    calculator: &dyn Calculator,
    structure: &Structure,
    options: &PrepareOptions,
) -> Result<PrepareReport, Error> {
    options.validate()?;
    let mut working = structure.clone();
    let mut steps = Vec::with_capacity(8);
    let mut optimized = None;
    let mut removed_sites = 0;
    let mut computed_charges: Option<Charges> = None;
    let mut ewald_value = None;
    let mut frequencies = Vec::new();
    let energy;

    if options.merge_duplicates {
        let (next, removed) = cleanup::merge_duplicates(&working)?;
        removed_sites += removed;
        steps.push(PrepareStep::new(
            "duplicates",
            true,
            format!("merged {removed} duplicate site(s)"),
        ));
        working = next;
    } else {
        steps.push(PrepareStep::new(
            "duplicates",
            false,
            "skipped: duplicate sites left as given".to_owned(),
        ));
    }
    if options.order_disorder {
        let (next, removed) = cleanup::order_disorder(&working)?;
        removed_sites += removed;
        steps.push(PrepareStep::new(
            "disorder",
            true,
            format!("chose the disorder components, removed {removed} site(s)"),
        ));
        working = next;
    } else {
        steps.push(PrepareStep::new(
            "disorder",
            false,
            "skipped: partial occupancy left as given".to_owned(),
        ));
    }
    if options.remove_solvent {
        let (next, removed) = cleanup::remove_solvent(&working)?;
        removed_sites += removed;
        steps.push(PrepareStep::new(
            "solvent",
            true,
            format!("removed {removed} solvent molecule(s)"),
        ));
        working = next;
    } else {
        steps.push(PrepareStep::new(
            "solvent",
            false,
            "skipped: solvent left in the pores".to_owned(),
        ));
    }

    if options.hydrogens {
        let filled =
            xtalforge_build::add_hydrogens(&working, &xtalforge_build::HydrogenOptions::default())?;
        let added = filled.sites().len().saturating_sub(working.sites().len());
        steps.push(PrepareStep::new(
            "hydrogens",
            true,
            format!("added {added} hydrogen site(s)"),
        ));
        working = filled;
    } else {
        steps.push(PrepareStep::new(
            "hydrogens",
            false,
            "skipped: hydrogen sites left as given".to_owned(),
        ));
    }

    if options.minimize {
        let result = fire(calculator, &working, &options.optimize)?;
        steps.push(PrepareStep::new(
            "minimize",
            true,
            format!(
                "FIRE: {} steps, energy {:.6} eV, max force {:.6} eV/A, converged {}",
                result.steps, result.energy, result.max_force, result.converged
            ),
        ));
        energy = result.energy;
        working = result.apply_to(&working)?;
        optimized = Some(result);
    } else {
        steps.push(PrepareStep::new(
            "minimize",
            false,
            "skipped: geometry left as given".to_owned(),
        ));
        energy = calculator.energy(&working)?;
    }

    if options.charges {
        let charges = qeq_charges(&working, &options.charge_options)?;
        let total = charges.total;
        let n_atoms = charges.len();
        working = charges.apply_to(&working)?;
        steps.push(PrepareStep::new(
            "charges",
            true,
            format!("QEq: {n_atoms} P1 charges, total {total:.6} e"),
        ));
        computed_charges = Some(charges);
    } else {
        steps.push(PrepareStep::new(
            "charges",
            false,
            "skipped: site charges left as given".to_owned(),
        ));
    }

    if options.ewald {
        let cell = working.p1()?;
        let per_atom = p1_charges(&working, &cell, computed_charges.as_ref());
        let value = ewald_energy(&working, &per_atom, &options.ewald_options)?;
        steps.push(PrepareStep::new(
            "ewald",
            true,
            format!("Ewald: {} P1 charges, energy {value:.9} eV", per_atom.len()),
        ));
        ewald_value = Some(value);
    } else {
        steps.push(PrepareStep::new(
            "ewald",
            false,
            "skipped: Coulomb energy not requested".to_owned(),
        ));
    }

    if options.vibrations {
        let p1 = p1_view(&working)?;
        let hessian_options = HessianOptions::default();
        let matrix = hessian(calculator, &p1, &hessian_options)?;
        let modes = vibrational_modes(&matrix, &p1, &hessian_options)?;
        steps.push(PrepareStep::new(
            "vibrations",
            true,
            format!(
                "Gamma Hessian: {} atoms, {} modes, {} imaginary",
                p1.sites().len(),
                modes.frequencies.len(),
                modes.imaginary.len()
            ),
        ));
        frequencies = modes.frequencies;
    } else {
        steps.push(PrepareStep::new(
            "vibrations",
            false,
            "skipped: frequencies not requested".to_owned(),
        ));
    }

    let total_charge = computed_charges
        .as_ref()
        .map_or_else(|| working.total_charge(), |charges| charges.total);
    let mut close_contacts = Vec::new();
    if options.check_close_contacts {
        close_contacts = cleanup::close_contacts(&working)?;
        steps.push(PrepareStep::new(
            "close_contacts",
            true,
            format!("{} close contact(s)", close_contacts.len()),
        ));
    } else {
        steps.push(PrepareStep::new(
            "close_contacts",
            false,
            "skipped: close contacts not checked".to_owned(),
        ));
    }

    let mut charge_balance = None;
    if options.check_charge_balance {
        let balance = cleanup::charge_balance(&working)?;
        steps.push(PrepareStep::new(
            "charge_balance",
            true,
            if balance.neutral {
                "cell is neutral".to_owned()
            } else {
                format!("cell is not neutral: total {:.6} e", balance.total_charge)
            },
        ));
        charge_balance = Some(balance);
    } else {
        steps.push(PrepareStep::new(
            "charge_balance",
            false,
            "skipped: charge balance not checked".to_owned(),
        ));
    }

    Ok(PrepareReport {
        structure: working,
        steps,
        energy,
        total_charge,
        ewald_energy: ewald_value,
        frequencies,
        optimized,
        removed_sites,
        close_contacts,
        charge_balance,
    })
}

/// One charge per P1 atom, from the solver when available and from the sites
/// otherwise.
fn p1_charges(structure: &Structure, cell: &P1Cell, computed: Option<&Charges>) -> Vec<f64> {
    if let Some(charges) = computed
        && charges.len() == cell.n_atoms()
    {
        return charges.per_atom.clone();
    }
    cell.atoms
        .iter()
        .map(|atom| {
            structure
                .sites()
                .get(atom.site_idx)
                .and_then(Site::charge)
                .unwrap_or(0.0)
        })
        .collect()
}

/// A P1 copy of `structure` carrying each parent site's stored charge.
///
/// The Hessian and mode analysis need a P1 structure whose site count equals
/// its P1 atom count; rebuilding here also carries the charges assigned by the
/// previous step without changing the caller's structure.
fn p1_view(structure: &Structure) -> Result<Structure, Error> {
    let cell = structure.p1()?;
    let mut builder = StructureBuilder::new(structure.name(), cell.lattice);
    for (index, atom) in cell.atoms.iter().enumerate() {
        let mut site = Site::new(
            format!("{}{}", atom.element, index + 1),
            &atom.element,
            atom.frac,
        )?
        .with_occupancy(atom.occupancy)?;
        if let Some(charge) = structure.sites().get(atom.site_idx).and_then(Site::charge) {
            site = site.with_charge(charge)?;
        }
        builder = builder.site(site);
    }
    for (key, value) in structure.metadata() {
        builder = builder.metadata(key.clone(), value.clone());
    }
    builder.build().map_err(Error::from)
}
