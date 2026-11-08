//! Native Universal Force Field and its UFF4MOF extension.
//!
//! The implementation follows Crystal-Builder `xtal/ff/uff/`:
//! [`params`] holds the 218-row table, [`typer`] decides which row an atom
//! gets, [`terms`], [`torsion`], [`inversion`] and [`vdw`] evaluate the
//! energy expressions and their analytic gradients, and [`model`] combines
//! them over a fixed P1 topology. [`calculator::UffCalculator`] is the
//! registry-facing engine.
//!
//! Units inside are the reference units (Angstrom, degrees, kcal/mol); the
//! public [`crate::Calculator`] boundary is eV and eV/Angstrom.

pub mod calculator;
pub mod inversion;
pub mod model;
pub mod params;
pub mod terms;
pub mod torsion;
pub mod typer;
pub mod vdw;

pub(crate) mod coordination;
pub(crate) mod geometry;
pub(crate) mod math;
pub(crate) mod perception;
pub(crate) mod rules;
pub(crate) mod topology;

pub use calculator::UffCalculator;
pub use inversion::{InversionParameters, InversionTerm, inversion_parameters};
pub use model::UffModel;
pub use params::{ParameterSet, UffParams};
pub use terms::{
    AngleForm, AngleTerm, BondTerm, angle_force_constant, bond_force_constant, natural_bond_length,
    vdw_pair,
};
pub use torsion::{TorsionParameters, TorsionTerm, hybridisation, torsion_parameters};
pub use typer::{AtomType, Confidence, Typing};
pub use vdw::VanDerWaalsTerm;

/// Default van der Waals cutoff in Angstrom, the reference value.
pub const DEFAULT_VDW_CUTOFF: f64 = 12.0;
