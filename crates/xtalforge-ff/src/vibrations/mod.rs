//! Vibrations, IR/Raman spectra, phonons, thermochemistry and quasi-harmonic
//! free energy (Tasks 4.8 and 19.4).
//!
//! The core path is a Gamma-point Hessian from central differences of analytic
//! forces, mass-weighted and diagonalized into normal modes.  The same
//! frequency list feeds five consumers:
//!
//! * [`ir_spectrum`] — fixed-QEq-charge dipole derivatives and a broadened
//!   stick spectrum (Task 4.8.2; the fixed-charge approximation is documented
//!   in the module).
//! * [`ir_spectrum_with_born`] — the same spectrum from per-atom Born
//!   effective-charge tensors, supplied by an engine or derived from QEq/EQeq
//!   (Task 19.4.1).
//! * [`raman_spectrum`] — powder-averaged Raman activity from engine-provided
//!   per-mode tensors; no native polarizability model (Task 19.4.2).
//! * [`phonons`] — the Gamma-point dynamical matrix and q-mesh (Task 4.8.3),
//!   with [`phonon_dos`] and [`inelastic_neutron_dos`] as consumers (Task
//!   19.4.3).
//! * [`thermochemistry`] and [`ideal_gas_thermochemistry`] — harmonic
//!   oscillator terms, with translation and rigid rotation for a molecule
//!   (Task 4.8.4).
//!
//! [`quasi_harmonic_free_energy`] evaluates the free energy at caller-supplied
//! volumes and refines the minimum with a parabola; volume relaxation and
//! thermal expansion remain Task 4.8.5 follow-up work.
//!
//! Units at this boundary are cm^-1 for frequencies, eV for energies, eV/K for
//! entropy and Angstrom for geometry.  Where a result depends on an
//! approximation, the module comment says which one.

mod constants;
mod hessian;
mod ir;
mod ir_born;
mod modes;
mod phonon_dos;
mod phonons;
mod qha;
mod qmesh;
mod raman;
mod rigid;
mod spectrum;
mod tensors;
mod thermo;

pub use constants::{CM_PER_SQRT_EV_PER_A2_AMU, CM_TO_EV, E2_PER_AMU_TO_KM_PER_MOL};
pub use hessian::{DEFAULT_HESSIAN_STEP, DisplacementMode, Hessian, HessianOptions, hessian};
pub use ir::{
    DEFAULT_IR_FWHM, DEFAULT_IR_STEP, IrOptions, IrPeak, IrSpectrum, MAX_IR_POINTS, ir_spectrum,
};
pub use ir_born::{BornCharges, ir_spectrum_with_born};
pub use modes::{IMAGINARY_BELOW_CM, VibrationalModes, ZERO_MODE_TOLERANCE_CM, vibrational_modes};
pub use phonon_dos::{
    InelasticNeutronDos, NeutronWeights, inelastic_neutron_dos, mode_neutron_weights, phonon_dos,
};
pub use phonons::{DosOptions, PhononDos, PhononOptions, PhononResult, phonons};
pub use qha::{QhaFreeEnergy, QhaPoint, QhaResult, quasi_harmonic_free_energy};
pub use raman::{
    DEFAULT_RAMAN_FWHM, DEFAULT_RAMAN_STEP, MAX_RAMAN_POINTS, RamanOptions, RamanPeak,
    RamanSpectrum, RamanTensorSet, powder_average, raman_spectrum,
};
pub use spectrum::Broadening;
pub use tensors::{TensorFile, parse_tensor_file};
pub use thermo::{
    LINEAR_MOMENT_AMU_ANGSTROM2, ThermoOptions, ThermoResult, ideal_gas_thermochemistry,
    thermochemistry,
};
