//! Physical constants and unit conversions used by vibrations and
//! thermochemistry.
//!
//! CODATA 2018 values, kept together so the frequency, intensity and
//! thermochemistry code cannot drift apart.  All public results stay in the
//! crate's boundary units: wavenumbers in cm^-1, energies in eV, forces in
//! eV/Angstrom and masses in unified atomic mass units.

/// Planck constant, J s.
pub(crate) const PLANCK_J_S: f64 = 6.626_070_15e-34;

/// Boltzmann constant, J/K.
pub(crate) const BOLTZMANN_J_K: f64 = 1.380_649e-23;

/// Boltzmann constant, eV/K.
pub(crate) const BOLTZMANN_EV_K: f64 = 8.617_333_262_145_179e-5;

/// Unified atomic mass unit, kg.
pub(crate) const AMU_KG: f64 = 1.660_539_066_60e-27;

/// One bar in pascal.
pub(crate) const BAR_PA: f64 = 1.0e5;

/// One wavenumber as a photon energy, eV (`h c` in eV cm).
pub const CM_TO_EV: f64 = 1.239_841_984_332_003e-4;

/// `sqrt(eV / (Angstrom^2 amu))` converted to cm^-1.
///
/// This is the factor that turns a mass-weighted Hessian eigenvalue in
/// eV/(Angstrom^2 amu) into a harmonic frequency in wavenumbers.
pub const CM_PER_SQRT_EV_PER_A2_AMU: f64 = 521.470_898_372_506_6;

/// `(e / sqrt(amu))^2` to km/mol for an IR intensity.
///
/// `42.2561` converts `(D/Angstrom)^2/amu` to km/mol and `4.80320` is the
/// Debye per elementary-charge Angstrom, squared here because the intensity is
/// proportional to the square of the dipole derivative.
pub const E2_PER_AMU_TO_KM_PER_MOL: f64 = 42.2561 * 4.803_20 * 4.803_20;

/// One amu Angstrom^2 in kg m^2, for rotational partition functions.
pub(crate) const AMU_ANGSTROM2_TO_KG_M2: f64 = AMU_KG * 1.0e-20;
