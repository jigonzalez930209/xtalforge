//! Harmonic thermochemistry and the ideal-gas rigid-rotor extension
//! (Task 4.8.4).
//!
//! [`thermochemistry`] sums the harmonic oscillator partition function over a
//! frequency list.  Frequencies at or below [`ZERO_MODE_TOLERANCE_CM`] and
//! imaginary modes are excluded: their harmonic entropy diverges at zero, and
//! a saddle point has no equilibrium partition function.  The exclusion is
//! documented and counted nowhere else, because the requested
//! [`ThermoResult`] carries the thermodynamic values rather than the mode
//! bookkeeping; use [`VibrationalModes::imaginary`] when the counts matter.
//!
//! **Standard state and units.**  Energies are eV per cell (or per molecule),
//! entropy is eV/K, temperature is K and pressure is bar.  The vibrational
//! result uses the Helmholtz convention of a solid: `thermal_u` is the thermal
//! vibrational energy above the zero-point energy, `thermal_h` adds the
//! ideal-gas `k_B T` term, and `gibbs = zpe + thermal_h - T S`.  At the
//! standard state of 1 bar the `k_B T` term is about 0.026 eV at 298.15 K; for
//! a crystal it is negligible compared with typical cell energies, but it is
//! reported rather than silently dropped.
//!
//! [`ideal_gas_thermochemistry`] adds translation and rigid rotation for a
//! molecule.  The rotational symmetry number is taken as 1 and electronic
//! degeneracy is ignored; both are documented approximations.  The
//! vibrational-only entry point never uses the pressure, so a crystal run has
//! no spurious ideal-gas contribution.
//!
//! **Documented approximations.**  Harmonic frequencies are used without
//! quasi-RRHO or anharmonic corrections.  A soft mode below the zero-mode
//! tolerance is omitted rather than replaced by a free rotor; the reference
//! implementation's quasi-RRHO option is future work.

use serde::{Deserialize, Serialize};
use xtalforge_core::{Structure, element_info};

use crate::Error;
use crate::math::symmetric_eigen;

use super::constants::{
    AMU_ANGSTROM2_TO_KG_M2, AMU_KG, BAR_PA, BOLTZMANN_EV_K, BOLTZMANN_J_K, CM_TO_EV, PLANCK_J_S,
};
use super::modes::{VibrationalModes, ZERO_MODE_TOLERANCE_CM};

/// A principal moment below this value (amu Angstrom^2) means a linear
/// molecule.
pub const LINEAR_MOMENT_AMU_ANGSTROM2: f64 = 1.0e-3;

/// Options for one thermochemistry calculation.
#[derive(Debug, Clone, Copy, PartialEq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct ThermoOptions {
    /// Temperature in kelvin; must be positive.
    pub temperature: f64,
    /// Standard-state pressure in bar; used by the ideal-gas extension.
    pub pressure: f64,
}

impl Default for ThermoOptions {
    fn default() -> Self {
        Self {
            temperature: 298.15,
            pressure: 1.0,
        }
    }
}

impl ThermoOptions {
    /// Checks that the options are usable.
    ///
    /// # Errors
    ///
    /// Returns [`Error::InvalidInput`] for a non-positive or non-finite
    /// temperature or pressure.
    pub fn validate(&self) -> Result<(), Error> {
        if !self.temperature.is_finite() || self.temperature <= 0.0 {
            return Err(Error::InvalidInput(format!(
                "thermochemistry temperature must be positive and finite, got {}",
                self.temperature
            )));
        }
        if !self.pressure.is_finite() || self.pressure <= 0.0 {
            return Err(Error::InvalidInput(format!(
                "thermochemistry pressure must be positive and finite, got {}",
                self.pressure
            )));
        }
        Ok(())
    }
}

/// Thermal corrections from one frequency list.
///
/// All energies are eV per cell (or per molecule); `entropy` is eV/K.
/// `thermal_u` and `thermal_h` exclude the zero-point energy, which is
/// reported separately in `zpe`.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct ThermoResult {
    /// Zero-point energy in eV.
    pub zpe: f64,
    /// Thermal vibrational energy above the zero point in eV.
    pub thermal_u: f64,
    /// Thermal enthalpy above the zero point in eV.
    pub thermal_h: f64,
    /// Entropy in eV/K.
    pub entropy: f64,
    /// Gibbs energy in eV, `zpe + thermal_h - T S`.
    pub gibbs: f64,
}

impl ThermoResult {
    /// Helmholtz energy in eV, `zpe + thermal_u - T S`.
    ///
    /// This is the solid-state free energy; [`ThermoResult::gibbs`] adds the
    /// ideal-gas `k_B T` term already contained in `thermal_h`.
    #[must_use]
    pub fn helmholtz(&self, temperature: f64) -> f64 {
        self.zpe + self.thermal_u - temperature * self.entropy
    }
}

impl VibrationalModes {
    /// Thermochemistry of these frequencies at the given temperature.
    ///
    /// # Errors
    ///
    /// Returns the errors of [`thermochemistry`].
    pub fn thermochemistry(&self, options: &ThermoOptions) -> Result<ThermoResult, Error> {
        thermochemistry(&self.frequencies, options)
    }
}

/// Harmonic thermochemistry of one frequency list.
///
/// The list may come from a Gamma Hessian, a q-mesh with duplicate entries or
/// an external engine.  Positive modes above
/// [`ZERO_MODE_TOLERANCE_CM`] are summed; imaginary, zero and non-finite modes
/// do not enter the partition function.
///
/// # Errors
///
/// Returns [`Error::InvalidInput`] for invalid options or a non-finite
/// frequency.
pub fn thermochemistry(
    frequencies: &[f64],
    options: &ThermoOptions,
) -> Result<ThermoResult, Error> {
    options.validate()?;
    let temperature = options.temperature;
    let kt = BOLTZMANN_EV_K * temperature;
    let mut zpe = 0.0;
    let mut thermal_u = 0.0;
    let mut entropy = 0.0;
    for (index, frequency) in frequencies.iter().enumerate() {
        if !frequency.is_finite() {
            return Err(Error::InvalidInput(format!(
                "frequency {index} is not finite: {frequency}"
            )));
        }
        if *frequency <= ZERO_MODE_TOLERANCE_CM {
            continue;
        }
        let energy = CM_TO_EV * frequency;
        zpe += 0.5 * energy;
        let x = energy / kt;
        if x <= 0.0 {
            continue;
        }
        let occupation = if x > 700.0 { 0.0 } else { 1.0 / x.exp_m1() };
        thermal_u += energy * occupation;
        entropy += BOLTZMANN_EV_K * harmonic_entropy(x);
    }
    let thermal_h = thermal_u + kt;
    let gibbs = zpe + thermal_h - temperature * entropy;
    Ok(ThermoResult {
        zpe,
        thermal_u,
        thermal_h,
        entropy,
        gibbs,
    })
}

/// S / k_B of one harmonic oscillator at `x = h nu / k T`.
fn harmonic_entropy(x: f64) -> f64 {
    if x > 700.0 {
        return 0.0;
    }
    if x < 1.0e-9 {
        // x / (e^x - 1) -> 1 and -ln(1 - e^-x) -> -ln(x).
        return 1.0 - x.ln();
    }
    let occupation = 1.0 / x.exp_m1();
    let log_term = (-x).exp_m1().abs().ln();
    x * occupation - log_term
}

/// Ideal-gas translation and rigid rotation on top of the vibrational terms.
///
/// This entry point is for an isolated molecule.  It adds the translational
/// and rotational partition functions, uses [`ThermoOptions::pressure`] for the
/// translational volume and assumes a rotational symmetry number of one and no
/// electronic degeneracy.  A single atom has no rotational term; a linear
/// molecule uses the two-dimensional rotor expression.  The structure must be
/// P1 and match the modes.
///
/// # Errors
///
/// Returns [`Error::InvalidInput`] for a non-P1 structure, a mode/atom
/// mismatch, a dummy atom, or a non-positive mass; and the errors of
/// [`thermochemistry`] for bad options or frequencies.
pub fn ideal_gas_thermochemistry(
    structure: &Structure,
    modes: &VibrationalModes,
    options: &ThermoOptions,
) -> Result<ThermoResult, Error> {
    let mut result = thermochemistry(&modes.frequencies, options)?;
    if !structure.is_p1() {
        return Err(Error::InvalidInput(
            "ideal-gas thermochemistry needs a P1 structure".to_owned(),
        ));
    }
    let n_atoms = structure.sites().len();
    if modes.n_atoms() != n_atoms {
        return Err(Error::InvalidInput(format!(
            "modes describe {} atoms but the structure has {n_atoms}",
            modes.n_atoms()
        )));
    }
    if n_atoms == 0 {
        return Err(Error::InvalidInput(
            "an ideal gas needs at least one atom".to_owned(),
        ));
    }

    let lattice = *structure.lattice();
    let mut masses = Vec::with_capacity(n_atoms);
    let mut cart = Vec::with_capacity(n_atoms);
    for site in structure.sites() {
        let element = element_info(site.element())?;
        if element.is_dummy() || element.mass <= 0.0 {
            return Err(Error::InvalidInput(format!(
                "element {} has no positive mass",
                site.element()
            )));
        }
        masses.push(element.mass);
        cart.push(lattice.to_cart(site.frac()));
    }
    let total_mass: f64 = masses.iter().sum();
    if total_mass <= 0.0 {
        return Err(Error::InvalidInput(
            "an ideal gas needs a positive total mass".to_owned(),
        ));
    }
    let mut centre = [0.0; 3];
    for (mass, point) in masses.iter().zip(&cart) {
        centre[0] += mass * point[0];
        centre[1] += mass * point[1];
        centre[2] += mass * point[2];
    }
    for value in &mut centre {
        *value /= total_mass;
    }

    let temperature = options.temperature;
    let kt_j = BOLTZMANN_J_K * temperature;
    let h = PLANCK_J_S;
    let mass_kg = total_mass * AMU_KG;
    let volume = kt_j / (options.pressure * BAR_PA);
    let q_trans = (2.0 * std::f64::consts::PI * mass_kg * kt_j / (h * h)).powf(1.5) * volume;
    let s_trans = BOLTZMANN_EV_K * (q_trans.ln() + 2.5);
    let u_trans = 1.5 * BOLTZMANN_EV_K * temperature;
    result.thermal_u += u_trans;
    result.thermal_h += u_trans + BOLTZMANN_EV_K * temperature;
    result.entropy += s_trans;

    if n_atoms > 1 {
        let inertia = inertia_tensor(&masses, &cart, centre);
        let (moments, _) = symmetric_eigen(&inertia).ok_or_else(|| {
            Error::InvalidInput("the inertia tensor could not be diagonalized".to_owned())
        })?;
        let moments: Vec<f64> = moments
            .iter()
            .map(|value| (value * AMU_ANGSTROM2_TO_KG_M2).max(0.0))
            .collect();
        let linear = moments
            .first()
            .map(|value| *value / AMU_ANGSTROM2_TO_KG_M2 < LINEAR_MOMENT_AMU_ANGSTROM2)
            .unwrap_or(false);
        if linear {
            let moment = moments.last().copied().unwrap_or(0.0);
            if moment > 0.0 {
                let q_rot =
                    8.0 * std::f64::consts::PI * std::f64::consts::PI * moment * kt_j / (h * h);
                let s_rot = BOLTZMANN_EV_K * (q_rot.ln() + 1.0);
                let u_rot = BOLTZMANN_EV_K * temperature;
                result.thermal_u += u_rot;
                result.thermal_h += u_rot;
                result.entropy += s_rot;
            }
        } else {
            let product: f64 = moments.iter().product();
            if product > 0.0 {
                let q_rot = std::f64::consts::PI.sqrt()
                    * (8.0 * std::f64::consts::PI * std::f64::consts::PI * kt_j / (h * h))
                        .powf(1.5)
                    * product.sqrt();
                let s_rot = BOLTZMANN_EV_K * (q_rot.ln() + 1.5);
                let u_rot = 1.5 * BOLTZMANN_EV_K * temperature;
                result.thermal_u += u_rot;
                result.thermal_h += u_rot;
                result.entropy += s_rot;
            }
        }
    }
    result.gibbs = result.zpe + result.thermal_h - temperature * result.entropy;
    Ok(result)
}

fn inertia_tensor(masses: &[f64], cart: &[[f64; 3]], centre: [f64; 3]) -> Vec<Vec<f64>> {
    let mut inertia = vec![vec![0.0; 3]; 3];
    for (mass, point) in masses.iter().zip(cart) {
        let r = [
            point[0] - centre[0],
            point[1] - centre[1],
            point[2] - centre[2],
        ];
        let squared = r[0] * r[0] + r[1] * r[1] + r[2] * r[2];
        for (row, values) in inertia.iter_mut().enumerate() {
            for (column, slot) in values.iter_mut().enumerate() {
                let diagonal = if row == column { squared } else { 0.0 };
                let first = r.get(row).copied().unwrap_or(0.0);
                let second = r.get(column).copied().unwrap_or(0.0);
                *slot += mass * (diagonal - first * second);
            }
        }
    }
    inertia
}
