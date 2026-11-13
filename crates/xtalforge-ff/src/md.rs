//! Molecular dynamics with the UFF calculator (Task 14.2.5).
//!
//! Four ensembles are available:
//!
//! * [`Ensemble::Nve`]: velocity Verlet, constant energy.
//! * [`Ensemble::NvtLangevin`]: velocity Verlet with an Ornstein-Uhlenbeck
//!   Langevin thermostat, constant temperature.
//! * [`Ensemble::NvtNoseHoover`]: velocity Verlet with a single Nose-Hoover
//!   thermostat (not a chain), constant temperature.
//! * [`Ensemble::Npt`]: Langevin thermostat plus an Andersen Monte-Carlo
//!   volume move, constant temperature and pressure.
//!
//! **Ensembles and what they mean.**  `NVE` conserves the total energy up to
//! the integrator error.  The two thermostats sample the canonical ensemble;
//! the Nose-Hoover thermostat is deterministic apart from the seeded initial
//! velocities, the Langevin thermostat uses the seed for its noise.  `NPT`
//! proposes isotropic volume moves at a fixed interval and accepts them with
//! the Andersen criterion
//!
//! ```text
//! min(1, exp(-(dE + P_ext dV - (N + 1) kT ln(V'/V)) / kT)),
//! ```
//!
//! which is a correct constant-pressure sampler; velocities and positions are
//! scaled by `(V'/V)^(1/3)` on acceptance.
//!
//! **Units.**  Angstrom, eV, atomic mass units, femtoseconds and kelvin.  The
//! pressure in [`MdSample`] and [`Ensemble::Npt`] is GPa.
//!
//! **Documented limitations.**
//!
//! * The structure must be P1; expand it with `Structure::as_p1()` first.  A
//!   symmetric structure would need the same force projection the optimizer
//!   does, which molecular dynamics does not implement.
//! * Interactions use the calculator as given.  UFF's Lennard-Jones term uses
//!   the minimum image inside the cell and has no long-range or tail
//!   correction, so [`MdSample::pressure_gpa`] is the molecular virial
//!   pressure, not a barostat-quality equation of state.
//! * A single Nose-Hoover thermostat is not ergodic for every system; use the
//!   Langevin thermostat when in doubt.
//! * Determinism: the same seed, options, crate and platform give the same
//!   trajectory.  Transcendental functions and the Metropolis decision are
//!   not claimed to be bit-identical across platforms, and the thermostat
//!   noise is an Irwin-Hall sum of twelve uniforms rather than a true
//!   Gaussian, so a run is reproducible but not a reference normal sampler.
//! * Velocities start from a Maxwell-Boltzmann draw at
//!   [`MdOptions::initial_temperature`] and the centre of mass is put at rest.
//!   The temperature uses `3N - 3` degrees of freedom.

use serde::{Deserialize, Serialize};
use xtalforge_core::Structure;

use crate::Error;
use crate::calculator::Calculator;
use crate::run::RunControl;
use crate::trajectory::{Trajectory, TrajectoryOptions};

mod barostat;
mod integrator;
mod random;
mod sampling;
mod thermostat;

pub use sampling::MdSample;

/// Boltzmann constant in eV per kelvin.
pub const KB_EV_PER_K: f64 = 8.617_333_262_145e-5;

/// One atomic mass unit in `eV fs^2 / Angstrom^2`: `1 amu = 103.6427
/// eV fs^2/Angstrom^2`.
///
/// Multiplying a mass in amu by a squared speed in `Angstrom^2/fs^2` and by
/// this constant gives an energy in eV, and a force in `eV/Angstrom` on a mass
/// in amu gives an acceleration in `Angstrom/fs^2` after dividing by it.
pub const AMU_TO_EV_FS2_PER_A2: f64 = 1.0 / 9.648_533_2e-3;

/// Which thermodynamic ensemble a run samples.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
#[serde(tag = "kind", rename_all = "snake_case", deny_unknown_fields)]
pub enum Ensemble {
    /// Constant energy; velocity Verlet without a thermostat.
    Nve,
    /// Constant temperature through a Langevin thermostat.
    NvtLangevin {
        /// Target temperature in kelvin.
        temperature: f64,
        /// Friction in inverse femtoseconds.
        friction: f64,
    },
    /// Constant temperature through a single Nose-Hoover thermostat.
    NvtNoseHoover {
        /// Target temperature in kelvin.
        temperature: f64,
        /// Thermostat relaxation time in femtoseconds.
        relaxation_fs: f64,
    },
    /// Constant temperature and pressure through a Langevin thermostat and an
    /// Andersen Monte-Carlo volume move.
    Npt {
        /// Target temperature in kelvin.
        temperature: f64,
        /// External pressure in GPa.
        pressure_gpa: f64,
        /// Langevin friction in inverse femtoseconds.
        friction: f64,
        /// Steps between volume moves.
        barostat_interval: usize,
        /// Relative standard deviation of one volume proposal.
        volume_move: f64,
    },
}

impl Ensemble {
    /// Stable short name: `"nve"`, `"nvt-langevin"`, `"nvt-nose-hoover"` or
    /// `"npt"`.
    #[must_use]
    pub const fn label(&self) -> &'static str {
        match self {
            Self::Nve => "nve",
            Self::NvtLangevin { .. } => "nvt-langevin",
            Self::NvtNoseHoover { .. } => "nvt-nose-hoover",
            Self::Npt { .. } => "npt",
        }
    }

    /// Target temperature, when the ensemble has one.
    #[must_use]
    pub const fn temperature(&self) -> Option<f64> {
        match self {
            Self::Nve => None,
            Self::NvtLangevin { temperature, .. }
            | Self::NvtNoseHoover { temperature, .. }
            | Self::Npt { temperature, .. } => Some(*temperature),
        }
    }

    /// External pressure in GPa, when the ensemble has one.
    #[must_use]
    pub const fn pressure_gpa(&self) -> Option<f64> {
        match self {
            Self::Npt { pressure_gpa, .. } => Some(*pressure_gpa),
            _ => None,
        }
    }

    /// Checks the ensemble parameters.
    ///
    /// # Errors
    ///
    /// Returns [`Error::InvalidInput`] for a non-positive temperature, a
    /// negative friction, a non-positive relaxation time, a non-finite
    /// pressure, a zero barostat interval or a volume move outside `(0, 1]`.
    pub fn validate(&self) -> Result<(), Error> {
        match self {
            Self::Nve => Ok(()),
            Self::NvtLangevin {
                temperature,
                friction,
            } => {
                positive_temperature(*temperature)?;
                non_negative_friction(*friction)
            }
            Self::NvtNoseHoover {
                temperature,
                relaxation_fs,
            } => {
                positive_temperature(*temperature)?;
                if !relaxation_fs.is_finite() || *relaxation_fs <= 0.0 {
                    return Err(Error::InvalidInput(format!(
                        "Nose-Hoover relaxation time must be positive and finite, got {relaxation_fs}"
                    )));
                }
                Ok(())
            }
            Self::Npt {
                temperature,
                pressure_gpa,
                friction,
                barostat_interval,
                volume_move,
            } => {
                positive_temperature(*temperature)?;
                non_negative_friction(*friction)?;
                if !pressure_gpa.is_finite() {
                    return Err(Error::InvalidInput(format!(
                        "external pressure must be finite, got {pressure_gpa}"
                    )));
                }
                if *barostat_interval == 0 {
                    return Err(Error::InvalidInput(
                        "the barostat interval must be at least one step".to_owned(),
                    ));
                }
                if !volume_move.is_finite() || *volume_move <= 0.0 || *volume_move > 1.0 {
                    return Err(Error::InvalidInput(format!(
                        "volume move must be in (0, 1], got {volume_move}"
                    )));
                }
                Ok(())
            }
        }
    }
}

/// Options for one molecular-dynamics run.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct MdOptions {
    /// Ensemble to sample.
    pub ensemble: Ensemble,
    /// Integration time step in femtoseconds.
    pub time_step_fs: f64,
    /// Number of integration steps.
    pub steps: usize,
    /// Seed for the initial velocities and every thermostat/barostat draw.
    ///
    /// Two runs with the same seed, options, structure and calculator produce
    /// the same trajectory.
    pub seed: u64,
    /// Temperature of the initial velocity draw, in kelvin; zero starts at
    /// rest.
    pub initial_temperature: f64,
    /// Record a sample every this many steps; `0` keeps only the first and
    /// last sample.
    pub sample_every: usize,
    /// Trajectory recording cadence; off by default.
    pub trajectory: TrajectoryOptions,
    /// Wrap fractional coordinates into `[0, 1)` after every step.
    pub wrap: bool,
}

impl Default for MdOptions {
    fn default() -> Self {
        Self {
            ensemble: Ensemble::Nve,
            time_step_fs: 0.5,
            steps: 1000,
            seed: 1,
            initial_temperature: 300.0,
            sample_every: 10,
            trajectory: TrajectoryOptions::default(),
            wrap: true,
        }
    }
}

impl MdOptions {
    /// Checks that the options are usable.
    ///
    /// # Errors
    ///
    /// Returns [`Error::InvalidInput`] for a non-positive or non-finite time
    /// step, zero steps, a negative or non-finite initial temperature, or
    /// invalid ensemble parameters.
    pub fn validate(&self) -> Result<(), Error> {
        if !self.time_step_fs.is_finite() || self.time_step_fs <= 0.0 {
            return Err(Error::InvalidInput(format!(
                "time step must be positive and finite, got {}",
                self.time_step_fs
            )));
        }
        if self.steps == 0 {
            return Err(Error::InvalidInput(
                "molecular dynamics needs at least one step".to_owned(),
            ));
        }
        if !self.initial_temperature.is_finite() || self.initial_temperature < 0.0 {
            return Err(Error::InvalidInput(format!(
                "initial temperature must be non-negative and finite, got {}",
                self.initial_temperature
            )));
        }
        self.ensemble.validate()
    }
}

/// What one molecular-dynamics run produced.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct MdResult {
    /// Final structure.
    pub structure: Structure,
    /// Ensemble that was sampled.
    pub ensemble: Ensemble,
    /// Recorded samples, in step order.
    pub samples: Vec<MdSample>,
    /// Total energy of the initial state in eV.
    pub initial_energy: f64,
    /// Total energy of the final state in eV.
    pub final_energy: f64,
    /// Mean temperature over the samples, in kelvin.
    pub mean_temperature: f64,
    /// Mean pressure over the samples, in GPa.
    pub mean_pressure_gpa: f64,
    /// Completed integration steps.
    pub steps: usize,
    /// Whether a [`RunControl`] stopped the run.
    pub cancelled: bool,
    /// Recorded trajectory frames, empty unless a cadence was set.
    pub trajectory: Trajectory,
}

impl MdResult {
    /// Potential energies of the samples, in eV.
    #[must_use]
    pub fn potential_energies(&self) -> Vec<f64> {
        self.samples
            .iter()
            .map(|sample| sample.potential_energy)
            .collect()
    }

    /// Total energies of the samples, in eV.
    #[must_use]
    pub fn total_energies(&self) -> Vec<f64> {
        self.samples
            .iter()
            .map(|sample| sample.total_energy)
            .collect()
    }
}

/// Runs molecular dynamics with the default control (no cancellation, no
/// progress callback).
///
/// # Errors
///
/// Returns the errors of [`molecular_dynamics_with`].
pub fn molecular_dynamics(
    calculator: &dyn Calculator,
    structure: &Structure,
    options: &MdOptions,
) -> Result<MdResult, Error> {
    molecular_dynamics_with(calculator, structure, options, &mut RunControl::new())
}

/// Runs molecular dynamics, reporting progress and honouring cancellation.
///
/// # Errors
///
/// Returns [`Error::InvalidInput`] for invalid options, a non-P1 structure, a
/// cell without positive volume or a site without a mass, and any calculator
/// error raised at a step.
pub fn molecular_dynamics_with(
    calculator: &dyn Calculator,
    structure: &Structure,
    options: &MdOptions,
    control: &mut RunControl<'_>,
) -> Result<MdResult, Error> {
    integrator::run(calculator, structure, options, control)
}

fn positive_temperature(temperature: f64) -> Result<(), Error> {
    if !temperature.is_finite() || temperature <= 0.0 {
        return Err(Error::InvalidInput(format!(
            "thermostat temperature must be positive and finite, got {temperature}"
        )));
    }
    Ok(())
}

fn non_negative_friction(friction: f64) -> Result<(), Error> {
    if !friction.is_finite() || friction < 0.0 {
        return Err(Error::InvalidInput(format!(
            "thermostat friction must be non-negative and finite, got {friction}"
        )));
    }
    Ok(())
}

/// Degrees of freedom for the temperature: `3N - 3`, or `3` for one atom.
fn degrees_of_freedom(n_atoms: usize) -> f64 {
    if n_atoms > 1 {
        (3 * n_atoms - 3) as f64
    } else {
        3.0
    }
}
