//! Langevin and Nose-Hoover thermostats for molecular dynamics.

use crate::md::random::SeededRng;
use crate::md::{AMU_TO_EV_FS2_PER_A2, Ensemble, KB_EV_PER_K};

/// The thermostat attached to one run.
pub(crate) enum Thermostat {
    /// No thermostat.
    None,
    /// Ornstein-Uhlenbeck (Langevin) friction and noise.
    Langevin { temperature: f64, friction: f64 },
    /// A single Nose-Hoover thermostat.
    NoseHoover {
        temperature: f64,
        mass: f64,
        dof: f64,
        xi: f64,
    },
}

impl Thermostat {
    /// Builds the thermostat of an ensemble.
    pub(crate) fn new(ensemble: &Ensemble, n_atoms: usize) -> Self {
        match ensemble {
            Ensemble::Nve => Self::None,
            Ensemble::NvtLangevin {
                temperature,
                friction,
            } => Self::Langevin {
                temperature: *temperature,
                friction: *friction,
            },
            Ensemble::NvtNoseHoover {
                temperature,
                relaxation_fs,
            } => {
                let dof = super::degrees_of_freedom(n_atoms);
                Self::NoseHoover {
                    temperature: *temperature,
                    mass: dof * KB_EV_PER_K * temperature * relaxation_fs * relaxation_fs,
                    dof,
                    xi: 0.0,
                }
            }
            Ensemble::Npt {
                temperature,
                friction,
                ..
            } => Self::Langevin {
                temperature: *temperature,
                friction: *friction,
            },
        }
    }

    /// Applies the first half of the thermostat before the drift step.
    pub(crate) fn pre(
        &mut self,
        velocities: &mut [[f64; 3]],
        masses: &[f64],
        dt: f64,
        rng: &mut SeededRng,
    ) {
        match self {
            Self::None => {}
            Self::Langevin {
                temperature,
                friction,
            } => {
                let decay = (-*friction * dt).exp();
                let noise = (1.0 - decay * decay).sqrt();
                for (velocity, mass) in velocities.iter_mut().zip(masses) {
                    if *mass <= 0.0 {
                        continue;
                    }
                    let sigma =
                        (KB_EV_PER_K * *temperature / (AMU_TO_EV_FS2_PER_A2 * *mass)).sqrt();
                    for value in velocity.iter_mut() {
                        *value = decay * *value + noise * sigma * rng.normal();
                    }
                }
            }
            Self::NoseHoover { xi, .. } => {
                let scale = (-0.5 * dt * *xi).exp();
                for velocity in velocities.iter_mut() {
                    for value in velocity.iter_mut() {
                        *value *= scale;
                    }
                }
            }
        }
    }

    /// Applies the second half of the thermostat after the second half kick.
    ///
    /// `kinetic` is the kinetic energy in eV before the second half.
    pub(crate) fn post(&mut self, velocities: &mut [[f64; 3]], dt: f64, kinetic: f64) {
        match self {
            Self::None | Self::Langevin { .. } => {}
            Self::NoseHoover {
                temperature,
                mass,
                dof,
                xi,
            } => {
                if *mass <= 0.0 {
                    return;
                }
                let drive = 2.0 * kinetic - *dof * KB_EV_PER_K * *temperature;
                *xi += dt * drive / *mass;
                let scale = (-0.5 * dt * *xi).exp();
                for velocity in velocities.iter_mut() {
                    for value in velocity.iter_mut() {
                        *value *= scale;
                    }
                }
            }
        }
    }
}
