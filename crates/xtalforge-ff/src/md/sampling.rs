//! Samples and trajectory frames recorded during molecular dynamics.

use serde::{Deserialize, Serialize};
use xtalforge_core::Lattice;

use crate::md::MdOptions;
use crate::trajectory::{Trajectory, TrajectoryFrame, TrajectoryOptions};

/// `1 GPa` in `eV/Angstrom^3`, matching the optimizer boundary.
const GPA_TO_EV_PER_A3: f64 = 1.0 / 160.217_663_4;

/// One recorded molecular-dynamics sample.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct MdSample {
    /// Integration step.
    pub step: usize,
    /// Simulated time in femtoseconds.
    pub time_fs: f64,
    /// Potential energy in eV.
    pub potential_energy: f64,
    /// Kinetic energy in eV.
    pub kinetic_energy: f64,
    /// Total energy in eV.
    pub total_energy: f64,
    /// Instantaneous temperature in kelvin.
    pub temperature: f64,
    /// Cell volume in Angstrom^3.
    pub volume: f64,
    /// Molecular virial pressure in GPa.
    pub pressure_gpa: f64,
    /// `potential_energy + pressure * volume` in eV.
    pub enthalpy: f64,
}

/// Everything needed to write one sample.
pub(crate) struct SampleInput {
    /// Integration step.
    pub(crate) step: usize,
    /// Simulated time in femtoseconds.
    pub(crate) time_fs: f64,
    /// Potential energy in eV.
    pub(crate) potential: f64,
    /// Kinetic energy in eV.
    pub(crate) kinetic: f64,
    /// Instantaneous temperature in kelvin.
    pub(crate) temperature: f64,
    /// Cell volume in Angstrom^3.
    pub(crate) volume: f64,
    /// Molecular virial pressure in GPa.
    pub(crate) pressure_gpa: f64,
    /// Largest force on any atom in eV/Angstrom.
    pub(crate) max_force: f64,
    /// Cell at this step.
    pub(crate) lattice: Lattice,
}

/// Collects samples, trajectory frames and their means.
pub(crate) struct Recorder {
    sample_every: usize,
    cadence: TrajectoryOptions,
    external_pressure_gpa: f64,
    samples: Vec<MdSample>,
    trajectory: Trajectory,
    temperature_sum: f64,
    pressure_sum: f64,
    counted: usize,
}

impl Recorder {
    /// A recorder for one run.
    pub(crate) fn new(options: &MdOptions) -> Self {
        Self {
            sample_every: options.sample_every,
            cadence: options.trajectory,
            external_pressure_gpa: options.ensemble.pressure_gpa().unwrap_or(0.0),
            samples: Vec::new(),
            trajectory: Trajectory::new(),
            temperature_sum: 0.0,
            pressure_sum: 0.0,
            counted: 0,
        }
    }

    /// Whether the run records a sample at `step`.
    pub(crate) fn sample_due(&self, step: usize, last: usize) -> bool {
        step == 0
            || step == last
            || (self.sample_every != 0 && step.is_multiple_of(self.sample_every))
    }

    /// Records one sample and its trajectory frame.
    pub(crate) fn record(&mut self, input: SampleInput, frac: &[[f64; 3]]) -> MdSample {
        let enthalpy =
            input.potential + self.external_pressure_gpa * GPA_TO_EV_PER_A3 * input.volume;
        let sample = MdSample {
            step: input.step,
            time_fs: input.time_fs,
            potential_energy: input.potential,
            kinetic_energy: input.kinetic,
            total_energy: input.potential + input.kinetic,
            temperature: input.temperature,
            volume: input.volume,
            pressure_gpa: input.pressure_gpa,
            enthalpy,
        };
        self.samples.push(sample.clone());
        self.temperature_sum += input.temperature;
        self.pressure_sum += input.pressure_gpa;
        self.counted += 1;
        if self.cadence.wants(input.step) {
            self.trajectory.push(TrajectoryFrame {
                step: input.step,
                point: None,
                time_fs: input.time_fs,
                energy: input.potential,
                enthalpy,
                max_force: input.max_force,
                volume: input.volume,
                lattice: input.lattice,
                positions: self.cadence.positions.then(|| frac.to_vec()),
            });
        }
        sample
    }

    /// Mean temperature and mean pressure over the recorded samples.
    pub(crate) fn means(&self) -> (f64, f64) {
        if self.counted == 0 {
            return (0.0, 0.0);
        }
        (
            self.temperature_sum / self.counted as f64,
            self.pressure_sum / self.counted as f64,
        )
    }

    /// Recorded samples.
    pub(crate) fn samples(&self) -> &[MdSample] {
        &self.samples
    }

    /// Recorded trajectory.
    pub(crate) fn into_trajectory(self) -> Trajectory {
        self.trajectory
    }
}
