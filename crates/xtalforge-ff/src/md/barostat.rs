//! Andersen Monte-Carlo barostat used by the NPT ensemble.
//!
//! Every `interval` steps the run proposes an isotropic volume change
//! `V' = V (1 + amplitude * normal)`, evaluates the potential energy in the
//! scaled cell and accepts the move with
//!
//! ```text
//! min(1, exp(-(dE + P_ext dV - (N + 1) kT ln(V'/V)) / kT)).
//! ```
//!
//! The `(N + 1)` term is the Jacobian of the scaled coordinates plus the
//! momentum rescaling; positions and velocities are scaled by `(V'/V)^(1/3)`
//! when the move is accepted.  This samples the isothermal-isobaric ensemble
//! without an auxiliary piston mass; `amplitude` is the only tuning knob and a
//! few percent is a good starting point for a molecular crystal.

use crate::md::random::SeededRng;
use crate::md::{Ensemble, KB_EV_PER_K};

/// `1 GPa` in `eV/Angstrom^3`, matching the optimizer boundary.
const GPA_TO_EV_PER_A3: f64 = 1.0 / 160.217_663_4;

/// Isotropic volume moves at a fixed interval.
pub(crate) struct Barostat {
    interval: usize,
    amplitude: f64,
    temperature: f64,
    /// External pressure in eV/Angstrom^3.
    pressure: f64,
}

/// One proposed volume change.
pub(crate) struct VolumeProposal {
    /// `V_new / V_old`.
    pub(crate) factor: f64,
    /// Proposed volume in Angstrom^3.
    pub(crate) new_volume: f64,
    /// Current volume in Angstrom^3.
    pub(crate) old_volume: f64,
    /// Current potential energy in eV.
    pub(crate) potential_energy: f64,
}

impl Barostat {
    /// The barostat of an ensemble, when it has one.
    pub(crate) fn new(ensemble: &Ensemble) -> Option<Self> {
        match ensemble {
            Ensemble::Npt {
                temperature,
                pressure_gpa,
                barostat_interval,
                volume_move,
                ..
            } => Some(Self {
                interval: *barostat_interval,
                amplitude: *volume_move,
                temperature: *temperature,
                pressure: pressure_gpa * GPA_TO_EV_PER_A3,
            }),
            _ => None,
        }
    }

    /// Steps between volume moves.
    pub(crate) fn interval(&self) -> usize {
        self.interval
    }

    /// Proposes one volume change; the caller evaluates the new energy.
    pub(crate) fn propose(
        &self,
        rng: &mut SeededRng,
        volume: f64,
        potential_energy: f64,
    ) -> VolumeProposal {
        let factor = (1.0 + self.amplitude * rng.normal()).clamp(0.25, 4.0);
        VolumeProposal {
            factor,
            new_volume: volume * factor,
            old_volume: volume,
            potential_energy,
        }
    }

    /// Metropolis decision for a proposal.
    pub(crate) fn accepted(
        &self,
        rng: &mut SeededRng,
        proposal: &VolumeProposal,
        new_potential: f64,
        n_atoms: usize,
    ) -> bool {
        let thermal = KB_EV_PER_K * self.temperature;
        if thermal <= 0.0 {
            return false;
        }
        let delta = new_potential - proposal.potential_energy
            + self.pressure * (proposal.new_volume - proposal.old_volume)
            - (n_atoms as f64 + 1.0) * thermal * proposal.factor.ln();
        if !delta.is_finite() {
            return false;
        }
        if delta <= 0.0 {
            return true;
        }
        rng.uniform() < (-delta / thermal).exp()
    }
}
