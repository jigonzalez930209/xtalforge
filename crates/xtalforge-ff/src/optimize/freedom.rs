//! Which cell quantities a variable-cell optimization may change (Task 14.2.1).
//!
//! The optimizer's six strain variables act on the base cell.  A [`CellFreedom`]
//! decides which of the six cell parameters survive that strain and whether the
//! volume may change at all:
//!
//! * a held length is copied back from the base cell,
//! * a held angle is copied back from the base cell,
//! * a held volume is met by scaling the free lengths.
//!
//! Scaling the free lengths is exact because the volume is `a b c` times a
//! function of the angles alone, so no held length or angle moves when the
//! volume is restored.  A held volume therefore needs at least one free
//! length; [`CellFreedom::validate`] refuses the alternative.
//!
//! The corrected cell is rebuilt from the six parameters and mapped back into
//! the base cell's orientation, so a rotated input cell keeps its orientation.
//! Because the correction is applied at every energy evaluation, a held
//! quantity is held exactly and not only to first order.
//!
//! When the cell gradient comes from an analytic stress (UFF provides one), a
//! constrained cell uses central differences of the constrained objective
//! instead; the stress of the unconstrained cell is not the gradient of the
//! constrained one.  That costs twelve extra energy evaluations per optimizer
//! iteration and keeps the held quantities consistent with what is reported.
//! [`CellFreedom::full`] keeps the analytic path unchanged.

use serde::{Deserialize, Serialize};
use xtalforge_core::Lattice;

use crate::Error;
use crate::math::mat3_mul;

/// Which cell lengths, angles and volume an optimization may change.
///
/// `lengths[i]` is `a`, `b`, `c`; `angles[i]` is `alpha`, `beta`, `gamma`;
/// `volume` says whether the cell volume may change.  The default is
/// [`CellFreedom::full`], which reproduces the variable-cell behaviour of
/// earlier releases.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct CellFreedom {
    /// Whether `a`, `b` and `c` may change.
    pub lengths: [bool; 3],
    /// Whether `alpha`, `beta` and `gamma` may change.
    pub angles: [bool; 3],
    /// Whether the volume may change; otherwise the cell is rescaled to the
    /// base volume after every update.
    pub volume: bool,
}

impl Default for CellFreedom {
    fn default() -> Self {
        Self::full()
    }
}

impl CellFreedom {
    /// Every length, angle and the volume are variables (the old behaviour).
    #[must_use]
    pub const fn full() -> Self {
        Self {
            lengths: [true; 3],
            angles: [true; 3],
            volume: true,
        }
    }

    /// Nothing may change; a variable-cell run with this freedom is an error.
    #[must_use]
    pub const fn none() -> Self {
        Self {
            lengths: [false; 3],
            angles: [false; 3],
            volume: false,
        }
    }

    /// Lengths are variables and the angles are held, as in an orthorhombic
    /// relaxation; the volume follows the lengths.
    #[must_use]
    pub const fn lengths_only() -> Self {
        Self {
            lengths: [true; 3],
            angles: [false; 3],
            volume: true,
        }
    }

    /// Every length and angle is a variable but the volume is held, so only
    /// the cell shape moves.
    #[must_use]
    pub const fn fixed_volume() -> Self {
        Self {
            lengths: [true; 3],
            angles: [true; 3],
            volume: false,
        }
    }

    /// Whether the six strain variables pass through untouched.
    #[must_use]
    pub fn is_full(&self) -> bool {
        self.volume && self.lengths.iter().all(|free| *free) && self.angles.iter().all(|free| *free)
    }

    /// Allows the length `a`, `b` or `c` (`index` 0, 1 or 2) to move.
    #[must_use]
    pub fn with_length(mut self, index: usize, free: bool) -> Self {
        if let Some(slot) = self.lengths.get_mut(index) {
            *slot = free;
        }
        self
    }

    /// Allows the angle `alpha`, `beta` or `gamma` (`index` 0, 1 or 2) to move.
    #[must_use]
    pub fn with_angle(mut self, index: usize, free: bool) -> Self {
        if let Some(slot) = self.angles.get_mut(index) {
            *slot = free;
        }
        self
    }

    /// Holds one length at its base value.
    #[must_use]
    pub fn hold_length(self, index: usize) -> Self {
        self.with_length(index, false)
    }

    /// Holds one angle at its base value.
    #[must_use]
    pub fn hold_angle(self, index: usize) -> Self {
        self.with_angle(index, false)
    }

    /// Holds the volume at its base value.
    #[must_use]
    pub fn hold_volume(mut self) -> Self {
        self.volume = false;
        self
    }

    /// Checks that the freedom can drive a variable-cell run.
    ///
    /// # Errors
    ///
    /// Returns [`Error::InvalidInput`] when nothing may change, which would
    /// make `cell: true` a no-op, or when the volume is held without a free
    /// length to absorb the rescale.
    pub fn validate(&self) -> Result<(), Error> {
        if self.is_full() {
            return Ok(());
        }
        let free_lengths = self.lengths.iter().any(|free| *free);
        let free_angles = self.angles.iter().any(|free| *free);
        if !self.volume && !free_lengths && !free_angles {
            return Err(Error::InvalidInput(
                "the cell is a variable but every cell quantity is held".to_owned(),
            ));
        }
        if !self.volume && !free_lengths {
            return Err(Error::InvalidInput(
                "holding the volume needs at least one free length".to_owned(),
            ));
        }
        Ok(())
    }

    /// Applies the freedom to a strained candidate cell.
    ///
    /// Held parameters are taken from `base`.  A held volume is restored by
    /// scaling the free lengths, which is exact because the volume is
    /// `a b c` times a function of the angles alone; held lengths and angles
    /// therefore stay exact too.  The result has the same orientation as
    /// `base`.
    ///
    /// # Errors
    ///
    /// Returns a core error when the corrected parameters are geometrically
    /// impossible or the corrected cell is singular.
    pub(crate) fn apply(&self, base: &Lattice, candidate: Lattice) -> Result<Lattice, Error> {
        self.apply_keeping(base, candidate, None)
    }

    /// [`CellFreedom::apply`] that leaves one free length out of the volume
    /// rescale.
    ///
    /// The relaxed-scan driver uses `keep_length` for the cell parameter it is
    /// scanning: the target must survive a held volume, so the other free
    /// lengths absorb the rescale.
    ///
    /// # Errors
    ///
    /// Returns the errors of [`CellFreedom::apply`], or
    /// [`Error::InvalidInput`] when the held volume leaves no length to scale.
    pub(crate) fn apply_keeping(
        &self,
        base: &Lattice,
        candidate: Lattice,
        keep_length: Option<usize>,
    ) -> Result<Lattice, Error> {
        if self.is_full() {
            return Ok(candidate);
        }
        let base_parameters = base.parameters();
        let mut parameters = candidate.parameters();
        for index in 0..3 {
            if self.lengths.get(index).copied().unwrap_or(false) {
                continue;
            }
            if let (Some(slot), Some(value)) = (
                parameters.get_mut(index),
                base_parameters.get(index).copied(),
            ) {
                *slot = value;
            }
        }
        for index in 0..3 {
            if self.angles.get(index).copied().unwrap_or(false) {
                continue;
            }
            if let (Some(slot), Some(value)) = (
                parameters.get_mut(3 + index),
                base_parameters.get(3 + index).copied(),
            ) {
                *slot = value;
            }
        }
        let corrected = lattice_of(&parameters)?;
        let corrected = if self.volume {
            corrected
        } else {
            // Volume = a*b*c*f(angles), so scaling the k free lengths by
            // s = (target / current)^(1/k) reaches the base volume exactly.
            let scaling: Vec<usize> = (0..3)
                .filter(|index| self.lengths.get(*index).copied().unwrap_or(false))
                .filter(|index| Some(*index) != keep_length)
                .collect();
            let current = corrected.volume().abs();
            let target = base.volume().abs();
            if scaling.is_empty() || current <= 0.0 || target <= 0.0 {
                return Err(Error::InvalidInput(
                    "holding the volume needs at least one free length".to_owned(),
                ));
            }
            let scale = (target / current).powf(1.0 / scaling.len() as f64);
            for index in scaling {
                if let Some(slot) = parameters.get_mut(index) {
                    *slot *= scale;
                }
            }
            lattice_of(&parameters)?
        };
        let base_standard = lattice_of(&base_parameters)?;
        // `corrected * base_standard^-1 * base` puts the corrected shape into
        // the orientation the base cell had.
        let map = mat3_mul(corrected.matrix(), base_standard.inverse());
        Ok(base.transform(map)?)
    }
}

/// Standard-orientation cell from six parameters.
fn lattice_of(parameters: &[f64; 6]) -> Result<Lattice, Error> {
    Ok(Lattice::from_parameters(
        parameters.first().copied().unwrap_or(1.0),
        parameters.get(1).copied().unwrap_or(1.0),
        parameters.get(2).copied().unwrap_or(1.0),
        parameters.get(3).copied().unwrap_or(90.0),
        parameters.get(4).copied().unwrap_or(90.0),
        parameters.get(5).copied().unwrap_or(90.0),
    )?)
}
