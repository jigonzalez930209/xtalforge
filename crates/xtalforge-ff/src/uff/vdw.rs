//! UFF van der Waals term: Lennard-Jones 12-6 with a shifted cutoff
//! (Task 4.2.3).
//!
//! Ported from Crystal-Builder `xtal/ff/uff/terms.py`, eq 20:
//! `E = D [(x/r)^12 - 2 (x/r)^6]`. The pair energy at the cutoff is
//! subtracted so the potential is exactly zero there; without it a rebuilt
//! pair list would move the total by a step. Bonded 1-2 and geminal 1-3
//! pairs are excluded by the topology, not by this term.

use crate::Error;
use crate::uff::math::{EPS, add, merge_virial, negate, norm, scale, shift_vector, sub};
use crate::uff::parallel::{Part, evaluate};

/// One Lennard-Jones pair, with the image of `j`.
#[derive(Debug, Clone, Default)]
pub struct VanDerWaalsTerm {
    pub(crate) i: Vec<usize>,
    pub(crate) j: Vec<usize>,
    pub(crate) shift: Vec<[i32; 3]>,
    pub(crate) x: Vec<f64>,
    pub(crate) d: Vec<f64>,
    pub(crate) cutoff: f64,
    pub(crate) shift_energy: Vec<f64>,
}

impl VanDerWaalsTerm {
    /// Number of pairs.
    #[must_use]
    pub fn len(&self) -> usize {
        self.i.len()
    }

    /// Whether the term has no pairs.
    #[must_use]
    pub fn is_empty(&self) -> bool {
        self.i.is_empty()
    }

    /// Cutoff in Angstrom.
    #[must_use]
    pub fn cutoff(&self) -> f64 {
        self.cutoff
    }

    /// Energy in kcal/mol at these positions.
    ///
    /// # Errors
    ///
    /// Returns [`Error::Topology`] for a position list that does not match.
    pub fn energy(&self, positions: &[[f64; 3]], matrix: [[f64; 3]; 3]) -> Result<f64, Error> {
        Ok(self.energy_and_gradient(positions, matrix, None)?.0)
    }

    /// Energy and gradient `dE/dx` in kcal/mol and kcal/mol/Angstrom.
    ///
    /// # Errors
    ///
    /// Returns [`Error::Topology`] for a position list that does not match.
    pub fn energy_and_gradient(
        &self,
        positions: &[[f64; 3]],
        matrix: [[f64; 3]; 3],
        virial: Option<&mut [[f64; 3]; 3]>,
    ) -> Result<(f64, Vec<[f64; 3]>), Error> {
        self.energy_and_gradient_with(positions, matrix, virial, true)
    }

    /// Energy and gradient, optionally evaluated without worker threads.
    ///
    /// # Errors
    ///
    /// Returns [`Error::Topology`] for a position list that does not match.
    pub(crate) fn energy_and_gradient_with(
        &self,
        positions: &[[f64; 3]],
        matrix: [[f64; 3]; 3],
        virial: Option<&mut [[f64; 3]; 3]>,
        parallel: bool,
    ) -> Result<(f64, Vec<[f64; 3]>), Error> {
        let mut gradient = vec![[0.0; 3]; positions.len()];
        let mut energy = 0.0;
        let mut virial_local = [[0.0; 3]; 3];
        evaluate(
            self.len(),
            parallel,
            |index| self.part(index, positions, matrix),
            &mut energy,
            &mut gradient,
            &mut virial_local,
        )?;
        if let Some(target) = virial {
            merge_virial(target, virial_local);
        }
        Ok((energy, gradient))
    }

    fn part(
        &self,
        index: usize,
        positions: &[[f64; 3]],
        matrix: [[f64; 3]; 3],
    ) -> Result<Part, Error> {
        let mut part = Part::zero();
        let (Some(atom_i), Some(atom_j), Some(shift), Some(x), Some(depth)) = (
            self.i.get(index),
            self.j.get(index),
            self.shift.get(index),
            self.x.get(index),
            self.d.get(index),
        ) else {
            return Err(Error::Topology(
                "van der Waals arrays have inconsistent lengths".to_owned(),
            ));
        };
        let start = positions
            .get(*atom_i)
            .ok_or_else(|| Error::Topology("van der Waals atom out of range".to_owned()))?;
        let end = positions
            .get(*atom_j)
            .ok_or_else(|| Error::Topology("van der Waals atom out of range".to_owned()))?;
        let displacement = sub(add(*end, shift_vector(*shift, matrix)), *start);
        let distance = norm(displacement);
        if distance <= EPS || distance >= self.cutoff {
            return Ok(part);
        }
        let ratio6 = (x / distance).powi(6);
        let ratio12 = ratio6 * ratio6;
        let shift_energy = self.shift_energy.get(index).copied().unwrap_or(0.0);
        part.energy = depth * (ratio12 - 2.0 * ratio6) - shift_energy;
        let factor = 12.0 * depth * (ratio6 - ratio12) / (distance * distance);
        let gradient_vector = scale(displacement, factor);
        part.push_virial(displacement, gradient_vector);
        part.push_scatter(*atom_j, gradient_vector);
        part.push_scatter(*atom_i, negate(gradient_vector));
        Ok(part)
    }
}
