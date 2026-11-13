//! Collective scan coordinates (Task 14.2.3).
//!
//! Two general-purpose variables complement the internal distance, angle and
//! torsion coordinates:
//!
//! * [`CollectiveVariable::CentroidDistance`] is the distance between the
//!   centroids of two atom groups.  With one atom per group it is an ordinary
//!   distance; with larger groups it describes a whole fragment moving away
//!   from another.
//! * [`CollectiveVariable::Linear`] is `sum_i w_i (r_i . axis)`, the linear
//!   combination used for proton transfer, layer separation and other
//!   reaction coordinates that are linear in the Cartesian positions.
//!
//! Both values are Angstrom, both gradients are analytic, and
//! [`CollectiveVariable::set`] moves the geometry to a target by translating
//! whole groups, which is the same "seed, then relax" contract as the internal
//! coordinates.  Centroid distances use the minimum image; the linear
//! combination reads the raw Cartesian positions, so its caller should work on
//! a whole (non-wrapped) fragment.

use serde::{Deserialize, Serialize};
use xtalforge_core::Lattice;

use crate::Error;
use crate::math::{add, dot, scale, sub, unit};

/// A coordinate defined over atom groups rather than one bonded path.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
#[serde(tag = "kind", rename_all = "snake_case", deny_unknown_fields)]
pub enum CollectiveVariable {
    /// Distance between the centroids of two atom groups.
    CentroidDistance {
        /// Atoms of the first group.
        first: Vec<usize>,
        /// Atoms of the second group.
        second: Vec<usize>,
    },
    /// `sum_i weights_i * (position_i . axis)`, a linear combination.
    Linear {
        /// Atoms that enter the sum.
        atoms: Vec<usize>,
        /// Coefficient of each atom, same length as `atoms`.
        weights: Vec<f64>,
        /// Direction the combination projects on; not normalised by the crate.
        axis: [f64; 3],
    },
}

impl CollectiveVariable {
    /// Centroid distance between two non-empty, disjoint groups.
    ///
    /// # Errors
    ///
    /// Returns [`Error::InvalidInput`] for an empty group or an atom that
    /// appears in both groups.
    pub fn centroid_distance(first: Vec<usize>, second: Vec<usize>) -> Result<Self, Error> {
        if first.is_empty() || second.is_empty() {
            return Err(Error::InvalidInput(
                "a centroid distance needs two non-empty atom groups".to_owned(),
            ));
        }
        if first.iter().any(|index| second.contains(index)) {
            return Err(Error::InvalidInput(
                "the two groups of a centroid distance must not share an atom".to_owned(),
            ));
        }
        Ok(Self::CentroidDistance { first, second })
    }

    /// Linear combination of Cartesian positions along `axis`.
    ///
    /// # Errors
    ///
    /// Returns [`Error::InvalidInput`] for an empty or mismatched atom/weight
    /// list, non-finite numbers, a zero axis, or weights whose squares sum to
    /// zero (which no translation could move).
    pub fn linear(atoms: Vec<usize>, weights: Vec<f64>, axis: [f64; 3]) -> Result<Self, Error> {
        if atoms.is_empty() {
            return Err(Error::InvalidInput(
                "a linear collective variable needs at least one atom".to_owned(),
            ));
        }
        if atoms.len() != weights.len() {
            return Err(Error::InvalidInput(format!(
                "a linear collective variable needs one weight per atom, got {} atoms and {} weights",
                atoms.len(),
                weights.len()
            )));
        }
        if weights.iter().any(|weight| !weight.is_finite()) {
            return Err(Error::InvalidInput(
                "collective-variable weights must be finite".to_owned(),
            ));
        }
        if weights.iter().map(|weight| weight * weight).sum::<f64>() <= 1e-30 {
            return Err(Error::InvalidInput(
                "collective-variable weights are all zero".to_owned(),
            ));
        }
        if !axis.iter().all(|value| value.is_finite()) {
            return Err(Error::InvalidInput(
                "collective-variable axis must be finite".to_owned(),
            ));
        }
        if dot(axis, axis) <= 1e-30 {
            return Err(Error::InvalidInput(
                "collective-variable axis has zero length".to_owned(),
            ));
        }
        Ok(Self::Linear {
            atoms,
            weights,
            axis,
        })
    }

    /// Human-readable label.
    #[must_use]
    pub fn label(&self) -> String {
        match self {
            Self::CentroidDistance { first, second } => {
                format!("centroid distance {first:?}-{second:?}")
            }
            Self::Linear { atoms, .. } => format!("linear {atoms:?}"),
        }
    }

    /// Units of the value; always `"A"`.
    #[must_use]
    pub fn units(&self) -> &'static str {
        "A"
    }

    /// Atom indices the variable reads.
    #[must_use]
    pub fn indices(&self) -> Vec<usize> {
        match self {
            Self::CentroidDistance { first, second } => {
                let mut indices = first.clone();
                indices.extend(second.iter().copied());
                indices
            }
            Self::Linear { atoms, .. } => atoms.clone(),
        }
    }

    /// Checks the indices against a site count.
    ///
    /// # Errors
    ///
    /// Returns [`Error::InvalidInput`] for an index outside `0..n_sites`.
    pub fn validate(&self, n_sites: usize) -> Result<(), Error> {
        for index in self.indices() {
            if index >= n_sites {
                return Err(Error::InvalidInput(format!(
                    "{} refers to atom {index}, but the structure has {n_sites} sites",
                    self.label()
                )));
            }
        }
        Ok(())
    }

    /// Value in Angstrom.
    ///
    /// # Errors
    ///
    /// Returns [`Error::InvalidInput`] for an out-of-range index, coincident
    /// centroids or a non-finite result.
    pub fn value(&self, positions: &[[f64; 3]], lattice: &Lattice) -> Result<f64, Error> {
        match self {
            Self::CentroidDistance { first, second } => {
                let difference = centroid_difference(positions, lattice, first, second)?;
                let distance = crate::math::norm(difference);
                if distance <= 1e-12 {
                    return Err(Error::InvalidInput(
                        "the two centroids of a centroid distance coincide".to_owned(),
                    ));
                }
                Ok(distance)
            }
            Self::Linear {
                atoms,
                weights,
                axis,
            } => {
                let mut total = 0.0;
                for (index, weight) in atoms.iter().zip(weights) {
                    let position = positions.get(*index).ok_or_else(|| {
                        Error::InvalidInput(format!(
                            "atom index {index} is outside the position list"
                        ))
                    })?;
                    total += weight * dot(*position, *axis);
                }
                if !total.is_finite() {
                    return Err(Error::InvalidInput(
                        "the collective variable is not finite".to_owned(),
                    ));
                }
                Ok(total)
            }
        }
    }

    /// Cartesian gradient `d(value)/d(position)` for every site.
    ///
    /// # Errors
    ///
    /// Returns the same errors as [`CollectiveVariable::value`].
    pub fn gradient(
        &self,
        positions: &[[f64; 3]],
        lattice: &Lattice,
    ) -> Result<Vec<[f64; 3]>, Error> {
        let mut result = vec![[0.0; 3]; positions.len()];
        match self {
            Self::CentroidDistance { first, second } => {
                let difference = centroid_difference(positions, lattice, first, second)?;
                let Some(direction) = unit(difference) else {
                    return Err(Error::InvalidInput(
                        "the two centroids of a centroid distance coincide".to_owned(),
                    ));
                };
                let first_count = first.len() as f64;
                let second_count = second.len() as f64;
                let first_part = scale(direction, -1.0 / first_count);
                let second_part = scale(direction, 1.0 / second_count);
                scatter(&mut result, first, first_part);
                scatter(&mut result, second, second_part);
            }
            Self::Linear {
                atoms,
                weights,
                axis,
            } => {
                let direction = *axis;
                for (index, weight) in atoms.iter().zip(weights) {
                    if let Some(slot) = result.get_mut(*index) {
                        *slot = add(*slot, scale(direction, *weight));
                    }
                }
            }
        }
        Ok(result)
    }

    /// Moves the geometry so the variable reads `target`.
    ///
    /// A centroid distance translates the second group; a linear combination
    /// distributes the change over the atoms in proportion to their weights,
    /// which is the minimum-norm translation for that coordinate.
    ///
    /// # Errors
    ///
    /// Returns [`Error::InvalidInput`] for a degenerate starting geometry or a
    /// non-finite target.
    pub fn set(
        &self,
        positions: &mut [[f64; 3]],
        lattice: &Lattice,
        target: f64,
    ) -> Result<(), Error> {
        if !target.is_finite() {
            return Err(Error::InvalidInput("scan target must be finite".to_owned()));
        }
        match self {
            Self::CentroidDistance { first, second } => {
                let difference = centroid_difference(positions, lattice, first, second)?;
                let Some(direction) = unit(difference) else {
                    return Err(Error::InvalidInput(
                        "the two centroids of a centroid distance coincide".to_owned(),
                    ));
                };
                let change = target - crate::math::norm(difference);
                // Moving every atom of the second group by the same vector
                // moves its centroid by that vector.
                let shift = scale(direction, change);
                scatter(&mut *positions, second, shift);
            }
            Self::Linear {
                atoms,
                weights,
                axis,
            } => {
                let current = self.value(positions, lattice)?;
                let denominator: f64 = weights.iter().map(|weight| weight * weight).sum();
                if denominator <= 1e-30 {
                    return Err(Error::InvalidInput(
                        "collective-variable weights are all zero".to_owned(),
                    ));
                }
                for (index, weight) in atoms.iter().zip(weights) {
                    if let Some(position) = positions.get_mut(*index) {
                        *position = add(
                            *position,
                            scale(*axis, weight * (target - current) / denominator),
                        );
                    }
                }
            }
        }
        Ok(())
    }
}

/// Minimum-image vector from the first centroid to the second.
fn centroid_difference(
    positions: &[[f64; 3]],
    lattice: &Lattice,
    first: &[usize],
    second: &[usize],
) -> Result<[f64; 3], Error> {
    let first_centroid = centroid(positions, first)?;
    let second_centroid = centroid(positions, second)?;
    let difference = sub(second_centroid, first_centroid);
    let frac = lattice.to_frac(difference);
    let (wrapped, _) = lattice.minimum_image(frac);
    Ok(lattice.to_cart(wrapped))
}

/// Arithmetic mean of a group's positions.
fn centroid(positions: &[[f64; 3]], group: &[usize]) -> Result<[f64; 3], Error> {
    if group.is_empty() {
        return Err(Error::InvalidInput(
            "a centroid needs at least one atom".to_owned(),
        ));
    }
    let mut total = [0.0; 3];
    for index in group {
        let position = positions.get(*index).ok_or_else(|| {
            Error::InvalidInput(format!("atom index {index} is outside the position list"))
        })?;
        total = add(total, *position);
    }
    Ok(scale(total, 1.0 / group.len() as f64))
}

/// Adds one contribution to every atom of a group.
fn scatter(gradient: &mut [[f64; 3]], group: &[usize], value: [f64; 3]) {
    for index in group {
        if let Some(slot) = gradient.get_mut(*index) {
            *slot = add(*slot, value);
        }
    }
}
