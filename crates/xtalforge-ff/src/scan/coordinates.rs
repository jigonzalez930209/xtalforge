//! Scan coordinates: internal distances, angles, torsions, cell parameters and
//! collective variables (Tasks 4.5.1 and 14.2.3).
//!
//! [`ScanCoordinate`] names the quantity a scan steps through:
//!
//! * [`ScanCoordinate::Distance`], [`ScanCoordinate::Angle`] and
//!   [`ScanCoordinate::Torsion`] are internal geometric coordinates addressed
//!   by site index; values are Angstrom or degrees.
//! * [`ScanCoordinate::CellParameter`] varies `a`, `b`, `c`, `alpha`, `beta` or
//!   `gamma` of the cell (Angstrom or degrees).
//! * [`ScanCoordinate::Collective`] wraps a [`CollectiveVariable`], a coordinate
//!   over atom groups such as a centroid distance or a linear combination.
//!
//! Scans operate on P1 structures, so the indices are P1 atom indices.  Values
//! and gradients are defined for every kind; a cell parameter has a zero
//! Cartesian gradient because it is imposed on the lattice rather than relaxed
//! through the sites, which is why it must be moved with
//! [`ScanCoordinate::set_in`] and a [`ScanGeometry`].

use serde::{Deserialize, Serialize};
use xtalforge_core::{Lattice, Structure};

use crate::Error;
use crate::structure_util::with_frac;

use super::cell_parameter::CellParameter;
use super::collective::CollectiveVariable;
use super::internal;

/// A mutable lattice and Cartesian positions that a scan coordinate acts on.
#[derive(Debug, Clone, PartialEq)]
pub struct ScanGeometry {
    /// Cell of the geometry.
    pub lattice: Lattice,
    /// Cartesian positions, one per site.
    pub positions: Vec<[f64; 3]>,
}

impl ScanGeometry {
    /// Reads the geometry of one structure.
    #[must_use]
    pub fn from_structure(structure: &Structure) -> Self {
        Self {
            lattice: *structure.lattice(),
            positions: structure.cart_sites(),
        }
    }

    /// Rebuilds `base` with this geometry's lattice and positions.
    ///
    /// # Errors
    ///
    /// Returns [`Error::InvalidInput`] when the position count differs, or a
    /// core error when the structure cannot be rebuilt.
    pub fn into_structure(self, base: &Structure) -> Result<Structure, Error> {
        if self.positions.len() != base.sites().len() {
            return Err(Error::InvalidInput(format!(
                "scan geometry has {} positions for a structure with {} sites",
                self.positions.len(),
                base.sites().len()
            )));
        }
        let frac: Vec<[f64; 3]> = self
            .positions
            .iter()
            .map(|position| self.lattice.to_frac(*position))
            .collect();
        with_frac(base, self.lattice, &frac)
    }
}

/// One coordinate a scan can step through.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
#[serde(tag = "kind", rename_all = "snake_case", deny_unknown_fields)]
pub enum ScanCoordinate {
    /// Distance between sites `i` and `j`.
    Distance {
        /// First site index.
        i: usize,
        /// Second site index.
        j: usize,
    },
    /// Angle at site `j` between `i` and `k`, in degrees.
    Angle {
        /// First endpoint.
        i: usize,
        /// Vertex.
        j: usize,
        /// Second endpoint.
        k: usize,
    },
    /// Torsion across `i-j-k-l`, in degrees in `(-180, 180]`.
    Torsion {
        /// First atom.
        i: usize,
        /// Second atom.
        j: usize,
        /// Third atom.
        k: usize,
        /// Fourth atom.
        l: usize,
    },
    /// One cell parameter, imposed on the lattice.
    CellParameter {
        /// Parameter being scanned.
        parameter: CellParameter,
    },
    /// A collective variable over atom groups.
    Collective {
        /// Variable definition.
        variable: CollectiveVariable,
    },
}

impl ScanCoordinate {
    /// Distance coordinate.
    ///
    /// # Errors
    ///
    /// Returns [`Error::InvalidInput`] when the two indices are equal.
    pub fn distance(i: usize, j: usize) -> Result<Self, Error> {
        if i == j {
            return Err(Error::InvalidInput(
                "a distance needs two different atoms".to_owned(),
            ));
        }
        Ok(Self::Distance { i, j })
    }

    /// Angle coordinate, vertex in the middle.
    ///
    /// # Errors
    ///
    /// Returns [`Error::InvalidInput`] when two indices are equal.
    pub fn angle(i: usize, j: usize, k: usize) -> Result<Self, Error> {
        if i == j || j == k || i == k {
            return Err(Error::InvalidInput(
                "an angle needs three different atoms".to_owned(),
            ));
        }
        Ok(Self::Angle { i, j, k })
    }

    /// Torsion coordinate.
    ///
    /// # Errors
    ///
    /// Returns [`Error::InvalidInput`] when two indices are equal.
    pub fn torsion(i: usize, j: usize, k: usize, l: usize) -> Result<Self, Error> {
        let indices = [i, j, k, l];
        for (first, value) in indices.iter().enumerate() {
            for other in indices.iter().skip(first + 1) {
                if value == other {
                    return Err(Error::InvalidInput(
                        "a torsion needs four different atoms".to_owned(),
                    ));
                }
            }
        }
        Ok(Self::Torsion { i, j, k, l })
    }

    /// Cell-parameter coordinate.
    #[must_use]
    pub fn cell_parameter(parameter: CellParameter) -> Self {
        Self::CellParameter { parameter }
    }

    /// Collective-variable coordinate.
    #[must_use]
    pub fn collective(variable: CollectiveVariable) -> Self {
        Self::Collective { variable }
    }

    /// The cell parameter, when this is one.
    #[must_use]
    pub fn parameter(&self) -> Option<CellParameter> {
        match self {
            Self::CellParameter { parameter } => Some(*parameter),
            _ => None,
        }
    }

    /// Whether the coordinate moves the cell instead of the sites.
    #[must_use]
    pub fn is_cell_parameter(&self) -> bool {
        matches!(self, Self::CellParameter { .. })
    }

    /// Human-readable label.
    #[must_use]
    pub fn label(&self) -> String {
        match self {
            Self::Distance { i, j } => format!("distance {i}-{j}"),
            Self::Angle { i, j, k } => format!("angle {i}-{j}-{k}"),
            Self::Torsion { i, j, k, l } => format!("torsion {i}-{j}-{k}-{l}"),
            Self::CellParameter { parameter } => format!("cell {}", parameter.name()),
            Self::Collective { variable } => variable.label(),
        }
    }

    /// Units of the value: `"A"` or `"deg"`.
    #[must_use]
    pub fn units(&self) -> &'static str {
        match self {
            Self::Distance { .. } | Self::Collective { .. } => "A",
            Self::Angle { .. } | Self::Torsion { .. } => "deg",
            Self::CellParameter { parameter } => parameter.units(),
        }
    }

    /// Site indices this coordinate reads; empty for a cell parameter.
    #[must_use]
    pub fn indices(&self) -> Vec<usize> {
        match self {
            Self::Distance { i, j } => vec![*i, *j],
            Self::Angle { i, j, k } => vec![*i, *j, *k],
            Self::Torsion { i, j, k, l } => vec![*i, *j, *k, *l],
            Self::CellParameter { .. } => Vec::new(),
            Self::Collective { variable } => variable.indices(),
        }
    }

    /// Checks the coordinate against a site count.
    ///
    /// # Errors
    ///
    /// Returns [`Error::InvalidInput`] for an index outside `0..n_sites` or a
    /// malformed collective variable.
    pub fn validate(&self, n_sites: usize) -> Result<(), Error> {
        match self {
            Self::CellParameter { .. } => Ok(()),
            Self::Collective { variable } => variable.validate(n_sites),
            _ => {
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
        }
    }

    /// Coordinate value for Cartesian site positions and a lattice.
    ///
    /// For a cell parameter the value is read from `lattice` and `positions`
    /// is ignored.
    ///
    /// # Errors
    ///
    /// Returns [`Error::InvalidInput`] for an index outside the position
    /// array, a degenerate geometry, or a non-finite result.
    pub fn value(&self, positions: &[[f64; 3]], lattice: &Lattice) -> Result<f64, Error> {
        match self {
            Self::Distance { i, j } => internal::distance_value(positions, lattice, *i, *j),
            Self::Angle { i, j, k } => internal::angle_value(positions, lattice, *i, *j, *k),
            Self::Torsion { i, j, k, l } => {
                internal::torsion_value(positions, lattice, *i, *j, *k, *l)
            }
            Self::CellParameter { parameter } => Ok(parameter.value(lattice)),
            Self::Collective { variable } => variable.value(positions, lattice),
        }
    }

    /// Cartesian gradient `d(value)/d(position)` for every site.
    ///
    /// A cell parameter has a zero gradient, because the sites do not move it.
    ///
    /// # Errors
    ///
    /// Returns the same errors as [`ScanCoordinate::value`].
    pub fn gradient(
        &self,
        positions: &[[f64; 3]],
        lattice: &Lattice,
    ) -> Result<Vec<[f64; 3]>, Error> {
        match self {
            Self::Distance { i, j } => internal::distance_gradient(positions, lattice, *i, *j),
            Self::Angle { i, j, k } => internal::angle_gradient(positions, lattice, *i, *j, *k),
            Self::Torsion { i, j, k, l } => {
                internal::torsion_gradient(positions, lattice, *i, *j, *k, *l)
            }
            Self::CellParameter { .. } => Ok(vec![[0.0; 3]; positions.len()]),
            Self::Collective { variable } => variable.gradient(positions, lattice),
        }
    }

    /// Value on a full geometry.
    ///
    /// # Errors
    ///
    /// Returns the errors of [`ScanCoordinate::value`].
    pub fn value_in(&self, geometry: &ScanGeometry) -> Result<f64, Error> {
        self.value(&geometry.positions, &geometry.lattice)
    }

    /// Moves the geometry so the coordinate reads `target`.
    ///
    /// This is a seed, not a constraint: internal coordinates move one atom or
    /// group, a collective variable translates whole groups and a cell
    /// parameter rebuilds the lattice while keeping every fractional
    /// coordinate fixed, so the atoms follow the cell.  The optimizer relaxes
    /// the rest.
    ///
    /// # Errors
    ///
    /// Returns [`Error::InvalidInput`] for a degenerate starting geometry, a
    /// non-finite target or an impossible cell.
    pub fn set_in(&self, geometry: &mut ScanGeometry, target: f64) -> Result<(), Error> {
        match self {
            Self::CellParameter { parameter } => {
                let frac: Vec<[f64; 3]> = geometry
                    .positions
                    .iter()
                    .map(|position| geometry.lattice.to_frac(*position))
                    .collect();
                let lattice = parameter.with_value(&geometry.lattice, target)?;
                geometry.positions = frac.iter().map(|point| lattice.to_cart(*point)).collect();
                geometry.lattice = lattice;
                Ok(())
            }
            _ => self.set(&mut geometry.positions, &geometry.lattice, target),
        }
    }

    /// Moves Cartesian positions so the coordinate reads `target`.
    ///
    /// # Errors
    ///
    /// Returns [`Error::InvalidInput`] for a degenerate starting geometry, a
    /// non-finite target, or a cell parameter, which needs
    /// [`ScanCoordinate::set_in`] to rebuild the lattice.
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
            Self::Distance { i, j } => internal::set_distance(positions, lattice, *i, *j, target),
            Self::Angle { i, j, k } => internal::set_angle(positions, lattice, *i, *j, *k, target),
            Self::Torsion { i, j, k, l } => {
                internal::set_torsion(positions, lattice, *i, *j, *k, *l, target)
            }
            Self::CellParameter { parameter } => Err(Error::InvalidInput(format!(
                "cell parameter {} is part of the lattice; use ScanCoordinate::set_in",
                parameter.name()
            ))),
            Self::Collective { variable } => variable.set(positions, lattice, target),
        }
    }
}
