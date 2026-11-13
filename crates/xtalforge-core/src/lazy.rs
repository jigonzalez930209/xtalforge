//! Lazy symmetry expansion over the asymmetric unit (Task 15.2.4).
//!
//! Some operations only need to know how many atoms each asymmetric site
//! generates, which elements they are, and how much occupancy they carry.
//! [`LazyExpansion`] answers those questions by expanding one site at a time
//! with the same wrapping and coincidence rules as [`crate::P1Cell::expand`],
//! so a large cell never exists as a whole `Vec<P1Atom>`.
//!
//! # Where a lazy expansion cannot apply
//!
//! A result must be expressible as a sum over independent site orbits. That
//! covers the atom count, composition, occupancy sum, charge sum, per-site
//! multiplicity and orbit images. It does **not** cover anything that pairs
//! atoms from different sites or indexes the materialised cell:
//!
//! - bond perception, neighbor lists, rings, connectivity and UFF,
//! - PXRD and other structure factors that need every position at once
//!   (they may use the orbits, but the current implementations take a
//!   `Structure`),
//! - supercells, slabs and `Structure::as_p1`, which must return a cell,
//! - operations that address atoms by their final P1 index.
//!
//! Callers with those needs must materialise the P1 cell first.

use std::collections::BTreeMap;

use crate::p1::{P1Atom, site_atoms};
use crate::spatial::perpendicular_widths;
use crate::{Error, Structure, SymOp};

/// Summary of one asymmetric site's orbit.
#[derive(Debug, Clone, PartialEq)]
pub struct Orbit {
    /// Index of the parent asymmetric site.
    pub site_idx: usize,
    /// Number of distinct images within the expansion tolerance.
    pub multiplicity: usize,
    /// Element symbol of the parent site.
    pub element: String,
    /// Occupancy of every image.
    pub occupancy: f64,
    /// Disorder group inherited from the parent site, when present.
    pub disorder_group: Option<u32>,
}

/// Lazy, per-site view of a structure's P1 expansion.
///
/// Construction validates the tolerance and the cell's perpendicular widths;
/// every query expands one site at a time.
#[derive(Debug, Clone)]
pub struct LazyExpansion<'a> {
    structure: &'a Structure,
    widths: [f64; 3],
    tolerance: f64,
    operations: Vec<SymOp>,
}

impl<'a> LazyExpansion<'a> {
    /// Prepares a lazy expansion with an explicit coincidence tolerance.
    ///
    /// # Errors
    ///
    /// Returns [`Error::InvalidInput`] for a negative or non-finite tolerance,
    /// or the lattice error of the perpendicular-width calculation.
    pub fn new(structure: &'a Structure, tolerance: f64) -> Result<Self, Error> {
        if !tolerance.is_finite() || tolerance < 0.0 {
            return Err(Error::InvalidInput(
                "lazy expansion tolerance must be non-negative and finite".to_owned(),
            ));
        }
        let widths = perpendicular_widths(structure.lattice())?;
        let operations = if structure.operations().is_empty() {
            vec![SymOp::identity()]
        } else {
            structure.operations().to_vec()
        };
        Ok(Self {
            structure,
            widths,
            tolerance,
            operations,
        })
    }

    /// Structure this expansion reads.
    #[must_use]
    pub fn structure(&self) -> &'a Structure {
        self.structure
    }

    /// Tolerance this expansion uses.
    #[must_use]
    pub fn tolerance(&self) -> f64 {
        self.tolerance
    }

    /// Number of asymmetric sites.
    #[must_use]
    pub fn n_sites(&self) -> usize {
        self.structure.sites().len()
    }

    /// Orbit of one asymmetric site.
    ///
    /// # Errors
    ///
    /// Returns [`Error::InvalidInput`] when `site_idx` is out of range.
    pub fn orbit(&self, site_idx: usize) -> Result<Orbit, Error> {
        let site = self.structure.sites().get(site_idx).ok_or_else(|| {
            Error::InvalidInput(format!("asymmetric site {site_idx} is out of range"))
        })?;
        let images = self.site_images(site_idx)?;
        Ok(Orbit {
            site_idx,
            multiplicity: images.len(),
            element: site.element().to_owned(),
            occupancy: site.occupancy(),
            disorder_group: site.disorder_group(),
        })
    }

    /// Orbits of every asymmetric site, in site order.
    ///
    /// # Errors
    ///
    /// Returns the errors of [`LazyExpansion::orbit`].
    pub fn orbits(&self) -> Result<Vec<Orbit>, Error> {
        (0..self.n_sites()).map(|site| self.orbit(site)).collect()
    }

    /// Total number of P1 atoms.
    ///
    /// # Errors
    ///
    /// Returns the errors of [`LazyExpansion::orbits`].
    pub fn n_atoms(&self) -> Result<usize, Error> {
        Ok(self.orbits()?.iter().map(|orbit| orbit.multiplicity).sum())
    }

    /// Atom count per element in the P1 cell.
    ///
    /// Occupancies do not weight the count; use
    /// [`LazyExpansion::occupancy_sum`] for that.
    ///
    /// # Errors
    ///
    /// Returns the errors of [`LazyExpansion::orbits`].
    pub fn composition(&self) -> Result<BTreeMap<String, usize>, Error> {
        let mut counts = BTreeMap::new();
        for orbit in self.orbits()? {
            *counts.entry(orbit.element).or_insert(0) += orbit.multiplicity;
        }
        Ok(counts)
    }

    /// Sum of site occupancy over every generated atom.
    ///
    /// # Errors
    ///
    /// Returns the errors of [`LazyExpansion::orbits`].
    pub fn occupancy_sum(&self) -> Result<f64, Error> {
        Ok(self
            .orbits()?
            .iter()
            .map(|orbit| orbit.occupancy * orbit.multiplicity as f64)
            .sum())
    }

    /// Total charge over every generated atom, occupancy weighted like
    /// [`crate::Structure::total_charge`].
    ///
    /// # Errors
    ///
    /// Returns the errors of [`LazyExpansion::orbits`].
    pub fn charge_sum(&self) -> Result<f64, Error> {
        let mut total = 0.0;
        for orbit in self.orbits()? {
            let Some(site) = self.structure.sites().get(orbit.site_idx) else {
                continue;
            };
            if let Some(charge) = site.charge() {
                total += charge * orbit.occupancy * orbit.multiplicity as f64;
            }
        }
        Ok(total)
    }

    /// Images of one asymmetric site, without expanding any other site.
    ///
    /// # Errors
    ///
    /// Returns [`Error::InvalidInput`] when `site_idx` is out of range.
    pub fn site_images(&self, site_idx: usize) -> Result<Vec<P1Atom>, Error> {
        if self.structure.sites().get(site_idx).is_none() {
            return Err(Error::InvalidInput(format!(
                "asymmetric site {site_idx} is out of range"
            )));
        }
        Ok(site_atoms(
            self.structure,
            site_idx,
            &self.operations,
            self.widths,
            self.tolerance,
        ))
    }

    /// Visits every generated atom one orbit at a time.
    ///
    /// The visitor returns `false` to stop. The returned flag is `true` when
    /// every atom was visited and `false` when the visitor stopped early.
    ///
    /// # Errors
    ///
    /// Returns the errors of [`LazyExpansion::site_images`].
    pub fn for_each_atom(&self, mut visit: impl FnMut(&P1Atom) -> bool) -> Result<bool, Error> {
        for site_idx in 0..self.n_sites() {
            for atom in self.site_images(site_idx)? {
                if !visit(&atom) {
                    return Ok(false);
                }
            }
        }
        Ok(true)
    }
}
