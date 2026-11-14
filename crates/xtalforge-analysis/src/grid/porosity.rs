//! Probe-occupiable volume, accessible surface area and the typed report
//! (Task 17.1.2).

use std::collections::BTreeMap;

use serde::{Deserialize, Serialize};
use xtalforge_core::{Error, Structure};

use super::{DistanceGrid, cross_length, distance_grid};

/// Cubic centimetres per cubic Angstrom.
const CM3_PER_A3: f64 = 1.0e-24;
/// Square metres per square Angstrom.
const M2_PER_A2: f64 = 1.0e-20;
/// Grams per atomic mass unit (CODATA 2018).
const GRAMS_PER_AMU: f64 = 1.660_539_066_60e-24;

/// Inputs for [`PorosityOptions::report`] and [`porosity_report`].
#[derive(Debug, Clone, Copy, PartialEq, Serialize, Deserialize)]
pub struct PorosityOptions {
    /// Probe radius in Angstrom.
    pub probe_radius: f64,
    /// Grid spacing in Angstrom.
    pub spacing: f64,
}

impl Default for PorosityOptions {
    fn default() -> Self {
        Self {
            probe_radius: 1.2,
            spacing: 0.5,
        }
    }
}

impl PorosityOptions {
    /// Builds the distance grid and reads every porosity quantity from it.
    ///
    /// Cell mass comes from [`Structure::cell_mass`], which expands symmetry
    /// to P1 and weights each atom by its occupancy before using the element
    /// table masses.
    ///
    /// # Errors
    ///
    /// Returns [`Error::InvalidInput`] for a non-positive or non-finite
    /// spacing, a negative or non-finite probe radius, a missing element
    /// radius, a grid above the sample limit, a failed P1 expansion, or a cell
    /// with no positive mass (for which the per-mass quantities are
    /// undefined).
    pub fn report(
        &self,
        structure: &Structure,
        radii: &BTreeMap<String, f64>,
    ) -> Result<PorosityReport, Error> {
        if !self.spacing.is_finite() || self.spacing <= 0.0 {
            return Err(Error::InvalidInput(format!(
                "grid spacing must be positive and finite, got {}",
                self.spacing
            )));
        }
        if !self.probe_radius.is_finite() || self.probe_radius < 0.0 {
            return Err(Error::InvalidInput(format!(
                "probe radius must be non-negative and finite, got {}",
                self.probe_radius
            )));
        }
        let grid = distance_grid(structure, radii, self.spacing)?;
        self.report_grid(structure, &grid)
    }

    /// Reads every porosity quantity from an existing distance grid.
    ///
    /// This avoids building the grid twice when the caller also needs the
    /// zero-level surface area or the largest included diameter. The grid must
    /// have been built from the same structure and with the same spacing.
    ///
    /// # Errors
    ///
    /// Returns the validation and mass errors of [`PorosityOptions::report`].
    pub fn report_grid(
        &self,
        structure: &Structure,
        grid: &DistanceGrid,
    ) -> Result<PorosityReport, Error> {
        if !self.probe_radius.is_finite() || self.probe_radius < 0.0 {
            return Err(Error::InvalidInput(format!(
                "probe radius must be non-negative and finite, got {}",
                self.probe_radius
            )));
        }
        let mass = structure.cell_mass()?;
        if !mass.is_finite() || mass <= 0.0 {
            return Err(Error::InvalidInput(format!(
                "the cell mass must be positive and finite for per-mass \
                 quantities, got {mass}"
            )));
        }
        let accessible_volume = grid.accessible_volume(self.probe_radius)?;
        let probe_occupiable_volume = accessible_volume;
        let surface_area = grid.accessible_surface_area(self.probe_radius)?;
        let percolation = grid.percolation(self.probe_radius)?;
        let pore_limiting_diameter = grid.pore_limiting_diameter()?;
        let grams = mass * GRAMS_PER_AMU;
        Ok(PorosityReport {
            probe_radius: self.probe_radius,
            spacing: self.spacing,
            accessible_volume,
            probe_occupiable_volume,
            pov_per_mass: probe_occupiable_volume * CM3_PER_A3 / grams,
            surface_area,
            asa_per_mass: surface_area * M2_PER_A2 / grams,
            percolating_directions: percolation.directions,
            pore_limiting_diameter,
        })
    }
}

/// Zeo++-style porosity quantities read from one distance grid.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct PorosityReport {
    /// Probe radius in Angstrom.
    pub probe_radius: f64,
    /// Grid spacing in Angstrom.
    pub spacing: f64,
    /// Volume where the probe fits, in cubic Angstrom per cell.
    ///
    /// Every sample with a surface distance of at least `probe_radius` is
    /// counted, whether its component reaches the cell boundary or not.
    pub accessible_volume: f64,
    /// Probe-occupiable volume in cubic Angstrom per cell.
    ///
    /// In this voxel model every sample where the probe fits can be occupied,
    /// so it equals [`PorosityReport::accessible_volume`]. Zeo++'s `POAV`
    /// additionally unions the probe spheres around the accessible samples and
    /// is not reproduced here; the field is kept so a report maps directly to
    /// the Zeo++ pore record.
    pub probe_occupiable_volume: f64,
    /// Probe-occupiable volume per cell mass, in cubic centimetres per gram.
    pub pov_per_mass: f64,
    /// Area of the probe isosurface in square Angstrom per cell.
    pub surface_area: f64,
    /// Accessible surface area per cell mass, in square metres per gram.
    pub asa_per_mass: f64,
    /// Whether the accessible void percolates along `a`, `b` and `c`.
    pub percolating_directions: [bool; 3],
    /// Pore-limiting diameter in Angstrom: twice the largest probe radius at
    /// which an accessible component spans at least one lattice direction.
    pub pore_limiting_diameter: f64,
}

/// Builds a distance grid and returns its Zeo++-style porosity quantities.
///
/// This is the four-argument form of [`PorosityOptions::report`]; it validates
/// `spacing` and `probe_radius` before touching the structure.
///
/// # Errors
///
/// Returns the errors of [`PorosityOptions::report`].
pub fn porosity_report(
    structure: &Structure,
    radii: &BTreeMap<String, f64>,
    spacing: f64,
    probe_radius: f64,
) -> Result<PorosityReport, Error> {
    PorosityOptions {
        probe_radius,
        spacing,
    }
    .report(structure, radii)
}

impl DistanceGrid {
    /// Area of the isosurface where the surface distance equals
    /// `probe_radius`, in square Angstrom.
    ///
    /// Grid faces whose two periodic ends are on opposite sides of the
    /// threshold are counted once and multiplied by the face area. A negative
    /// probe radius is treated as zero. This is the accessible surface area
    /// (ASA) of the probe, not the van der Waals area of
    /// [`DistanceGrid::surface_area`], which is the `0` level.
    ///
    /// # Errors
    ///
    /// Returns [`Error::InvalidInput`] when the probe radius is not finite.
    pub fn accessible_surface_area(&self, probe_radius: f64) -> Result<f64, Error> {
        if !probe_radius.is_finite() {
            return Err(Error::InvalidInput(format!(
                "probe radius must be finite, got {probe_radius}"
            )));
        }
        let threshold = probe_radius.max(0.0);
        let [n_x, n_y, n_z] = self.grid.shape;
        let step = self.grid.step;
        let area_x = cross_length(step[1], step[2]);
        let area_y = cross_length(step[0], step[2]);
        let area_z = cross_length(step[0], step[1]);
        let value = |x: usize, y: usize, z: usize| {
            self.grid
                .get(x % n_x, y % n_y, z % n_z)
                .unwrap_or(f64::NEG_INFINITY)
        };
        let mut area = 0.0;
        for z in 0..n_z {
            for y in 0..n_y {
                for x in 0..n_x {
                    let free = value(x, y, z) >= threshold;
                    if free != (value(x + 1, y, z) >= threshold) {
                        area += area_x;
                    }
                    if free != (value(x, y + 1, z) >= threshold) {
                        area += area_y;
                    }
                    if free != (value(x, y, z + 1) >= threshold) {
                        area += area_z;
                    }
                }
            }
        }
        Ok(area)
    }
}
