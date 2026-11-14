//! Periodic distance grids and the porosity numbers read off them (Tasks 5.2,
//! 5.3 and 17.1.2).
//!
//! The field stored in [`DistanceGrid`] is the distance from every sample to
//! the *surface* of the nearest atom, `centre distance - radius`, with the
//! minimum-image convention. It is negative inside an atom, zero on the van
//! der Waals surface, and its `r` level is where the centre of a probe of
//! radius `r` can sit.
//!
//! The grid itself lives in this module. The Zeo++-style extras live beside it:
//! [`percolation`](DistanceGrid::percolation) and
//! [`percolation_radius`](DistanceGrid::percolation_radius) in the private
//! `percolation` submodule, [`accessible_surface_area`](DistanceGrid::accessible_surface_area)
//! and the typed [`PorosityReport`] in the private `porosity` submodule. Their
//! public types and the [`porosity_report`] entry point are re-exported here,
//! so the public path remains `xtalforge_analysis::grid::...`.
//!
//! # Approximations (deliberate, documented here)
//!
//! * Every grid sample is compared with every atom. This is exact but
//!   `O(points * atoms)`; it is intended for the small and medium cells of a
//!   test suite, not for a 5000-atom framework.
//! * The accessible volume is a voxel count times the volume of one sample.
//!   It therefore converges from below as the spacing shrinks and can differ
//!   from a continuous integral by a boundary-layer error of about one
//!   spacing.
//! * The surface area counts grid faces whose two ends have opposite sign
//!   (outside versus inside) and multiplies by the face area. This is the
//!   area of the voxelised isosurface and overestimates a smooth surface at
//!   coarse spacing.
//! * The pore-size histogram bins `2 * distance` (the diameter of the largest
//!   probe that fits at that sample) with a caller-supplied bin width. It is
//!   a geometric pore-size distribution, not a Zeo++ result.
//!
//! All loops are ordered, so the output is bit-identical between runs.

use std::collections::BTreeMap;

use serde::{Deserialize, Serialize};
use xtalforge_core::elements::element_by_symbol;
use xtalforge_core::{Error, Structure};
use xtalforge_io::grid::Grid3;

mod percolation;
mod porosity;

pub use percolation::Percolation;
pub use porosity::{PorosityOptions, PorosityReport, porosity_report};

/// Grid samples above which a distance grid is refused.
const MAX_GRID_POINTS: usize = 8_000_000;
/// Largest number of histogram bins accepted by [`DistanceGrid::pore_diameter`].
const MAX_HISTOGRAM_BINS: usize = 1_000_000;

/// A periodic field of distances to the nearest atom surface.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct DistanceGrid {
    /// Sampled scalar field in Angstrom, one value per grid point.
    pub grid: Grid3,
    /// Radii in Angstrom used for each element symbol.
    pub radii: BTreeMap<String, f64>,
}

/// One bin of the geometric pore-size histogram.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct PoreSizeBin {
    /// Lower diameter bound in Angstrom (inclusive).
    pub diameter_min: f64,
    /// Upper diameter bound in Angstrom (exclusive, except the last bin).
    pub diameter_max: f64,
    /// Number of grid samples in the bin.
    pub count: usize,
    /// Volume represented by the samples in the bin, in cubic Angstrom.
    pub volume: f64,
}

/// Default van der Waals radii from the element table.
///
/// # Errors
///
/// Returns [`Error::InvalidInput`] when a site element is missing from the
/// reference table.
pub fn default_radii(structure: &Structure) -> Result<BTreeMap<String, f64>, Error> {
    let mut radii = BTreeMap::new();
    for site in structure.sites() {
        let element = element_by_symbol(site.element())?;
        radii.insert(element.symbol.clone(), element.vdw_radius);
    }
    Ok(radii)
}

/// Samples the distance-to-nearest-atom-surface field of a structure.
///
/// The cell is sampled at approximately `spacing` Angstrom, using the cell
/// lengths so that anisotropic cells get a similar resolution along every
/// axis. `radii` must contain a positive radius for every element present in
/// the P1 expansion.
///
/// # Errors
///
/// Returns [`Error::InvalidInput`] for a non-positive or non-finite spacing, a
/// missing or invalid radius, a grid above the size limit, or a structure
/// whose P1 expansion fails.
pub fn distance_grid(
    structure: &Structure,
    radii: &BTreeMap<String, f64>,
    spacing: f64,
) -> Result<DistanceGrid, Error> {
    if !spacing.is_finite() || spacing <= 0.0 {
        return Err(Error::InvalidInput(format!(
            "grid spacing must be positive and finite, got {spacing}"
        )));
    }
    let lattice = structure.lattice();
    let lengths = lattice.lengths();
    let shape = lengths.map(|length| {
        let count = (length / spacing).round() as usize;
        count.max(1)
    });
    let total = shape[0]
        .checked_mul(shape[1])
        .and_then(|value| value.checked_mul(shape[2]))
        .ok_or_else(|| Error::InvalidInput("grid dimensions overflow".to_owned()))?;
    if total > MAX_GRID_POINTS {
        return Err(Error::InvalidInput(format!(
            "a {}x{}x{} grid has {total} samples; use a larger spacing",
            shape[0], shape[1], shape[2]
        )));
    }

    let cell = structure.p1()?;
    let mut atoms: Vec<([f64; 3], f64)> = Vec::with_capacity(cell.atoms.len());
    for atom in &cell.atoms {
        if atom.occupancy <= 0.0 {
            continue;
        }
        let radius = radii.get(&atom.element).copied().ok_or_else(|| {
            Error::InvalidInput(format!("no radius supplied for element {:?}", atom.element))
        })?;
        if !radius.is_finite() || radius <= 0.0 {
            return Err(Error::InvalidInput(format!(
                "radius for element {:?} must be positive and finite, got {radius}",
                atom.element
            )));
        }
        atoms.push((atom.frac, radius));
    }

    let matrix = lattice.matrix();
    let [n_x, n_y, n_z] = shape;
    let step = [
        [
            matrix[0][0] / n_x as f64,
            matrix[0][1] / n_x as f64,
            matrix[0][2] / n_x as f64,
        ],
        [
            matrix[1][0] / n_y as f64,
            matrix[1][1] / n_y as f64,
            matrix[1][2] / n_y as f64,
        ],
        [
            matrix[2][0] / n_z as f64,
            matrix[2][1] / n_z as f64,
            matrix[2][2] / n_z as f64,
        ],
    ];

    let mut data = Vec::with_capacity(total);
    for z in 0..n_z {
        for y in 0..n_y {
            for x in 0..n_x {
                let point = [
                    x as f64 / n_x as f64,
                    y as f64 / n_y as f64,
                    z as f64 / n_z as f64,
                ];
                let mut nearest = f64::INFINITY;
                for (frac, radius) in &atoms {
                    let delta = [point[0] - frac[0], point[1] - frac[1], point[2] - frac[2]];
                    let distance = lattice.minimum_image(delta).1 - radius;
                    if distance < nearest {
                        nearest = distance;
                    }
                }
                data.push(nearest);
            }
        }
    }

    let grid = Grid3::new([0.0; 3], step, shape, data, true)?;
    Ok(DistanceGrid {
        grid,
        radii: radii.clone(),
    })
}

impl DistanceGrid {
    /// Number of samples.
    #[must_use]
    pub fn len(&self) -> usize {
        self.grid.len()
    }

    /// Whether the grid has no samples (never true for a valid grid).
    #[must_use]
    pub fn is_empty(&self) -> bool {
        self.grid.is_empty()
    }

    /// Volume represented by one sample in cubic Angstrom.
    #[must_use]
    pub fn point_volume(&self) -> f64 {
        self.cell_volume() / self.len() as f64
    }

    /// Volume of the periodic cell in cubic Angstrom.
    #[must_use]
    pub fn cell_volume(&self) -> f64 {
        let step = self.grid.step;
        let determinant = step[0][0] * (step[1][1] * step[2][2] - step[1][2] * step[2][1])
            - step[0][1] * (step[1][0] * step[2][2] - step[1][2] * step[2][0])
            + step[0][2] * (step[1][0] * step[2][1] - step[1][1] * step[2][0]);
        let [n_x, n_y, n_z] = self.grid.shape;
        determinant.abs() * n_x as f64 * n_y as f64 * n_z as f64
    }

    /// Effective spacing in Angstrom, `(cell volume / samples)^(1/3)`.
    #[must_use]
    pub fn spacing(&self) -> f64 {
        self.point_volume().cbrt()
    }

    fn value(&self, x: usize, y: usize, z: usize) -> f64 {
        self.grid.get(x, y, z).unwrap_or(f64::INFINITY)
    }

    /// Cartesian positions of the samples closest to an atom surface.
    ///
    /// A sample is included when `|distance| <= spacing / 2`.
    #[must_use]
    pub fn contact_surface(&self) -> Vec<[f64; 3]> {
        let [n_x, n_y, n_z] = self.grid.shape;
        let threshold = 0.5 * self.spacing();
        let mut surface = Vec::new();
        for z in 0..n_z {
            for y in 0..n_y {
                for x in 0..n_x {
                    if self.value(x, y, z).abs() <= threshold {
                        surface.push(self.grid.position(x, y, z));
                    }
                }
            }
        }
        surface
    }

    /// Volume where a probe of `probe_radius` fits, in cubic Angstrom.
    ///
    /// A sample is included when its surface distance is at least
    /// `probe_radius`. A negative probe radius is treated as zero.
    ///
    /// # Errors
    ///
    /// Returns [`Error::InvalidInput`] when the probe radius is not finite.
    pub fn accessible_volume(&self, probe_radius: f64) -> Result<f64, Error> {
        if !probe_radius.is_finite() {
            return Err(Error::InvalidInput(format!(
                "probe radius must be finite, got {probe_radius}"
            )));
        }
        Ok(self.count_at_least(probe_radius.max(0.0)) as f64 * self.point_volume())
    }

    /// Volume outside the van der Waals atoms, in cubic Angstrom.
    #[must_use]
    pub fn pore_volume(&self) -> f64 {
        self.count_at_least(0.0) as f64 * self.point_volume()
    }

    fn count_at_least(&self, threshold: f64) -> usize {
        self.grid
            .data
            .iter()
            .filter(|value| **value >= threshold)
            .count()
    }

    /// Approximate area of the van der Waals surface in square Angstrom.
    #[must_use]
    pub fn surface_area(&self) -> f64 {
        let [n_x, n_y, n_z] = self.grid.shape;
        let step = self.grid.step;
        let area_x = cross_length(step[1], step[2]);
        let area_y = cross_length(step[0], step[2]);
        let area_z = cross_length(step[0], step[1]);
        let mut area = 0.0;
        for z in 0..n_z {
            for y in 0..n_y {
                for x in 0..n_x {
                    let inside = self.value(x, y, z) < 0.0;
                    if inside != (self.value((x + 1) % n_x, y, z) < 0.0) {
                        area += area_x;
                    }
                    if inside != (self.value(x, (y + 1) % n_y, z) < 0.0) {
                        area += area_y;
                    }
                    if inside != (self.value(x, y, (z + 1) % n_z) < 0.0) {
                        area += area_z;
                    }
                }
            }
        }
        area
    }

    /// Diameter of the largest probe that fits anywhere on the grid.
    ///
    /// Returns `0` for a grid with no free sample and infinity for a grid
    /// without atoms.
    #[must_use]
    pub fn largest_included_diameter(&self) -> f64 {
        let mut maximum = 0.0_f64;
        for value in &self.grid.data {
            if *value > maximum {
                maximum = *value;
            }
        }
        2.0 * maximum
    }

    /// Geometric pore-size histogram from `2 * distance` at every free sample.
    ///
    /// # Errors
    ///
    /// Returns [`Error::InvalidInput`] for a non-positive or non-finite bin
    /// width, more than one million bins, or a grid without atoms (whose pore
    /// diameter is unbounded).
    pub fn pore_diameter(&self, bin_width: f64) -> Result<Vec<PoreSizeBin>, Error> {
        if !bin_width.is_finite() || bin_width <= 0.0 {
            return Err(Error::InvalidInput(format!(
                "pore-size bin width must be positive and finite, got {bin_width}"
            )));
        }
        let mut maximum = 0.0_f64;
        for value in &self.grid.data {
            if !value.is_finite() {
                return Err(Error::InvalidInput(
                    "pore diameter is unbounded for a grid without atoms".to_owned(),
                ));
            }
            if *value > maximum {
                maximum = *value;
            }
        }
        let largest = 2.0 * maximum;
        let bin_count = if largest <= 0.0 {
            1
        } else {
            (largest / bin_width).ceil() as usize
        };
        if bin_count > MAX_HISTOGRAM_BINS {
            return Err(Error::InvalidInput(format!(
                "pore-size histogram would have {bin_count} bins; use a larger bin width"
            )));
        }
        let volume = self.point_volume();
        let mut bins: Vec<PoreSizeBin> = (0..bin_count)
            .map(|index| PoreSizeBin {
                diameter_min: index as f64 * bin_width,
                diameter_max: (index + 1) as f64 * bin_width,
                count: 0,
                volume: 0.0,
            })
            .collect();
        for value in &self.grid.data {
            if *value < 0.0 {
                continue;
            }
            let diameter = 2.0 * value;
            let index = ((diameter / bin_width).floor() as usize).min(bin_count - 1);
            if let Some(bin) = bins.get_mut(index) {
                bin.count += 1;
                bin.volume += volume;
            }
        }
        Ok(bins)
    }
}

pub(super) fn cross_length(first: [f64; 3], second: [f64; 3]) -> f64 {
    let [a_x, a_y, a_z] = first;
    let [b_x, b_y, b_z] = second;
    let x = a_y * b_z - a_z * b_y;
    let y = a_z * b_x - a_x * b_z;
    let z = a_x * b_y - a_y * b_x;
    (x * x + y * y + z * z).sqrt()
}
