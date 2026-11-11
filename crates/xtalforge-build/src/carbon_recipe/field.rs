//! Convolution surface around a net (Task 6.8.1).
//!
//! The field is the sum of Gaussians over beads spaced along every net
//! edge, not a distance to the nearest edge: summed Gaussians are smooth
//! at a node, where several struts meet, while a distance field creases
//! there. The level is the value a distance `radius` from an endless row
//! of beads has, so along an edge the tube radius is the strut radius.

#![allow(clippy::indexing_slicing)] // bounded local arrays

use std::collections::HashMap;

use xtalforge_analysis::isosurface;
use xtalforge_core::{Error, Lattice};
use xtalforge_io::grid::Grid3;

use super::mesh::{Mesh, weld};
use super::net::{Skeleton, norm, transform};

/// The Gaussian width as a fraction of the strut radius.
pub(crate) const WIDTH_FRACTION: f64 = 0.5;
/// Points further than this many widths get nothing from a bead.
const CUTOFF_WIDTHS: f64 = 5.0;
/// Marching spacing as a fraction of the width, never finer than the floor.
const GRID_WIDTHS: f64 = 0.6;
const GRID_FLOOR: f64 = 1.0;

/// The summed-Gaussian field of a skeleton, sampled on a periodic grid.
#[derive(Debug, Clone)]
pub(crate) struct GaussianField {
    pub(crate) lattice: Lattice,
    pub(crate) width: f64,
    bead_frac: Vec<[f64; 3]>,
    bins: [usize; 3],
    bead_bins: HashMap<[i32; 3], Vec<usize>>,
}

impl GaussianField {
    /// Builds the field and the grid resolution it will be sampled at.
    pub(crate) fn new(skeleton: &Skeleton, width: f64) -> Result<Self, Error> {
        if !width.is_finite() || width <= 0.0 {
            return Err(Error::InvalidInput(
                "the carbon Gaussian width must be positive and finite".to_owned(),
            ));
        }
        let bead_frac: Vec<[f64; 3]> = skeleton
            .beads
            .iter()
            .map(|bead| {
                skeleton
                    .lattice
                    .to_frac(*bead)
                    .map(|value| value.rem_euclid(1.0))
            })
            .collect();
        let lengths = lattice_lengths(&skeleton.lattice);
        let cutoff = (CUTOFF_WIDTHS * width).max(0.5);
        let bins = [
            ((lengths[0] / cutoff).ceil() as usize).max(2),
            ((lengths[1] / cutoff).ceil() as usize).max(2),
            ((lengths[2] / cutoff).ceil() as usize).max(2),
        ];
        let mut bead_bins: HashMap<[i32; 3], Vec<usize>> = HashMap::new();
        for (index, bead) in bead_frac.iter().enumerate() {
            let bin = bin_of(*bead, bins);
            bead_bins.entry(bin).or_default().push(index);
        }
        Ok(Self {
            lattice: skeleton.lattice,
            width,
            bead_frac,
            bins,
            bead_bins,
        })
    }

    /// Beads and lattice shifts within the cutoff of a point.
    fn near(&self, frac: [f64; 3]) -> Vec<(usize, [i32; 3])> {
        let base = bin_of(frac, self.bins);
        let mut out = Vec::new();
        for shift_x in -1..=1i32 {
            for shift_y in -1..=1i32 {
                for shift_z in -1..=1i32 {
                    let raw = [base[0] + shift_x, base[1] + shift_y, base[2] + shift_z];
                    let mut bin = [0i32; 3];
                    let mut shift = [0i32; 3];
                    for axis in 0..3 {
                        let size = self.bins[axis] as i32;
                        bin[axis] = raw[axis].rem_euclid(size);
                        shift[axis] = raw[axis].div_euclid(size);
                    }
                    if let Some(beads) = self.bead_bins.get(&bin) {
                        for bead in beads {
                            out.push((*bead, shift));
                        }
                    }
                }
            }
        }
        out
    }

    /// Field value at one fractional point.
    pub(crate) fn value(&self, frac: [f64; 3]) -> f64 {
        let matrix = self.lattice.matrix();
        let cutoff2 = (CUTOFF_WIDTHS * self.width).powi(2);
        let two_width2 = 2.0 * self.width * self.width;
        let mut value = 0.0;
        for (index, shift) in self.near(frac) {
            let Some(bead) = self.bead_frac.get(index) else {
                continue;
            };
            let delta = [
                bead[0] + f64::from(shift[0]) - frac[0],
                bead[1] + f64::from(shift[1]) - frac[1],
                bead[2] + f64::from(shift[2]) - frac[2],
            ];
            let cart = transform(delta, matrix);
            let distance2 = cart[0] * cart[0] + cart[1] * cart[1] + cart[2] * cart[2];
            if distance2 <= cutoff2 {
                value += (-distance2 / two_width2).exp();
            }
        }
        value
    }

    /// Field value and Cartesian gradient at one fractional point.
    pub(crate) fn gaussian_value_gradient(&self, frac: [f64; 3]) -> (f64, [f64; 3]) {
        let matrix = self.lattice.matrix();
        let point = transform(frac, matrix);
        let cutoff2 = (CUTOFF_WIDTHS * self.width).powi(2);
        let two_width2 = 2.0 * self.width * self.width;
        let mut value = 0.0;
        let mut gradient = [0.0; 3];
        for (index, shift) in self.near(frac) {
            let Some(bead) = self.bead_frac.get(index) else {
                continue;
            };
            let delta = [
                bead[0] + f64::from(shift[0]) - frac[0],
                bead[1] + f64::from(shift[1]) - frac[1],
                bead[2] + f64::from(shift[2]) - frac[2],
            ];
            let offset = transform(delta, matrix);
            let distance2 = offset[0] * offset[0] + offset[1] * offset[1] + offset[2] * offset[2];
            if distance2 <= cutoff2 {
                let weight = (-distance2 / two_width2).exp();
                value += weight;
                let bead_cart = [
                    point[0] + offset[0],
                    point[1] + offset[1],
                    point[2] + offset[2],
                ];
                for axis in 0..3 {
                    gradient[axis] +=
                        weight * (bead_cart[axis] - point[axis]) / (self.width * self.width);
                }
            }
        }
        (value, gradient)
    }

    /// Marches the level set at an explicit spacing and welds it.
    pub(crate) fn march(&self, level: f64, spacing: f64) -> Result<Mesh, Error> {
        if !spacing.is_finite() || spacing <= 0.0 {
            return Err(Error::InvalidInput(
                "the marching spacing must be positive and finite".to_owned(),
            ));
        }
        let lengths = lattice_lengths(&self.lattice);
        let bins = [
            ((lengths[0] / spacing).ceil() as usize).max(4),
            ((lengths[1] / spacing).ceil() as usize).max(4),
            ((lengths[2] / spacing).ceil() as usize).max(4),
        ];
        let shape = [bins[0] + 1, bins[1] + 1, bins[2] + 1];
        let mut values = Vec::with_capacity(shape[0] * shape[1] * shape[2]);
        for k in 0..shape[2] {
            for j in 0..shape[1] {
                for i in 0..shape[0] {
                    let point = [
                        i as f64 / bins[0] as f64,
                        j as f64 / bins[1] as f64,
                        k as f64 / bins[2] as f64,
                    ];
                    values.push(self.value(point));
                }
            }
        }
        let matrix = self.lattice.matrix();
        let step = [
            [
                matrix[0][0] / bins[0] as f64,
                matrix[0][1] / bins[0] as f64,
                matrix[0][2] / bins[0] as f64,
            ],
            [
                matrix[1][0] / bins[1] as f64,
                matrix[1][1] / bins[1] as f64,
                matrix[1][2] / bins[1] as f64,
            ],
            [
                matrix[2][0] / bins[2] as f64,
                matrix[2][1] / bins[2] as f64,
                matrix[2][2] / bins[2] as f64,
            ],
        ];
        let grid = Grid3::new([0.0; 3], step, shape, values, true)?;
        let mesh = isosurface(&grid, level)?;
        weld(&mesh.vertices, &mesh.triangles, &self.lattice, shape)
    }

    /// Marches the level set on the field's natural grid.
    pub(crate) fn march_natural(&self, level: f64) -> Result<Mesh, Error> {
        let spacing = (GRID_WIDTHS * self.width).max(GRID_FLOOR);
        self.march(level, spacing)
    }
}

/// Field value of an endless row of beads at `radius`.
pub(crate) fn level_for(radius: f64, width: f64, spacing: f64) -> f64 {
    let reach = ((CUTOFF_WIDTHS * width / spacing).ceil() as i64) + 1;
    let mut value = 0.0;
    for step in -reach..=reach {
        let along = step as f64 * spacing;
        value += (-(radius * radius + along * along) / (2.0 * width * width)).exp();
    }
    value
}

fn bin_of(frac: [f64; 3], bins: [usize; 3]) -> [i32; 3] {
    [
        (frac[0] * bins[0] as f64).floor() as i32,
        (frac[1] * bins[1] as f64).floor() as i32,
        (frac[2] * bins[2] as f64).floor() as i32,
    ]
}

fn lattice_lengths(lattice: &Lattice) -> [f64; 3] {
    let matrix = lattice.matrix();
    [norm(matrix[0]), norm(matrix[1]), norm(matrix[2])]
}

/// A periodic scalar field whose level set is remeshed and relaxed.
pub(crate) trait LevelSet {
    /// Field value and Cartesian gradient at one fractional point.
    fn value_gradient(&self, frac: [f64; 3]) -> (f64, [f64; 3]);
}

impl LevelSet for GaussianField {
    fn value_gradient(&self, frac: [f64; 3]) -> (f64, [f64; 3]) {
        self.gaussian_value_gradient(frac)
    }
}
