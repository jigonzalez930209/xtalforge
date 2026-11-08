//! A scalar field on a regular grid (Task 3.5.1).

#![allow(clippy::indexing_slicing)] // fixed-format parser, lengths checked above
use serde::{Deserialize, Serialize};
use xtalforge_core::Error;

/// A scalar field sampled at `shape[0] * shape[1] * shape[2]` points.
///
/// `origin` is the Cartesian position of sample `(0, 0, 0)` and `step` holds
/// the three Cartesian step vectors (rows). Data is stored x-fastest.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct Grid3 {
    /// Cartesian origin in Angstrom.
    pub origin: [f64; 3],
    /// Step vectors in Angstrom.
    pub step: [[f64; 3]; 3],
    /// Number of samples along x, y and z.
    pub shape: [usize; 3],
    /// Sample values.
    pub data: Vec<f64>,
    /// Whether the grid wraps at its boundary.
    pub periodic: bool,
}

impl Grid3 {
    /// Creates a grid and checks the data length.
    ///
    /// # Errors
    ///
    /// Returns [`Error::InvalidInput`] for a zero dimension, non-finite values
    /// or a data length that does not match the shape.
    pub fn new(
        origin: [f64; 3],
        step: [[f64; 3]; 3],
        shape: [usize; 3],
        data: Vec<f64>,
        periodic: bool,
    ) -> Result<Self, Error> {
        if shape.contains(&0) {
            return Err(Error::InvalidInput(
                "a grid needs three positive dimensions".to_owned(),
            ));
        }
        if !origin
            .iter()
            .chain(step.iter().flatten())
            .all(|value| value.is_finite())
        {
            return Err(Error::InvalidInput(
                "grid origin and steps must be finite".to_owned(),
            ));
        }
        let expected = shape[0]
            .checked_mul(shape[1])
            .and_then(|value| value.checked_mul(shape[2]))
            .ok_or_else(|| Error::InvalidInput("grid dimensions overflow".to_owned()))?;
        if data.len() != expected {
            return Err(Error::InvalidInput(format!(
                "grid has {} samples for shape {shape:?}",
                data.len()
            )));
        }
        Ok(Self {
            origin,
            step,
            shape,
            data,
            periodic,
        })
    }

    /// Total number of samples.
    #[must_use]
    pub fn len(&self) -> usize {
        self.data.len()
    }

    /// Whether the grid has no samples (never true for a valid grid).
    #[must_use]
    pub fn is_empty(&self) -> bool {
        self.data.is_empty()
    }

    /// Linear index of a sample, or `None` outside the shape.
    #[must_use]
    pub fn index(&self, x: usize, y: usize, z: usize) -> Option<usize> {
        if x >= self.shape[0] || y >= self.shape[1] || z >= self.shape[2] {
            return None;
        }
        Some((z * self.shape[1] + y) * self.shape[0] + x)
    }

    /// Sample value, or `None` outside the shape.
    #[must_use]
    pub fn get(&self, x: usize, y: usize, z: usize) -> Option<f64> {
        self.index(x, y, z)
            .and_then(|index| self.data.get(index).copied())
    }

    /// Value at a fractional grid coordinate, with periodic wrapping when
    /// `periodic` is set.
    #[must_use]
    pub fn sample(&self, point: [f64; 3]) -> Option<f64> {
        let mut indices = [0usize; 3];
        for (axis, value) in point.iter().enumerate() {
            let scaled = value * self.shape[axis] as f64;
            let wrapped = if self.periodic {
                scaled.rem_euclid(self.shape[axis] as f64)
            } else if (0.0..self.shape[axis] as f64).contains(&scaled) {
                scaled
            } else {
                return None;
            };
            *indices.get_mut(axis)? = wrapped.floor() as usize;
        }
        self.get(indices[0], indices[1], indices[2])
    }

    /// Cartesian position of a sample.
    #[must_use]
    pub fn position(&self, x: usize, y: usize, z: usize) -> [f64; 3] {
        [
            self.origin[0]
                + x as f64 * self.step[0][0]
                + y as f64 * self.step[1][0]
                + z as f64 * self.step[2][0],
            self.origin[1]
                + x as f64 * self.step[0][1]
                + y as f64 * self.step[1][1]
                + z as f64 * self.step[2][1],
            self.origin[2]
                + x as f64 * self.step[0][2]
                + y as f64 * self.step[1][2]
                + z as f64 * self.step[2][2],
        ]
    }

    /// Minimum and maximum sample values.
    #[must_use]
    pub fn range(&self) -> (f64, f64) {
        let mut minimum = f64::INFINITY;
        let mut maximum = f64::NEG_INFINITY;
        for value in &self.data {
            minimum = minimum.min(*value);
            maximum = maximum.max(*value);
        }
        (minimum, maximum)
    }
}
