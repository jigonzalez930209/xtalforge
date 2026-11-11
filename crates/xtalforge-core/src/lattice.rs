//! Unit cell geometry: the lattice matrix and everything derived from it.
//!
//! Conventions (Task 1.2.1): the matrix **rows** are the lattice vectors
//! a, b, c in Angstrom; fractional to Cartesian is `frac @ M`; the standard
//! orientation puts a along +x and b in the xy plane; the reciprocal lattice
//! uses the crystallographic convention without 2 pi, so `a . a* = 1` and
//! `d(hkl) = 1 / |hkl @ reciprocal|`.

use serde::{Deserialize, Serialize};

use crate::Error;

mod functions;
mod math;
pub use functions::{lattice_info, lattice_minimum_image};
use math::{
    add, angle_between, determinant, inverse, multiply, multiply_row, norm, subtract, transpose,
};

const DEG: f64 = std::f64::consts::PI / 180.0;
const MIN_VOLUME: f64 = 1e-9;
const DEFAULT_WRAP_TOLERANCE: f64 = 1e-8;

/// A unit cell. The matrix rows are the vectors a, b and c.
#[derive(Debug, Clone, Copy, PartialEq, Serialize, Deserialize)]
#[serde(try_from = "[[f64; 3]; 3]", into = "[[f64; 3]; 3]")]
pub struct Lattice {
    matrix: [[f64; 3]; 3],
}

impl Lattice {
    /// Creates a lattice from a matrix whose rows are a, b and c.
    ///
    /// # Errors
    ///
    /// Returns [`Error::InvalidInput`] for non-finite entries or a singular
    /// matrix (zero volume). A negative determinant is allowed and reported by
    /// [`Lattice::is_right_handed`].
    pub fn new(matrix: [[f64; 3]; 3]) -> Result<Self, Error> {
        if !matrix.iter().flatten().all(|value| value.is_finite()) {
            return Err(Error::InvalidInput(
                "lattice matrix contains non-finite values".to_owned(),
            ));
        }
        let determinant = determinant(matrix);
        if determinant.abs() < MIN_VOLUME {
            return Err(Error::InvalidInput(format!(
                "lattice matrix is singular (zero volume): determinant={determinant}"
            )));
        }
        Ok(Self { matrix })
    }

    /// Builds the standard orientation from `a, b, c` (Angstrom) and
    /// `alpha, beta, gamma` (degrees).
    ///
    /// # Errors
    ///
    /// Returns [`Error::InvalidInput`] for non-positive lengths, angles
    /// outside `(0, 180)`, non-finite values or geometrically impossible
    /// angles.
    pub fn from_parameters(
        a: f64,
        b: f64,
        c: f64,
        alpha: f64,
        beta: f64,
        gamma: f64,
    ) -> Result<Self, Error> {
        let values = [a, b, c, alpha, beta, gamma];
        if !values.iter().all(|value| value.is_finite()) {
            return Err(Error::InvalidInput(
                "cell parameters must be finite numbers".to_owned(),
            ));
        }
        for (name, value) in [("a", a), ("b", b), ("c", c)] {
            if value <= 0.0 {
                return Err(Error::InvalidInput(format!(
                    "cell edge {name} must be positive, got {value}"
                )));
            }
        }
        for (name, value) in [("alpha", alpha), ("beta", beta), ("gamma", gamma)] {
            if !(0.0 < value && value < 180.0) {
                return Err(Error::InvalidInput(format!(
                    "cell angle {name} must be in (0, 180) degrees, got {value}"
                )));
            }
        }

        let cos_alpha = (alpha * DEG).cos();
        let cos_beta = (beta * DEG).cos();
        let cos_gamma = (gamma * DEG).cos();
        let sin_gamma = (gamma * DEG).sin();
        let cx = cos_beta;
        let cy = (cos_alpha - cos_beta * cos_gamma) / sin_gamma;
        let cz_squared = 1.0 - cx * cx - cy * cy;
        if cz_squared <= 0.0 || !cz_squared.is_finite() {
            return Err(Error::InvalidInput(format!(
                "cell angles are geometrically impossible: \
                 alpha={alpha}, beta={beta}, gamma={gamma}"
            )));
        }

        Self::new([
            [a, 0.0, 0.0],
            [b * cos_gamma, b * sin_gamma, 0.0],
            [c * cx, c * cy, c * cz_squared.sqrt()],
        ])
    }

    /// Cubic cell with edge `a`.
    ///
    /// # Errors
    ///
    /// Returns [`Error::InvalidInput`] when `a` is not positive and finite.
    pub fn cubic(a: f64) -> Result<Self, Error> {
        Self::from_parameters(a, a, a, 90.0, 90.0, 90.0)
    }

    /// Orthorhombic cell with edges `a, b, c`.
    ///
    /// # Errors
    ///
    /// Returns [`Error::InvalidInput`] when an edge is not positive and finite.
    pub fn orthorhombic(a: f64, b: f64, c: f64) -> Result<Self, Error> {
        Self::from_parameters(a, b, c, 90.0, 90.0, 90.0)
    }

    /// The matrix whose rows are a, b and c.
    #[must_use]
    pub fn matrix(&self) -> [[f64; 3]; 3] {
        self.matrix
    }

    /// Edge lengths `(a, b, c)` in Angstrom.
    #[must_use]
    pub fn lengths(&self) -> [f64; 3] {
        self.matrix.map(norm)
    }

    /// Angles `(alpha, beta, gamma)` in degrees.
    #[must_use]
    pub fn angles(&self) -> [f64; 3] {
        let [a, b, c] = self.lengths();
        let [row_a, row_b, row_c] = self.matrix;
        [
            angle_between(row_b, row_c, b * c),
            angle_between(row_a, row_c, a * c),
            angle_between(row_a, row_b, a * b),
        ]
    }

    /// `(a, b, c, alpha, beta, gamma)` in Angstrom and degrees.
    #[must_use]
    pub fn parameters(&self) -> [f64; 6] {
        let [a, b, c] = self.lengths();
        let [alpha, beta, gamma] = self.angles();
        [a, b, c, alpha, beta, gamma]
    }

    /// Cell volume in cubic Angstrom.
    #[must_use]
    pub fn volume(&self) -> f64 {
        determinant(self.matrix).abs()
    }

    /// Inverse matrix: Cartesian to fractional is `cart @ inverse`.
    #[must_use]
    pub fn inverse(&self) -> [[f64; 3]; 3] {
        inverse(self.matrix)
    }

    /// Metric tensor `M M^T`; `frac @ G @ frac` is a squared length.
    #[must_use]
    pub fn metric_tensor(&self) -> [[f64; 3]; 3] {
        multiply(self.matrix, transpose(self.matrix))
    }

    /// Reciprocal cell whose rows are a*, b*, c*.
    ///
    /// The crystallographic convention is the default (`a . a* = 1`);
    /// `two_pi = true` multiplies by 2 pi.
    ///
    /// # Errors
    ///
    /// Returns [`Error::InvalidInput`] if the reciprocal matrix is singular,
    /// which cannot happen for a valid lattice.
    pub fn reciprocal(&self, two_pi: bool) -> Result<Self, Error> {
        let mut reciprocal = transpose(self.inverse());
        if two_pi {
            let factor = 2.0 * std::f64::consts::PI;
            for row in &mut reciprocal {
                for value in row {
                    *value *= factor;
                }
            }
        }
        Self::new(reciprocal)
    }

    /// Interplanar spacing for a Miller index, in Angstrom.
    ///
    /// # Errors
    ///
    /// Returns [`Error::InvalidInput`] for `(0, 0, 0)` or a non-finite index.
    pub fn d_spacing(&self, hkl: [f64; 3]) -> Result<f64, Error> {
        let vector = self.plane_vector(hkl)?;
        Ok(1.0 / norm(vector))
    }

    /// Unit Cartesian normal to the planes `(hkl)`.
    ///
    /// # Errors
    ///
    /// Returns [`Error::InvalidInput`] for `(0, 0, 0)` or a non-finite index.
    pub fn plane_normal(&self, hkl: [f64; 3]) -> Result<[f64; 3], Error> {
        let vector = self.plane_vector(hkl)?;
        let length = norm(vector);
        let [x, y, z] = vector;
        Ok([x / length, y / length, z / length])
    }

    /// Fractional to Cartesian for one point.
    #[must_use]
    pub fn to_cart(&self, frac: [f64; 3]) -> [f64; 3] {
        multiply_row(frac, self.matrix)
    }

    /// Fractional to Cartesian for a block of points.
    #[must_use]
    pub fn to_cart_many(&self, frac: &[[f64; 3]]) -> Vec<[f64; 3]> {
        frac.iter().map(|point| self.to_cart(*point)).collect()
    }

    /// Cartesian to fractional for one point.
    #[must_use]
    pub fn to_frac(&self, cart: [f64; 3]) -> [f64; 3] {
        multiply_row(cart, self.inverse())
    }

    /// Cartesian to fractional for a block of points.
    #[must_use]
    pub fn to_frac_many(&self, cart: &[[f64; 3]]) -> Vec<[f64; 3]> {
        cart.iter().map(|point| self.to_frac(*point)).collect()
    }

    /// Wraps fractional coordinates into `[0, 1)`.
    #[must_use]
    pub fn wrap(frac: [f64; 3]) -> [f64; 3] {
        // `rem_euclid` rounds a tiny negative value such as -1e-22 up to
        // exactly 1.0, which is outside [0, 1) and duplicates the atom of the
        // neighbouring cell; that value maps to 0.
        frac.map(|value| {
            let wrapped = value.rem_euclid(1.0);
            if wrapped >= 1.0 { 0.0 } else { wrapped }
        })
    }

    /// Wraps fractional coordinates and maps values within `tolerance` of 1
    /// to 0.
    #[must_use]
    pub fn wrap_with_tolerance(frac: [f64; 3], tolerance: f64) -> [f64; 3] {
        let tolerance = tolerance.max(0.0);
        frac.map(|value| {
            let wrapped = value.rem_euclid(1.0);
            if wrapped >= 1.0 - tolerance || wrapped >= 1.0 {
                0.0
            } else {
                wrapped
            }
        })
    }

    /// Minimum-image fractional vector and its Cartesian length.
    ///
    /// The search is exact for any non-singular cell: it starts from the
    /// rounded fractional shift and enumerates every integer shift that could
    /// beat the current distance, using the bound
    /// `|x| <= |x M| * ||M^-1||_F`. A strongly skewed triclinic cell therefore
    /// searches more than the 27 neighbouring images when it has to.
    #[must_use]
    pub fn minimum_image(&self, delta: [f64; 3]) -> ([f64; 3], f64) {
        let inverted = self.inverse();
        let inverse_norm = inverted
            .iter()
            .flatten()
            .map(|value| value * value)
            .sum::<f64>()
            .sqrt();
        let rounded = delta.map(f64::round);
        let mut best = subtract(delta, rounded);
        let mut best_distance = norm(self.to_cart(best));
        if best_distance == 0.0 {
            return (best, 0.0);
        }
        let radius = best_distance * inverse_norm;
        let [(i_min, i_max), (j_min, j_max), (k_min, k_max)] = delta.map(|value| {
            (
                (value - radius).floor() as i64,
                (value + radius).ceil() as i64,
            )
        });
        for i in i_min..=i_max {
            for j in j_min..=j_max {
                for k in k_min..=k_max {
                    let [dx, dy, dz] = delta;
                    let candidate = [dx - i as f64, dy - j as f64, dz - k as f64];
                    let distance = norm(self.to_cart(candidate));
                    if distance < best_distance {
                        best_distance = distance;
                        best = candidate;
                    }
                }
            }
        }
        (best, best_distance)
    }

    /// Minimum-image distance between two fractional points.
    #[must_use]
    pub fn distance(&self, first: [f64; 3], second: [f64; 3]) -> f64 {
        self.minimum_image(subtract(first, second)).1
    }

    /// Applies a cell transformation: `M' = P M`.
    ///
    /// Rows of `P` are the new vectors in units of the old ones.
    ///
    /// # Errors
    ///
    /// Returns [`Error::InvalidInput`] for a non-finite matrix, a singular
    /// matrix or a negative determinant (which would flip handedness).
    pub fn transform(&self, p: [[f64; 3]; 3]) -> Result<Self, Error> {
        if !p.iter().flatten().all(|value| value.is_finite()) {
            return Err(Error::InvalidInput(
                "transformation matrix contains non-finite values".to_owned(),
            ));
        }
        let determinant = determinant(p);
        if determinant.abs() < MIN_VOLUME {
            return Err(Error::InvalidInput(
                "transformation matrix is singular".to_owned(),
            ));
        }
        if determinant < 0.0 {
            return Err(Error::InvalidInput(
                "transformation matrix has negative determinant (it would flip handedness)"
                    .to_owned(),
            ));
        }
        Self::new(multiply(p, self.matrix))
    }

    /// Applies a Cartesian strain tensor: `M' = M (I + eps)^T`.
    ///
    /// # Errors
    ///
    /// Returns [`Error::InvalidInput`] for non-finite values or a resulting
    /// singular cell.
    pub fn strained(&self, eps: [[f64; 3]; 3]) -> Result<Self, Error> {
        if !eps.iter().flatten().all(|value| value.is_finite()) {
            return Err(Error::InvalidInput(
                "strain tensor contains non-finite values".to_owned(),
            ));
        }
        let identity = [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]];
        let deformation = transpose(add(identity, eps));
        Self::new(multiply(self.matrix, deformation))
    }

    /// Isotropic strain convenience wrapper.
    ///
    /// # Errors
    ///
    /// Returns [`Error::InvalidInput`] for non-finite values or a singular
    /// resulting cell.
    pub fn strained_isotropic(&self, eps: f64) -> Result<Self, Error> {
        self.strained([[eps, 0.0, 0.0], [0.0, eps, 0.0], [0.0, 0.0, eps]])
    }

    /// Isotropically rescales the cell to a target volume.
    ///
    /// # Errors
    ///
    /// Returns [`Error::InvalidInput`] when the target is not positive and
    /// finite.
    pub fn scaled_to_volume(&self, volume: f64) -> Result<Self, Error> {
        if !volume.is_finite() || volume <= 0.0 {
            return Err(Error::InvalidInput(format!(
                "target volume must be positive and finite, got {volume}"
            )));
        }
        let factor = (volume / self.volume()).cbrt();
        Self::new(self.matrix.map(|row| row.map(|value| value * factor)))
    }

    /// Copy with selected parameters replaced, in standard orientation.
    ///
    /// # Errors
    ///
    /// Returns the errors of [`Lattice::from_parameters`].
    pub fn with_parameters(
        &self,
        a: Option<f64>,
        b: Option<f64>,
        c: Option<f64>,
        alpha: Option<f64>,
        beta: Option<f64>,
        gamma: Option<f64>,
    ) -> Result<Self, Error> {
        let [old_a, old_b, old_c, old_alpha, old_beta, old_gamma] = self.parameters();
        Self::from_parameters(
            a.unwrap_or(old_a),
            b.unwrap_or(old_b),
            c.unwrap_or(old_c),
            alpha.unwrap_or(old_alpha),
            beta.unwrap_or(old_beta),
            gamma.unwrap_or(old_gamma),
        )
    }

    /// Same parameters rebuilt with a along +x.
    ///
    /// # Errors
    ///
    /// Returns the errors of [`Lattice::from_parameters`].
    pub fn standard_orientation(&self) -> Result<Self, Error> {
        let [a, b, c, alpha, beta, gamma] = self.parameters();
        Self::from_parameters(a, b, c, alpha, beta, gamma)
    }

    /// Whether the matrix determinant is positive.
    #[must_use]
    pub fn is_right_handed(&self) -> bool {
        determinant(self.matrix) > 0.0
    }

    /// Element-wise comparison with an absolute tolerance.
    #[must_use]
    pub fn almost_equal(&self, other: &Self, tolerance: f64) -> bool {
        self.matrix
            .iter()
            .flatten()
            .zip(other.matrix.iter().flatten())
            .all(|(first, second)| (first - second).abs() <= tolerance)
    }

    /// Default-wrap tolerance used by [`Lattice::wrap`] callers.
    #[must_use]
    pub fn default_wrap_tolerance() -> f64 {
        DEFAULT_WRAP_TOLERANCE
    }

    fn plane_vector(&self, hkl: [f64; 3]) -> Result<[f64; 3], Error> {
        if !hkl.iter().all(|value| value.is_finite()) {
            return Err(Error::InvalidInput(
                "Miller index must contain finite numbers".to_owned(),
            ));
        }
        let reciprocal = self.reciprocal(false)?;
        let vector = multiply_row(hkl, reciprocal.matrix);
        if norm(vector) == 0.0 {
            return Err(Error::InvalidInput("(0 0 0) is not a plane".to_owned()));
        }
        Ok(vector)
    }
}

impl TryFrom<[[f64; 3]; 3]> for Lattice {
    type Error = Error;

    fn try_from(matrix: [[f64; 3]; 3]) -> Result<Self, Self::Error> {
        Self::new(matrix)
    }
}

impl From<Lattice> for [[f64; 3]; 3] {
    fn from(lattice: Lattice) -> Self {
        lattice.matrix
    }
}
