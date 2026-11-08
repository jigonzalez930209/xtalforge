//! UFF bond and angle terms with analytic gradients (Task 4.2.3).
//!
//! Ported from Crystal-Builder `xtal/ff/uff/terms.py`. The angle energy is
//! UFF eq 11/12, the cosine expansion whose curvature at the equilibrium
//! angle is the `K_IJK` of eq 13; near `theta0` that is the familiar
//! `0.5 * K * (theta - theta0)^2`, and the k formula is unchanged. The
//! special linear, trigonal and square-planar expansions of eq 12 are used
//! when the centre's `theta0` is exactly 180, 120 or 90 degrees.
//!
//! Units are the reference units: Angstrom, radians internally, kcal/mol and
//! kcal/mol/Angstrom. The calculator converts at its boundary.

use crate::Error;
use crate::uff::math::{
    EPS, add, add_virial, dot, negate, norm, scale, scatter, shift_vector, sub,
};
use crate::uff::params::{FORCE_CONSTANT, LAMBDA, get};

/// The angle expansion UFF uses for one equilibrium angle.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum AngleForm {
    /// `C0 + C1 cos(theta) + C2 cos(2 theta)`.
    General,
    /// `(1 + cos(theta))`, for a linear centre.
    Linear,
    /// `(1 - cos(3 theta)) / 9`, for a trigonal centre.
    Trigonal,
    /// `(1 - cos(4 theta)) / 16`, for a square-planar centre.
    Square,
}

impl AngleForm {
    /// Chooses the expansion for an equilibrium angle in degrees.
    #[must_use]
    pub fn for_theta0(theta0: f64) -> Self {
        for (form, ideal) in [
            (Self::Linear, 180.0),
            (Self::Trigonal, 120.0),
            (Self::Square, 90.0),
        ] {
            if (theta0 - ideal).abs() < 1e-6 {
                return form;
            }
        }
        Self::General
    }
}

/// UFF equilibrium bond length, eq 2-4: `r_i + r_j + r_BO - r_EN`.
///
/// The paper prints the electronegativity term with a plus sign; subtracting
/// it reproduces the paper's own worked C-H value and is the convention in
/// the reference implementation.
///
/// # Errors
///
/// Returns [`Error::UnknownType`] when either type is absent from the table.
pub fn natural_bond_length(a: &str, b: &str, order: f64) -> Result<f64, Error> {
    let first = get(a)?;
    let second = get(b)?;
    let total = first.r1 + second.r1;
    let r_bo = -LAMBDA * total * order.max(1e-6).ln();
    let numerator = first.r1 * second.r1 * (first.chi.sqrt() - second.chi.sqrt()).powi(2);
    let denominator = first.chi * first.r1 + second.chi * second.r1;
    let r_en = if denominator.abs() <= EPS {
        0.0
    } else {
        numerator / denominator
    };
    Ok(total + r_bo - r_en)
}

/// UFF bond force constant, eq 6: `664.12 Z*_i Z*_j / r0^3`, in
/// kcal/mol/Angstrom^2.
///
/// # Errors
///
/// Returns [`Error::UnknownType`] when either type is absent from the table.
pub fn bond_force_constant(a: &str, b: &str, r0: f64) -> Result<f64, Error> {
    let first = get(a)?;
    let second = get(b)?;
    if r0.abs() <= EPS {
        return Err(Error::Topology(format!(
            "zero equilibrium bond length for {a}-{b}"
        )));
    }
    Ok(FORCE_CONSTANT * first.z1 * second.z1 / r0.powi(3))
}

/// `(x_ij, D_ij)` for the Lennard-Jones 12-6 term: geometric means of the
/// two atoms' well distances and depths, eq 22.
///
/// # Errors
///
/// Returns [`Error::UnknownType`] when either type is absent from the table.
pub fn vdw_pair(a: &str, b: &str) -> Result<(f64, f64), Error> {
    let first = get(a)?;
    let second = get(b)?;
    Ok(((first.x1 * second.x1).sqrt(), (first.d1 * second.d1).sqrt()))
}

/// UFF angle force constant, eq 13, in kcal/mol/radian^2.
///
/// # Errors
///
/// Returns [`Error::UnknownType`] or a topology error for degenerate lengths.
pub fn angle_force_constant(
    a: &str,
    b: &str,
    c: &str,
    theta0: f64,
    order_ab: f64,
    order_bc: f64,
) -> Result<f64, Error> {
    let r_ab = natural_bond_length(a, b, order_ab)?;
    let r_bc = natural_bond_length(b, c, order_bc)?;
    let cos0 = theta0.to_radians().cos();
    let r_ac2 = r_ab.powi(2) + r_bc.powi(2) - 2.0 * r_ab * r_bc * cos0;
    let r_ac = r_ac2.sqrt();
    if r_ac <= EPS {
        return Err(Error::Topology(format!(
            "degenerate angle {a}-{b}-{c} at {theta0} degrees"
        )));
    }
    let za = get(a)?.z1;
    let zc = get(c)?.z1;
    Ok(FORCE_CONSTANT * za * zc / r_ac.powi(5)
        * (3.0 * r_ab * r_bc * (1.0 - cos0.powi(2)) - r_ac2 * cos0))
}

/// Coefficients of the general cosine expansion for `theta0`.
pub(crate) fn general_coefficients(theta0: f64) -> (f64, f64, f64) {
    let cos0 = theta0.to_radians().cos();
    let sin0_2 = (1.0 - cos0.powi(2)).max(1e-8);
    let k2 = 1.0 / (4.0 * sin0_2);
    let k1 = -4.0 * k2 * cos0;
    let k0 = k2 * (2.0 * cos0.powi(2) + 1.0);
    (k0, k1, k2)
}

fn vector(
    positions: &[[f64; 3]],
    i: usize,
    j: usize,
    shift: [i32; 3],
    matrix: [[f64; 3]; 3],
) -> Result<[f64; 3], Error> {
    let start = positions
        .get(i)
        .ok_or_else(|| Error::Topology(format!("atom index {i} is out of range")))?;
    let end = positions
        .get(j)
        .ok_or_else(|| Error::Topology(format!("atom index {j} is out of range")))?;
    Ok(sub(add(*end, shift_vector(shift, matrix)), *start))
}

/// Harmonic bond stretch, eq 1a: `E = 1/2 k (r - r0)^2`.
#[derive(Debug, Clone, Default)]
pub struct BondTerm {
    pub(crate) i: Vec<usize>,
    pub(crate) j: Vec<usize>,
    pub(crate) shift: Vec<[i32; 3]>,
    pub(crate) r0: Vec<f64>,
    pub(crate) k: Vec<f64>,
}

impl BondTerm {
    /// Number of bonds.
    #[must_use]
    pub fn len(&self) -> usize {
        self.i.len()
    }

    /// Whether the term has no bonds.
    #[must_use]
    pub fn is_empty(&self) -> bool {
        self.i.is_empty()
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
        let mut gradient = vec![[0.0; 3]; positions.len()];
        if self.is_empty() {
            return Ok((0.0, gradient));
        }
        let mut energy = 0.0;
        let mut virial_local = [[0.0; 3]; 3];
        for index in 0..self.len() {
            let (Some(atom_i), Some(atom_j), Some(shift), Some(r0), Some(k)) = (
                self.i.get(index),
                self.j.get(index),
                self.shift.get(index),
                self.r0.get(index),
                self.k.get(index),
            ) else {
                return Err(Error::Topology(
                    "bond term arrays have inconsistent lengths".to_owned(),
                ));
            };
            let displacement = vector(positions, *atom_i, *atom_j, *shift, matrix)?;
            let distance = norm(displacement);
            if distance <= EPS {
                continue;
            }
            let delta = distance - r0;
            energy += 0.5 * k * delta * delta;
            let gradient_vector = scale(displacement, k * delta / distance);
            add_virial(&mut virial_local, displacement, gradient_vector);
            scatter(&mut gradient, *atom_j, gradient_vector);
            scatter(&mut gradient, *atom_i, negate(gradient_vector));
        }
        if let Some(target) = virial {
            for row in 0..3 {
                for column in 0..3 {
                    if let (Some(slot), Some(value)) = (
                        target
                            .get_mut(row)
                            .and_then(|values| values.get_mut(column)),
                        virial_local.get(row).and_then(|values| values.get(column)),
                    ) {
                        *slot += value;
                    }
                }
            }
        }
        Ok((energy, gradient))
    }
}

/// Angle bend about atom `j`, eq 11/12.
#[derive(Debug, Clone, Default)]
pub struct AngleTerm {
    pub(crate) i: Vec<usize>,
    pub(crate) j: Vec<usize>,
    pub(crate) k: Vec<usize>,
    pub(crate) shift_i: Vec<[i32; 3]>,
    pub(crate) shift_k: Vec<[i32; 3]>,
    pub(crate) force: Vec<f64>,
    pub(crate) c0: Vec<f64>,
    pub(crate) c1: Vec<f64>,
    pub(crate) c2: Vec<f64>,
    pub(crate) form: Vec<AngleForm>,
}

impl AngleTerm {
    /// Number of angles.
    #[must_use]
    pub fn len(&self) -> usize {
        self.j.len()
    }

    /// Whether the term has no angles.
    #[must_use]
    pub fn is_empty(&self) -> bool {
        self.j.is_empty()
    }

    /// Energy in kcal/mol at these positions.
    ///
    /// # Errors
    ///
    /// Returns [`Error::Topology`] for a position list that does not match.
    pub fn energy(&self, positions: &[[f64; 3]], matrix: [[f64; 3]; 3]) -> Result<f64, Error> {
        Ok(self.energy_and_gradient(positions, matrix, None)?.0)
    }

    fn energy_of_cos(form: AngleForm, cos: f64, c0: f64, c1: f64, c2: f64) -> (f64, f64) {
        match form {
            AngleForm::General => (
                c0 + c1 * cos + c2 * (2.0 * cos * cos - 1.0),
                c1 + 4.0 * c2 * cos,
            ),
            AngleForm::Linear => (1.0 + cos, 1.0),
            AngleForm::Trigonal => {
                let cos3 = 4.0 * cos.powi(3) - 3.0 * cos;
                ((1.0 - cos3) / 9.0, -(12.0 * cos * cos - 3.0) / 9.0)
            }
            AngleForm::Square => {
                let cos4 = 8.0 * cos.powi(4) - 8.0 * cos * cos + 1.0;
                (
                    (1.0 - cos4) / 16.0,
                    -(32.0 * cos.powi(3) - 16.0 * cos) / 16.0,
                )
            }
        }
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
        let mut gradient = vec![[0.0; 3]; positions.len()];
        if self.is_empty() {
            return Ok((0.0, gradient));
        }
        let mut energy = 0.0;
        let mut virial_local = [[0.0; 3]; 3];
        for index in 0..self.len() {
            let (Some(atom_i), Some(atom_j), Some(atom_k), Some(shift_i), Some(shift_k)) = (
                self.i.get(index),
                self.j.get(index),
                self.k.get(index),
                self.shift_i.get(index),
                self.shift_k.get(index),
            ) else {
                return Err(Error::Topology(
                    "angle term arrays have inconsistent lengths".to_owned(),
                ));
            };
            let (Some(force), Some(c0), Some(c1), Some(c2), Some(form)) = (
                self.force.get(index),
                self.c0.get(index),
                self.c1.get(index),
                self.c2.get(index),
                self.form.get(index),
            ) else {
                return Err(Error::Topology(
                    "angle term parameters have inconsistent lengths".to_owned(),
                ));
            };
            let u = vector(positions, *atom_j, *atom_i, *shift_i, matrix)?;
            let v = vector(positions, *atom_j, *atom_k, *shift_k, matrix)?;
            let nu = norm(u);
            let nv = norm(v);
            if nu <= EPS || nv <= EPS {
                continue;
            }
            let cos = (dot(u, v) / (nu * nv)).clamp(-1.0, 1.0);
            let (energy_of_cos, dcos) = Self::energy_of_cos(*form, cos, *c0, *c1, *c2);
            energy += force * energy_of_cos;
            let dcos = force * dcos;
            let du = sub(scale(v, 1.0 / (nu * nv)), scale(u, cos / (nu * nu)));
            let dv = sub(scale(u, 1.0 / (nu * nv)), scale(v, cos / (nv * nv)));
            let gradient_i = scale(du, dcos);
            let gradient_k = scale(dv, dcos);
            add_virial(&mut virial_local, u, gradient_i);
            add_virial(&mut virial_local, v, gradient_k);
            scatter(&mut gradient, *atom_i, gradient_i);
            scatter(&mut gradient, *atom_k, gradient_k);
            scatter(&mut gradient, *atom_j, negate(add(gradient_i, gradient_k)));
        }
        if let Some(target) = virial {
            for row in 0..3 {
                for column in 0..3 {
                    if let (Some(slot), Some(value)) = (
                        target
                            .get_mut(row)
                            .and_then(|values| values.get_mut(column)),
                        virial_local.get(row).and_then(|values| values.get(column)),
                    ) {
                        *slot += value;
                    }
                }
            }
        }
        Ok((energy, gradient))
    }
}
