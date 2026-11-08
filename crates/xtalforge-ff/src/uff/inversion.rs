//! UFF inversion (out-of-plane) term, eq 19 (Task 4.2.3).
//!
//! Ported from Crystal-Builder `xtal/ff/uff/terms.py`. `omega` is the angle
//! between the bond to the out-of-plane atom and the plane of the other
//! three; the energy is written through `sin(omega)`, which is what the
//! geometry gives directly.

use crate::Error;
use crate::uff::math::{
    EPS, MIN_SINE, add, add_virial, cross, dot, merge_virial, negate, norm, scale, scatter,
    shift_vector, sine_from_cross, sub,
};

/// sp2 centres: total barrier shared by the three inversions at the atom.
const SP2_INVERSION: f64 = 6.0;
/// UFF's harder barrier for an sp2 carbon carrying a carbonyl oxygen.
const SP2_INVERSION_WITH_O2: f64 = 50.0;

fn is_sp2_centre(name: &str) -> bool {
    matches!(name, "C_2" | "C_R" | "N_2" | "N_R" | "O_2" | "O_R")
}

fn group_15(centre: &str) -> Option<(f64, f64)> {
    match centre {
        "P_3+3" => Some((84.4339, 22.0)),
        "As3+3" => Some((86.9735, 22.0)),
        "Sb3+3" => Some((87.7047, 22.0)),
        "Bi3+3" => Some((90.0, 22.0)),
        _ => None,
    }
}

/// Coefficients `(K, C0, C1, C2)` for one three-coordinate centre.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct InversionParameters {
    /// Force constant, already divided by three.
    pub force: f64,
    /// Constant term of the cosine expansion.
    pub c0: f64,
    /// `cos(omega)` coefficient.
    pub c1: f64,
    /// `cos(2 omega)` coefficient.
    pub c2: f64,
}

impl InversionParameters {
    /// No inversion term.
    #[must_use]
    pub fn none() -> Self {
        Self {
            force: 0.0,
            c0: 0.0,
            c1: 0.0,
            c2: 0.0,
        }
    }
}

/// Parameters for one three-coordinate centre.
///
/// `neighbours` are the UFF type names surrounding the centre; UFF4MOF's
/// `O_2` is the special case that raises an sp2 carbon's barrier.
#[must_use]
pub fn inversion_parameters(centre: &str, neighbours: &[&str]) -> InversionParameters {
    if is_sp2_centre(centre) {
        let total = if neighbours.contains(&"O_2") {
            SP2_INVERSION_WITH_O2
        } else {
            SP2_INVERSION
        };
        return InversionParameters {
            force: total / 3.0,
            c0: 1.0,
            c1: -1.0,
            c2: 0.0,
        };
    }
    if let Some((omega0, total)) = group_15(centre) {
        let omega = omega0.to_radians();
        let c2 = 1.0;
        let c1 = -4.0 * omega.cos();
        let c0 = -(c1 * omega.cos() + c2 * (2.0 * omega).cos());
        return InversionParameters {
            force: total / 3.0,
            c0,
            c1,
            c2,
        };
    }
    InversionParameters::none()
}

/// One out-of-plane term at `centre` with `j`, `k` in the plane and `l` out.
#[derive(Debug, Clone, Default)]
pub struct InversionTerm {
    pub(crate) centre: Vec<usize>,
    pub(crate) j: Vec<usize>,
    pub(crate) k: Vec<usize>,
    pub(crate) l: Vec<usize>,
    pub(crate) shift_j: Vec<[i32; 3]>,
    pub(crate) shift_k: Vec<[i32; 3]>,
    pub(crate) shift_l: Vec<[i32; 3]>,
    pub(crate) force: Vec<f64>,
    pub(crate) c0: Vec<f64>,
    pub(crate) c1: Vec<f64>,
    pub(crate) c2: Vec<f64>,
}

fn image_position(
    positions: &[[f64; 3]],
    index: usize,
    shift: [i32; 3],
    matrix: [[f64; 3]; 3],
) -> Result<[f64; 3], Error> {
    positions
        .get(index)
        .map(|position| add(*position, shift_vector(shift, matrix)))
        .ok_or_else(|| Error::Topology(format!("atom index {index} is out of range")))
}

impl InversionTerm {
    /// Number of inversions.
    #[must_use]
    pub fn len(&self) -> usize {
        self.centre.len()
    }

    /// Whether the term has no inversions.
    #[must_use]
    pub fn is_empty(&self) -> bool {
        self.centre.is_empty()
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
        let mut energy = 0.0;
        let mut virial_local = [[0.0; 3]; 3];
        for index in 0..self.len() {
            let (Some(centre), Some(atom_j), Some(atom_k), Some(atom_l)) = (
                self.centre.get(index),
                self.j.get(index),
                self.k.get(index),
                self.l.get(index),
            ) else {
                return Err(Error::Topology(
                    "inversion term arrays have inconsistent lengths".to_owned(),
                ));
            };
            let (Some(shift_j), Some(shift_k), Some(shift_l)) = (
                self.shift_j.get(index),
                self.shift_k.get(index),
                self.shift_l.get(index),
            ) else {
                return Err(Error::Topology(
                    "inversion shifts have inconsistent lengths".to_owned(),
                ));
            };
            let (Some(force), Some(c0), Some(c1), Some(c2)) = (
                self.force.get(index),
                self.c0.get(index),
                self.c1.get(index),
                self.c2.get(index),
            ) else {
                return Err(Error::Topology(
                    "inversion parameters have inconsistent lengths".to_owned(),
                ));
            };
            let origin = *positions
                .get(*centre)
                .ok_or_else(|| Error::Topology("inversion centre out of range".to_owned()))?;
            let a = sub(
                image_position(positions, *atom_j, *shift_j, matrix)?,
                origin,
            );
            let b = sub(
                image_position(positions, *atom_k, *shift_k, matrix)?,
                origin,
            );
            let c = sub(
                image_position(positions, *atom_l, *shift_l, matrix)?,
                origin,
            );

            let normal = cross(a, b);
            let nn = norm(normal);
            let nc = norm(c);
            if sine_from_cross(nn, a, b) <= MIN_SINE || nc <= EPS {
                continue;
            }
            let sin = (dot(normal, c) / (nn * nc)).clamp(-1.0, 1.0);
            let cos2 = (1.0 - sin * sin).clamp(0.0, 1.0);
            let cos = cos2.sqrt();
            energy += force * (c0 + c1 * cos + c2 * (2.0 * cos2 - 1.0));
            if cos <= EPS {
                continue;
            }
            let dsin = force * (-c1 * sin / cos - 4.0 * c2 * sin);

            let dn = sub(scale(c, 1.0 / (nn * nc)), scale(normal, sin / (nn * nn)));
            let dc = sub(scale(normal, 1.0 / (nn * nc)), scale(c, sin / (nc * nc)));
            let da = cross(b, dn);
            let db = cross(dn, a);
            let gradient_a = scale(da, dsin);
            let gradient_b = scale(db, dsin);
            let gradient_c = scale(dc, dsin);
            add_virial(&mut virial_local, a, gradient_a);
            add_virial(&mut virial_local, b, gradient_b);
            add_virial(&mut virial_local, c, gradient_c);
            scatter(&mut gradient, *atom_j, gradient_a);
            scatter(&mut gradient, *atom_k, gradient_b);
            scatter(&mut gradient, *atom_l, gradient_c);
            scatter(
                &mut gradient,
                *centre,
                negate(add(add(gradient_a, gradient_b), gradient_c)),
            );
        }
        if let Some(target) = virial {
            merge_virial(target, virial_local);
        }
        Ok((energy, gradient))
    }
}
