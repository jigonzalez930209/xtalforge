//! UFF torsion term, eq 15-18, and its analytic gradient (Task 4.2.3).
//!
//! Ported from Crystal-Builder `xtal/ff/uff/terms.py`. The energy is written
//! in `cos(phi)` through the Chebyshev expansion so a torsion passing through
//! planarity is an ordinary point. A torsion whose central bonds are nearly
//! collinear has no defined dihedral; it is dropped at `MIN_SINE` rather than
//! producing an infinite derivative.

use crate::Error;
use crate::uff::math::{
    MIN_SINE, add, cross, dot, merge_virial, negate, norm, scale, shift_vector, sine_from_cross,
    sub,
};
use crate::uff::parallel::{Part, evaluate};
use crate::uff::params::get;

const GROUP_16_BARRIER_HEAVY: f64 = 6.8;

fn group_16(element: &str) -> bool {
    matches!(element, "O" | "S" | "Se" | "Te" | "Po")
}

fn group_16_barrier(element: &str) -> f64 {
    if element == "O" {
        2.0
    } else {
        GROUP_16_BARRIER_HEAVY
    }
}

/// Hybridisation implied by a UFF type's geometry character.
///
/// # Errors
///
/// Returns [`Error::UnknownType`] for an unknown type name.
pub fn hybridisation(type_name: &str) -> Result<&'static str, Error> {
    Ok(match get(type_name)?.geometry_char() {
        Some('1') => "sp",
        Some('2' | 'R') => "sp2",
        Some('3') => "sp3",
        _ => "",
    })
}

/// The torsion a central bond contributes, eq 16-18.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct TorsionParameters {
    /// Barrier `V` in kcal/mol; zero means no torsion about this bond.
    pub barrier: f64,
    /// Periodicity `n`: 2, 3 or 6.
    pub periodicity: i32,
    /// `cos(n phi0)`, either `+1` or `-1`.
    pub cos_n_phi0: f64,
}

impl TorsionParameters {
    /// No torsion.
    #[must_use]
    pub fn none() -> Self {
        Self {
            barrier: 0.0,
            periodicity: 1,
            cos_n_phi0: 1.0,
        }
    }
}

/// Torsion parameters for the bond between `type_j` and `type_k`.
///
/// Both atoms must be sp2 or sp3 and carry at least two neighbours. The
/// barrier is computed for the bond and later split between the torsions
/// around it; a metal on either end disables the term, because a dihedral
/// about a bond to a six-coordinate centre says nothing and N-Zn-N at 180
/// degrees is exactly the geometry whose torsion gradient is undefined.
///
/// # Errors
///
/// Returns [`Error::UnknownType`] for an unknown type name.
pub fn torsion_parameters(
    type_j: &str,
    type_k: &str,
    order: f64,
    neighbors_j: usize,
    neighbors_k: usize,
) -> Result<TorsionParameters, Error> {
    let hybrid_j = hybridisation(type_j)?;
    let hybrid_k = hybridisation(type_k)?;
    if !matches!(hybrid_j, "sp2" | "sp3") || !matches!(hybrid_k, "sp2" | "sp3") {
        return Ok(TorsionParameters::none());
    }
    if neighbors_j < 2 || neighbors_k < 2 {
        return Ok(TorsionParameters::none());
    }
    let params_j = get(type_j)?;
    let params_k = get(type_k)?;
    if params_j.is_metal() || params_k.is_metal() {
        return Ok(TorsionParameters::none());
    }
    let order = order.max(1e-6);
    let sp2_barrier = 5.0 * (params_j.uj * params_k.uj).sqrt() * (1.0 + 4.18 * order.ln());
    if hybrid_j == "sp3" && hybrid_k == "sp3" {
        if group_16(&params_j.element) && group_16(&params_k.element) {
            return Ok(TorsionParameters {
                barrier: (group_16_barrier(&params_j.element)
                    * group_16_barrier(&params_k.element))
                .sqrt(),
                periodicity: 2,
                cos_n_phi0: -1.0,
            });
        }
        return Ok(TorsionParameters {
            barrier: (params_j.vi * params_k.vi).sqrt(),
            periodicity: 3,
            cos_n_phi0: -1.0,
        });
    }
    if hybrid_j == "sp2" && hybrid_k == "sp2" {
        return Ok(TorsionParameters {
            barrier: sp2_barrier,
            periodicity: 2,
            cos_n_phi0: 1.0,
        });
    }
    let sp3_element = if hybrid_j == "sp3" {
        &params_j.element
    } else {
        &params_k.element
    };
    if group_16(sp3_element) {
        return Ok(TorsionParameters {
            barrier: sp2_barrier,
            periodicity: 2,
            cos_n_phi0: -1.0,
        });
    }
    Ok(TorsionParameters {
        barrier: 1.0,
        periodicity: 6,
        cos_n_phi0: 1.0,
    })
}

/// One dihedral `i-j-k-l`, with every image relative to `j`.
#[derive(Debug, Clone, Default)]
pub struct TorsionTerm {
    pub(crate) i: Vec<usize>,
    pub(crate) j: Vec<usize>,
    pub(crate) k: Vec<usize>,
    pub(crate) l: Vec<usize>,
    pub(crate) shift_i: Vec<[i32; 3]>,
    pub(crate) shift_k: Vec<[i32; 3]>,
    pub(crate) shift_l: Vec<[i32; 3]>,
    pub(crate) barrier: Vec<f64>,
    pub(crate) periodicity: Vec<i32>,
    pub(crate) cos_n_phi0: Vec<f64>,
}

fn term_vector(
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

fn chebyshev(cos: f64, order: i32) -> Option<(f64, f64)> {
    match order {
        2 => Some((2.0 * cos * cos - 1.0, 4.0 * cos)),
        3 => Some((4.0 * cos.powi(3) - 3.0 * cos, 12.0 * cos * cos - 3.0)),
        6 => Some((
            32.0 * cos.powi(6) - 48.0 * cos.powi(4) + 18.0 * cos * cos - 1.0,
            192.0 * cos.powi(5) - 192.0 * cos.powi(3) + 36.0 * cos,
        )),
        _ => None,
    }
}

impl TorsionTerm {
    /// Number of torsions.
    #[must_use]
    pub fn len(&self) -> usize {
        self.j.len()
    }

    /// Whether the term has no torsions.
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
        self.energy_and_gradient_with(positions, matrix, virial, true)
    }

    pub(crate) fn energy_and_gradient_with(
        &self,
        positions: &[[f64; 3]],
        matrix: [[f64; 3]; 3],
        virial: Option<&mut [[f64; 3]; 3]>,
        parallel: bool,
    ) -> Result<(f64, Vec<[f64; 3]>), Error> {
        let mut gradient = vec![[0.0; 3]; positions.len()];
        let mut energy = 0.0;
        let mut virial_local = [[0.0; 3]; 3];
        evaluate(
            self.len(),
            parallel,
            |index| self.part(index, positions, matrix),
            &mut energy,
            &mut gradient,
            &mut virial_local,
        )?;
        if let Some(target) = virial {
            merge_virial(target, virial_local);
        }
        Ok((energy, gradient))
    }

    fn part(
        &self,
        index: usize,
        positions: &[[f64; 3]],
        matrix: [[f64; 3]; 3],
    ) -> Result<Part, Error> {
        let mut part = Part::zero();
        let (Some(atom_i), Some(atom_j), Some(atom_k), Some(atom_l)) = (
            self.i.get(index),
            self.j.get(index),
            self.k.get(index),
            self.l.get(index),
        ) else {
            return Err(Error::Topology(
                "torsion term arrays have inconsistent lengths".to_owned(),
            ));
        };
        let (Some(shift_i), Some(shift_k), Some(shift_l)) = (
            self.shift_i.get(index),
            self.shift_k.get(index),
            self.shift_l.get(index),
        ) else {
            return Err(Error::Topology(
                "torsion shifts have inconsistent lengths".to_owned(),
            ));
        };
        let (Some(barrier), Some(periodicity), Some(cos_n_phi0)) = (
            self.barrier.get(index),
            self.periodicity.get(index),
            self.cos_n_phi0.get(index),
        ) else {
            return Err(Error::Topology(
                "torsion parameters have inconsistent lengths".to_owned(),
            ));
        };
        let position_i = term_vector(positions, *atom_i, *shift_i, matrix)?;
        let position_j = term_vector(positions, *atom_j, [0, 0, 0], matrix)?;
        let position_k = term_vector(positions, *atom_k, *shift_k, matrix)?;
        let position_l = term_vector(positions, *atom_l, *shift_l, matrix)?;
        let b1 = sub(position_j, position_i);
        let b2 = sub(position_k, position_j);
        let b3 = sub(position_l, position_k);

        let n1 = cross(b1, b2);
        let n2 = cross(b2, b3);
        let m1 = norm(n1);
        let m2 = norm(n2);
        if sine_from_cross(m1, b1, b2) <= MIN_SINE || sine_from_cross(m2, b2, b3) <= MIN_SINE {
            return Ok(part);
        }
        let cos = (dot(n1, n2) / (m1 * m2)).clamp(-1.0, 1.0);
        let Some((value, slope)) = chebyshev(cos, *periodicity) else {
            return Err(Error::Topology(format!(
                "unsupported torsion periodicity {periodicity}"
            )));
        };
        let amplitude = 0.5 * barrier * cos_n_phi0;
        part.energy = 0.5 * barrier - amplitude * value;
        let dcos = -amplitude * slope;

        let w1 = sub(scale(n2, 1.0 / (m1 * m2)), scale(n1, cos / (m1 * m1)));
        let w2 = sub(scale(n1, 1.0 / (m1 * m2)), scale(n2, cos / (m2 * m2)));
        let db1 = cross(b2, w1);
        let db2 = add(cross(w1, b1), cross(b3, w2));
        let db3 = cross(w2, b2);

        let gradient_i = scale(db1, -dcos);
        let gradient_j = scale(sub(db1, db2), dcos);
        let gradient_k = scale(sub(db2, db3), dcos);
        let gradient_l = scale(db3, dcos);
        part.push_virial(b1, negate(gradient_i));
        part.push_virial(b2, scale(db2, dcos));
        part.push_virial(b3, gradient_l);
        part.push_scatter(*atom_i, gradient_i);
        part.push_scatter(*atom_j, gradient_j);
        part.push_scatter(*atom_k, gradient_k);
        part.push_scatter(*atom_l, gradient_l);
        Ok(part)
    }
}
