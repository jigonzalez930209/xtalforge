//! Ewald-summed EEQ2019 Coulomb matrix for periodic cells (audit 2026-11-11).
//!
//! The EEQ2019 interaction of two Gaussian charges is `erf(r / s_ij) / r` with
//! `s_ij = sqrt(rad_i^2 + rad_j^2)`. Summed over the images of a crystal it
//! converges only conditionally, and a plain real-space cutoff is not enough:
//! MOF-5 came out with Zn at +4.4. The periodic sum is split as
//!
//! ```text
//! sum_n' erf(r/s)/r = Ewald(1/r) - sum_n' erfc(r/s)/r,
//! ```
//!
//! where `Ewald(1/r)` is the standard Ewald potential of point charges with a
//! neutralising background (real space `erfc(a r)/r`, reciprocal space
//! `4 pi / V sum_k exp(-k^2 / 4 a^2) / k^2 cos(k . r)`, self term `-2 a /
//! sqrt(pi)`) and the second sum is short ranged. A constant background shift
//! does not change the charges, because the total-charge constraint absorbs
//! it. All quantities are in Bohr and Hartree.

#![allow(clippy::indexing_slicing)] // dense square matrix indexed by 0..n

use std::f64::consts::PI;

use crate::math::{cross, dot, norm};

/// Ewald splitting parameter in 1/Bohr.
const ALPHA: f64 = 0.3;
/// Real- and reciprocal-space accuracy target (`erfc(x) < 1e-8` at x = 4).
const REACH: f64 = 4.0;

/// Complementary error function (Numerical Recipes `erfcc`, 1.2e-7).
pub(crate) fn erfc(x: f64) -> f64 {
    let z = x.abs();
    let t = 1.0 / (1.0 + 0.5 * z);
    let value = t
        * (-z * z - 1.265_512_23
            + t * (1.000_023_68
                + t * (0.374_091_96
                    + t * (0.096_784_18
                        + t * (-0.186_288_06
                            + t * (0.278_868_07
                                + t * (-1.135_203_98
                                    + t * (1.488_515_87
                                        + t * (-0.822_152_23 + t * 0.170_872_77)))))))))
            .exp();
    if x >= 0.0 { value } else { 2.0 - value }
}

fn image_counts(lattice: [[f64; 3]; 3], radius: f64) -> [i32; 3] {
    let volume = dot(lattice[0], cross(lattice[1], lattice[2]))
        .abs()
        .max(1e-12);
    let widths = [
        volume / norm(cross(lattice[1], lattice[2])).max(1e-12),
        volume / norm(cross(lattice[2], lattice[0])).max(1e-12),
        volume / norm(cross(lattice[0], lattice[1])).max(1e-12),
    ];
    widths.map(|width| (radius / width).ceil() as i32 + 1)
}

fn combine(n: [i32; 3], lattice: [[f64; 3]; 3]) -> [f64; 3] {
    let [a, b, c] = n.map(f64::from);
    [0, 1, 2].map(|k| a * lattice[0][k] + b * lattice[1][k] + c * lattice[2][k])
}

/// The periodic `sum_n' erf(r / s_ij) / r` matrix, without the EEQ diagonal
/// self term, for Cartesian positions and a lattice in Bohr.
pub(crate) fn periodic_matrix(
    lattice: [[f64; 3]; 3],
    positions: &[[f64; 3]],
    widths: &[f64],
    dielectric: f64,
) -> Vec<Vec<f64>> {
    let n = positions.len();
    let mut matrix = vec![vec![0.0; n]; n];
    let volume = dot(lattice[0], cross(lattice[1], lattice[2])).abs();
    let largest_width = widths.iter().copied().fold(0.0_f64, f64::max);
    let real_cutoff = (REACH / ALPHA).max(REACH * largest_width * std::f64::consts::SQRT_2);
    let real = image_counts(lattice, real_cutoff);
    for na in -real[0]..=real[0] {
        for nb in -real[1]..=real[1] {
            for nc in -real[2]..=real[2] {
                let shift = combine([na, nb, nc], lattice);
                let home = na == 0 && nb == 0 && nc == 0;
                for i in 0..n {
                    for j in 0..n {
                        if home && i == j {
                            continue;
                        }
                        let d = [0, 1, 2].map(|k| positions[j][k] + shift[k] - positions[i][k]);
                        let r = norm(d);
                        if r > real_cutoff || r < 1e-10 {
                            continue;
                        }
                        let s = (widths[i] * widths[i] + widths[j] * widths[j]).sqrt();
                        matrix[i][j] += (erfc(ALPHA * r) - erfc(r / s)) / r;
                    }
                }
            }
        }
    }
    // Reciprocal space from the reciprocal basis b_i = 2 pi (a_j x a_k) / V.
    let reciprocal = [
        cross(lattice[1], lattice[2]),
        cross(lattice[2], lattice[0]),
        cross(lattice[0], lattice[1]),
    ]
    .map(|row| row.map(|value| 2.0 * PI * value / volume));
    let k_cutoff = 2.0 * ALPHA * REACH;
    let counts = image_counts(reciprocal, k_cutoff);
    let mut cosines = vec![0.0; n];
    let mut sines = vec![0.0; n];
    for ka in -counts[0]..=counts[0] {
        for kb in -counts[1]..=counts[1] {
            for kc in -counts[2]..=counts[2] {
                if ka == 0 && kb == 0 && kc == 0 {
                    continue;
                }
                let k = combine([ka, kb, kc], reciprocal);
                let k2 = dot(k, k);
                if k2 > k_cutoff * k_cutoff {
                    continue;
                }
                let factor = 4.0 * PI / volume * (-k2 / (4.0 * ALPHA * ALPHA)).exp() / k2;
                for (index, position) in positions.iter().enumerate() {
                    let phase = dot(k, *position);
                    cosines[index] = phase.cos();
                    sines[index] = phase.sin();
                }
                for i in 0..n {
                    for j in 0..n {
                        matrix[i][j] += factor * (cosines[i] * cosines[j] + sines[i] * sines[j]);
                    }
                }
            }
        }
    }
    let self_term = -2.0 * ALPHA / PI.sqrt();
    for (i, row) in matrix.iter_mut().enumerate() {
        row[i] += self_term;
        for value in row.iter_mut() {
            *value /= dielectric;
        }
    }
    matrix
}
