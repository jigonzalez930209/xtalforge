//! Commensurate twisted graphene multilayers (Task 6.5.3).
//!
//! A twist is commensurate when a superlattice vector of the bottom layer,
//! `L1 = i a1 + j a2`, is also a lattice vector of the rotated top layer. With
//! `a1 = (a, 0)` and `a2 = (a/2, a sqrt(3)/2)` the rotation that maps
//! `j a1 + i a2` onto `i a1 + j a2` has
//!
//! ```text
//! cos(theta) = (i^2 + 4 i j + j^2) / (2 (i^2 + i j + j^2))
//! ```
//!
//! and the moire cell spanned by `L1` and `L2 = R(60) L1` holds
//! `2 (i^2 + i j + j^2)` atoms per layer (Lopes dos Santos et al., Phys. Rev.
//! Lett. 99, 256802, 2007). `(i, j) = (2, 1)` gives 21.787 degrees and 14 atoms
//! per layer.

use xtalforge_core::{Error, Lattice, Site, Structure, StructureBuilder};

use crate::carbon::A_GRAPHENE;

/// Interlayer spacing of the stacked layers in Angstrom.
const INTERLAYER: f64 = 3.35;
/// Vacuum added above the stack in Angstrom.
const VACUUM: f64 = 15.0;
/// Largest moire index searched; `(20, 19)` already holds 2 282 atoms per layer.
const MAX_INDEX: u32 = 20;
/// Angle error, in degrees, accepted for the smallest matching cell.
const ANGLE_TOLERANCE: f64 = 0.05;

/// A commensurate moire cell.
#[derive(Debug, Clone, Copy)]
struct Moire {
    i: u32,
    j: u32,
    /// Twist angle in degrees, in `[0, 60]`.
    angle: f64,
}

impl Moire {
    fn atoms_per_layer(self) -> u32 {
        2 * (self.i * self.i + self.i * self.j + self.j * self.j)
    }
}

fn gcd(mut a: u32, mut b: u32) -> u32 {
    while b != 0 {
        (a, b) = (b, a % b);
    }
    a
}

fn moire_angle(i: u32, j: u32) -> f64 {
    let (i, j) = (f64::from(i), f64::from(j));
    let cosine = (i * i + 4.0 * i * j + j * j) / (2.0 * (i * i + i * j + j * j));
    cosine.clamp(-1.0, 1.0).acos().to_degrees()
}

/// Smallest commensurate cell within [`ANGLE_TOLERANCE`], else the closest one.
fn find_moire(theta: f64) -> Moire {
    // Rotations by 120 degrees and mirrors are lattice symmetries, so the
    // twist reduces to [0, 60]; 0 (AA) and 60 degrees (AA') stay distinct.
    let reduced = theta.rem_euclid(120.0);
    let target = if reduced > 60.0 {
        120.0 - reduced
    } else {
        reduced
    };
    if target <= ANGLE_TOLERANCE {
        return Moire {
            i: 1,
            j: 0,
            angle: 0.0,
        };
    }
    let mut candidates = Vec::new();
    for i in 0..=MAX_INDEX {
        for j in 0..=MAX_INDEX {
            if (i == 0 && j == 0) || gcd(i, j) != 1 {
                continue;
            }
            let angle = moire_angle(i, j);
            let error = (angle - target).abs();
            candidates.push((Moire { i, j, angle }, error));
        }
    }
    candidates.sort_by(|(a, ea), (b, eb)| {
        let within_a = *ea <= ANGLE_TOLERANCE;
        let within_b = *eb <= ANGLE_TOLERANCE;
        within_b
            .cmp(&within_a)
            .then_with(|| {
                if within_a && within_b {
                    a.atoms_per_layer().cmp(&b.atoms_per_layer())
                } else {
                    ea.total_cmp(eb)
                }
            })
            .then_with(|| ea.total_cmp(eb))
            .then_with(|| (b.i, b.j).cmp(&(a.i, a.j)))
    });
    candidates.first().map_or(
        Moire {
            i: 1,
            j: 0,
            angle: 0.0,
        },
        |(moire, _)| *moire,
    )
}

fn rotate(v: [f64; 2], angle: f64) -> [f64; 2] {
    let (s, c) = angle.sin_cos();
    [c * v[0] - s * v[1], s * v[0] + c * v[1]]
}

/// Fractional coordinates (in the moire basis) of one rotated honeycomb layer.
fn layer_positions(l1: [f64; 2], l2: [f64; 2], rotation: f64, expected: usize) -> Vec<[f64; 2]> {
    let a1 = rotate([A_GRAPHENE, 0.0], rotation);
    let a2 = rotate([A_GRAPHENE / 2.0, A_GRAPHENE * 3f64.sqrt() / 2.0], rotation);
    let det = l1[0] * l2[1] - l1[1] * l2[0];
    let span = i64::from(MAX_INDEX) * 3;
    let mut positions: Vec<[f64; 2]> = Vec::with_capacity(expected);
    for p in -span..=span {
        for q in -span..=span {
            for (u, v) in [(0.0, 0.0), (1.0 / 3.0, 1.0 / 3.0)] {
                let (pu, qv) = (p as f64 + u, q as f64 + v);
                let x = pu * a1[0] + qv * a2[0];
                let y = pu * a1[1] + qv * a2[1];
                let fa = ((x * l2[1] - y * l2[0]) / det).rem_euclid(1.0);
                let fb = ((l1[0] * y - l1[1] * x) / det).rem_euclid(1.0);
                let wrap = |f: f64| if f > 1.0 - 1e-6 { 0.0 } else { f };
                let point = [wrap(fa), wrap(fb)];
                let duplicate = positions.iter().any(|seen| {
                    let d0 = (seen[0] - point[0]).abs();
                    let d1 = (seen[1] - point[1]).abs();
                    d0.min(1.0 - d0) < 1e-6 && d1.min(1.0 - d1) < 1e-6
                });
                if !duplicate {
                    positions.push(point);
                }
            }
        }
    }
    positions.sort_by(|a, b| a[1].total_cmp(&b[1]).then(a[0].total_cmp(&b[0])));
    positions
}

/// Builds a commensurate twisted graphene multilayer.
///
/// The moire cell is the smallest commensurate cell whose twist angle lies
/// within 0.05 degrees of `theta_degrees` (indices up to 20); when none does,
/// the closest angle is used. The angle is first reduced to `[0, 60]` with the
/// lattice symmetries; 0 degrees is AA stacking in the two-atom cell. Even layers keep the reference orientation and
/// odd layers are rotated by the twist, 3.35 Angstrom apart with 15 Angstrom of
/// vacuum. The moire cell is then repeated `n x m` times in plane. The
/// metadata records `twist_requested`, the exact `twist_angle`, the moire
/// indices `twist_indices` and, for compatibility, `twist_approximation`
/// (equal to `twist_angle`).
///
/// # Errors
///
/// Returns [`Error::InvalidInput`] for zero repeats/layers or a non-finite
/// angle.
pub fn commensurate_twist(
    n: u32,
    m: u32,
    theta_degrees: f64,
    layers: u32,
) -> Result<Structure, Error> {
    if n == 0 || m == 0 || layers == 0 || !theta_degrees.is_finite() {
        return Err(Error::InvalidInput(
            "twist repeats, layers and angle must be valid".to_owned(),
        ));
    }
    let moire = find_moire(theta_degrees);
    let (fi, fj) = (f64::from(moire.i), f64::from(moire.j));
    let a1 = [A_GRAPHENE, 0.0];
    let a2 = [A_GRAPHENE / 2.0, A_GRAPHENE * 3f64.sqrt() / 2.0];
    let l1 = [fi * a1[0] + fj * a2[0], fi * a1[1] + fj * a2[1]];
    let l2 = rotate(l1, std::f64::consts::FRAC_PI_3);
    // Rotation taking j a1 + i a2 onto i a1 + j a2.
    let source = [fj * a1[0] + fi * a2[0], fj * a1[1] + fi * a2[1]];
    let rotation = if moire.angle > 0.0 {
        l1[1].atan2(l1[0]) - source[1].atan2(source[0])
    } else {
        0.0
    };
    let expected = moire.atoms_per_layer() as usize;
    let bottom = layer_positions(l1, l2, 0.0, expected);
    let top = layer_positions(l1, l2, rotation, expected);
    if bottom.len() != expected || top.len() != expected {
        return Err(Error::InvalidInput(format!(
            "twist cell ({}, {}) produced {} and {} atoms instead of {expected}",
            moire.i,
            moire.j,
            bottom.len(),
            top.len()
        )));
    }
    let height = f64::from(layers) * INTERLAYER + VACUUM;
    let (fn_, fm) = (f64::from(n), f64::from(m));
    let lattice = Lattice::new([
        [l1[0] * fn_, l1[1] * fn_, 0.0],
        [l2[0] * fm, l2[1] * fm, 0.0],
        [0.0, 0.0, height],
    ])?;
    let mut builder = StructureBuilder::new("twisted-graphene", lattice)
        .metadata("twist_requested", serde_json::json!(theta_degrees))
        .metadata("twist_angle", serde_json::json!(moire.angle))
        .metadata("twist_indices", serde_json::json!([moire.i, moire.j]))
        .metadata("twist_approximation", serde_json::json!(moire.angle));
    let mut index = 0usize;
    for layer in 0..layers {
        let z = (f64::from(layer) * INTERLAYER + 1.0) / height;
        let positions = if layer % 2 == 0 { &bottom } else { &top };
        for ra in 0..n {
            for rb in 0..m {
                for [x, y] in positions {
                    index += 1;
                    let frac = [(x + f64::from(ra)) / fn_, (y + f64::from(rb)) / fm, z];
                    builder = builder.site(Site::new(format!("C{index}"), "C", frac)?);
                }
            }
        }
    }
    builder.build()
}
