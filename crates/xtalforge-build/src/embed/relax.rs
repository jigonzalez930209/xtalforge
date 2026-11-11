//! Gradient-descent relaxation for the reduced distance-geometry embedder.

use super::geometry::unit;

/// Largest atom displacement in one relaxation step, Angstrom.
const MAX_STEP: f64 = 0.12;
/// Gradient descent step size.
const LEARNING_RATE: f64 = 0.35;
/// Stiffness of the bond-length term.
const BOND_STIFFNESS: f64 = 1.0;
/// Stiffness of the 1-3 (bond-angle) distance term.
const ANGLE_STIFFNESS: f64 = 0.4;
/// Stiffness of the non-bonded repulsion term.
const REPULSION_STIFFNESS: f64 = 0.35;

/// Gradient descent on bond lengths, bond angles and non-bonded repulsion.
pub(super) fn relax(
    points: &mut [[f64; 3]],
    bonds: &[(usize, usize, f64)],
    angles: &[(usize, usize, f64)],
    repulsion: &[(usize, usize, f64)],
    iterations: usize,
) {
    let mut forces = vec![[0.0; 3]; points.len()];
    for _ in 0..iterations {
        forces.fill([0.0; 3]);
        for &(first, second, target) in bonds {
            let (delta, distance) = pair_vector(points, first, second);
            let direction = unit(delta, first, second);
            let magnitude = BOND_STIFFNESS * (distance - target);
            apply_force(&mut forces, first, second, direction, magnitude);
        }
        for &(first, second, target) in angles {
            let (delta, distance) = pair_vector(points, first, second);
            let direction = unit(delta, first, second);
            let magnitude = ANGLE_STIFFNESS * (distance - target);
            apply_force(&mut forces, first, second, direction, magnitude);
        }
        for &(first, second, minimum) in repulsion {
            let (delta, distance) = pair_vector(points, first, second);
            if distance >= minimum {
                continue;
            }
            let direction = unit(delta, first, second);
            let magnitude = REPULSION_STIFFNESS * (distance - minimum);
            apply_force(&mut forces, first, second, direction, magnitude);
        }
        for (point, force) in points.iter_mut().zip(forces.iter()) {
            let [fx, fy, fz] = *force;
            let norm = (fx * fx + fy * fy + fz * fz).sqrt();
            let scale = if norm > MAX_STEP / LEARNING_RATE {
                MAX_STEP / LEARNING_RATE / norm
            } else {
                1.0
            };
            let [px, py, pz] = *point;
            *point = [
                px + LEARNING_RATE * fx * scale,
                py + LEARNING_RATE * fy * scale,
                pz + LEARNING_RATE * fz * scale,
            ];
        }
    }
}

fn pair_vector(points: &[[f64; 3]], first: usize, second: usize) -> ([f64; 3], f64) {
    let a = points.get(first).copied().unwrap_or([0.0; 3]);
    let b = points.get(second).copied().unwrap_or([0.0; 3]);
    let delta = [b[0] - a[0], b[1] - a[1], b[2] - a[2]];
    let distance = (delta[0] * delta[0] + delta[1] * delta[1] + delta[2] * delta[2]).sqrt();
    (delta, distance)
}

fn apply_force(
    forces: &mut [[f64; 3]],
    first: usize,
    second: usize,
    direction: [f64; 3],
    magnitude: f64,
) {
    let [dx, dy, dz] = direction;
    if let Some(force) = forces.get_mut(first) {
        let [fx, fy, fz] = *force;
        *force = [
            fx + dx * magnitude,
            fy + dy * magnitude,
            fz + dz * magnitude,
        ];
    }
    if let Some(force) = forces.get_mut(second) {
        let [fx, fy, fz] = *force;
        *force = [
            fx - dx * magnitude,
            fy - dy * magnitude,
            fz - dz * magnitude,
        ];
    }
}
