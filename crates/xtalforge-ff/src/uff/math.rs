//! Small vector helpers shared by the UFF terms.
//!
//! Everything is a plain `[f64; 3]`; the crate deliberately avoids a linear
//! algebra dependency for three-vectors. Functions are `pub(crate)` because
//! the public surface is the force field, not its arithmetic.

/// Below this length a vector has no direction and the term is skipped.
pub(crate) const EPS: f64 = 1e-9;
/// Minimum sine for a torsion or inversion; see `terms.py` and Task 4.2.3.
pub(crate) const MIN_SINE: f64 = 1e-2;

pub(crate) fn dot(first: [f64; 3], second: [f64; 3]) -> f64 {
    first.iter().zip(second.iter()).map(|(a, b)| a * b).sum()
}

pub(crate) fn cross(first: [f64; 3], second: [f64; 3]) -> [f64; 3] {
    let [ax, ay, az] = first;
    let [bx, by, bz] = second;
    [ay * bz - az * by, az * bx - ax * bz, ax * by - ay * bx]
}

pub(crate) fn norm(vector: [f64; 3]) -> f64 {
    dot(vector, vector).sqrt()
}

pub(crate) fn sub(first: [f64; 3], second: [f64; 3]) -> [f64; 3] {
    let [ax, ay, az] = first;
    let [bx, by, bz] = second;
    [ax - bx, ay - by, az - bz]
}

pub(crate) fn add(first: [f64; 3], second: [f64; 3]) -> [f64; 3] {
    let [ax, ay, az] = first;
    let [bx, by, bz] = second;
    [ax + bx, ay + by, az + bz]
}

pub(crate) fn scale(vector: [f64; 3], factor: f64) -> [f64; 3] {
    let [x, y, z] = vector;
    [x * factor, y * factor, z * factor]
}

pub(crate) fn negate(vector: [f64; 3]) -> [f64; 3] {
    let [x, y, z] = vector;
    [-x, -y, -z]
}

/// Cartesian vector of an integer lattice translation, `shift @ matrix`.
pub(crate) fn shift_vector(shift: [i32; 3], matrix: [[f64; 3]; 3]) -> [f64; 3] {
    let [[a, b, c], [d, e, f], [g, h, i]] = matrix;
    let [x, y, z] = shift.map(f64::from);
    [
        x * a + y * d + z * g,
        x * b + y * e + z * h,
        x * c + y * f + z * i,
    ]
}

/// Angle between two vectors in radians, clamped at the domain edges.
pub(crate) fn angle_between(first: [f64; 3], second: [f64; 3]) -> f64 {
    let lengths = norm(first) * norm(second);
    if lengths <= EPS {
        return 0.0;
    }
    (dot(first, second) / lengths).clamp(-1.0, 1.0).acos()
}

/// `sin(angle)` from a cross-product magnitude already taken.
pub(crate) fn sine_from_cross(cross_norm: f64, first: [f64; 3], second: [f64; 3]) -> f64 {
    let lengths = norm(first) * norm(second);
    cross_norm / lengths.max(EPS)
}

/// Adds `sum d (x) g` to a virial tensor.
pub(crate) fn add_virial(virial: &mut [[f64; 3]; 3], vectors: [f64; 3], gradients: [f64; 3]) {
    for row in 0..3 {
        for column in 0..3 {
            if let (Some(slot), Some(vector), Some(gradient)) = (
                virial
                    .get_mut(row)
                    .and_then(|values| values.get_mut(column)),
                vectors.get(row),
                gradients.get(column),
            ) {
                *slot += vector * gradient;
            }
        }
    }
}

/// Adds one local virial tensor into another.
pub(crate) fn merge_virial(target: &mut [[f64; 3]; 3], source: [[f64; 3]; 3]) {
    for row in 0..3 {
        for column in 0..3 {
            if let (Some(slot), Some(value)) = (
                target
                    .get_mut(row)
                    .and_then(|values| values.get_mut(column)),
                source.get(row).and_then(|values| values.get(column)),
            ) {
                *slot += value;
            }
        }
    }
}

/// `grad[indices] += values`, accumulating repeated indices.
pub(crate) fn scatter(forces: &mut [[f64; 3]], index: usize, value: [f64; 3]) {
    if let Some(slot) = forces.get_mut(index) {
        *slot = add(*slot, value);
    }
}

/// Symmetrised stress from a virial in kcal/mol/Angstrom^3.
pub(crate) fn virial_to_stress(virial: [[f64; 3]; 3], volume: f64) -> [[f64; 3]; 3] {
    let mut stress = [[0.0; 3]; 3];
    for row in 0..3 {
        for column in 0..3 {
            let forward = virial
                .get(row)
                .and_then(|values| values.get(column))
                .copied()
                .unwrap_or(0.0);
            let reverse = virial
                .get(column)
                .and_then(|values| values.get(row))
                .copied()
                .unwrap_or(0.0);
            let value =
                0.5 * (forward + reverse) / volume * crate::uff::params::KCAL_PER_MOL_A3_TO_GPA;
            if let Some(slot) = stress
                .get_mut(row)
                .and_then(|values| values.get_mut(column))
            {
                *slot = value;
            }
        }
    }
    stress
}
