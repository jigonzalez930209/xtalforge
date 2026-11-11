//! Kabsch orientation of a node block onto the edges of a vertex.
//!
//! The greedy orientation aligns only the mean connection direction with the
//! mean edge direction, which leaves the in-plane angle of a planar node (or
//! the twist of an octahedral one) arbitrary; the linkers then have to be
//! stretched or squeezed to close the gaps. Here every one-to-one assignment
//! of connection points to edges is tried (up to six connections, 720
//! permutations), the proper rotation that best maps the unit connection
//! directions onto the unit edge directions is found with Horn's quaternion
//! form of the Kabsch fit (the dominant eigenvector of a 4 x 4 matrix), and
//! the assignment with the smallest residual wins. Ties keep the first
//! permutation in lexicographic order, so the result is deterministic.

#![allow(clippy::indexing_slicing)] // permutation indices are bounded by the slice length

use crate::orient::{Slot, normalize};

/// Largest connection count searched exhaustively.
const MAX_CONNECTIONS: usize = 6;

/// A rotation and the connection assigned to every slot.
type Orientation = ([[f64; 3]; 3], Vec<Option<usize>>);

/// Best rotation and the connection assigned to every slot.
pub(crate) fn kabsch_orientation(connections: &[[f64; 3]], slots: &[Slot]) -> Option<Orientation> {
    let count = connections.len();
    if count == 0 || count != slots.len() || count > MAX_CONNECTIONS {
        return None;
    }
    let sources: Vec<[f64; 3]> = connections
        .iter()
        .map(|point| normalize(*point))
        .collect::<Option<_>>()?;
    let targets: Vec<[f64; 3]> = slots
        .iter()
        .map(|slot| normalize(slot.direction))
        .collect::<Option<_>>()?;
    let mut permutation: Vec<usize> = (0..count).collect();
    let mut best: Option<(f64, [[f64; 3]; 3], Vec<usize>)> = None;
    loop {
        let rotation = horn(&sources, &targets, &permutation);
        let residual: f64 = permutation
            .iter()
            .zip(&targets)
            .filter_map(|(source, target)| {
                sources.get(*source).map(|value| {
                    let moved = apply(rotation, *value);
                    (0..3)
                        .map(|axis| (moved[axis] - target[axis]).powi(2))
                        .sum::<f64>()
                })
            })
            .sum();
        if best
            .as_ref()
            .is_none_or(|(value, _, _)| residual < *value - 1e-12)
        {
            best = Some((residual, rotation, permutation.clone()));
        }
        if !next_permutation(&mut permutation) {
            break;
        }
    }
    let (_, rotation, permutation) = best?;
    // Round-off would otherwise leave 0.99999999 and 1e-17 in the axis
    // aligned fits; exact entries keep the cell scale and coordinates clean.
    let matrix = rotation.map(|row| {
        row.map(|value| {
            let nearest = value.round();
            if (value - nearest).abs() < 1e-9 {
                nearest
            } else {
                value
            }
        })
    });
    Some((matrix, permutation.into_iter().map(Some).collect()))
}

fn apply(matrix: [[f64; 3]; 3], vector: [f64; 3]) -> [f64; 3] {
    matrix.map(|row| row[0] * vector[0] + row[1] * vector[1] + row[2] * vector[2])
}

/// Proper rotation minimising `sum |R source[p(i)] - target[i]|^2`.
///
/// Horn (1987): the optimal unit quaternion is the eigenvector of the
/// largest eigenvalue of the symmetric 4 x 4 matrix built from the
/// correlation `S = sum source target^T`. The eigenvector is found by power
/// iteration on `N + c I`, with `c` large enough to make every eigenvalue
/// positive.
fn horn(sources: &[[f64; 3]], targets: &[[f64; 3]], permutation: &[usize]) -> [[f64; 3]; 3] {
    let mut s = [[0.0; 3]; 3];
    for (target, source) in targets.iter().zip(permutation) {
        if let Some(source) = sources.get(*source) {
            for a in 0..3 {
                for b in 0..3 {
                    s[a][b] += source[a] * target[b];
                }
            }
        }
    }
    let [[sxx, sxy, sxz], [syx, syy, syz], [szx, szy, szz]] = s;
    let n = [
        [sxx + syy + szz, syz - szy, szx - sxz, sxy - syx],
        [syz - szy, sxx - syy - szz, sxy + syx, szx + sxz],
        [szx - sxz, sxy + syx, -sxx + syy - szz, syz + szy],
        [sxy - syx, szx + sxz, syz + szy, -sxx - syy + szz],
    ];
    let shift = n.iter().flatten().map(|value| value.abs()).sum::<f64>() + 1.0;
    // Start away from any special axis so no eigenvector is missed.
    let mut q = [1.0, 0.31, 0.17, 0.07];
    for _ in 0..500 {
        let mut next = [0.0; 4];
        for (row, value) in next.iter_mut().enumerate() {
            *value = shift * q[row] + (0..4).map(|column| n[row][column] * q[column]).sum::<f64>();
        }
        let length = next.iter().map(|value| value * value).sum::<f64>().sqrt();
        if length < 1e-300 {
            break;
        }
        let next = next.map(|value| value / length);
        let change: f64 = (0..4).map(|index| (next[index] - q[index]).abs()).sum();
        q = next;
        if change < 1e-15 {
            break;
        }
    }
    let [w, x, y, z] = q;
    [
        [
            w * w + x * x - y * y - z * z,
            2.0 * (x * y - w * z),
            2.0 * (x * z + w * y),
        ],
        [
            2.0 * (x * y + w * z),
            w * w - x * x + y * y - z * z,
            2.0 * (y * z - w * x),
        ],
        [
            2.0 * (x * z - w * y),
            2.0 * (y * z + w * x),
            w * w - x * x - y * y + z * z,
        ],
    ]
}

/// Lexicographic next permutation; false after the last one.
fn next_permutation(values: &mut [usize]) -> bool {
    let Some(pivot) = (1..values.len())
        .rev()
        .find(|&index| values[index - 1] < values[index])
    else {
        return false;
    };
    let Some(successor) = (pivot..values.len())
        .rev()
        .find(|&index| values[index] > values[pivot - 1])
    else {
        return false;
    };
    values.swap(pivot - 1, successor);
    values[pivot..].reverse();
    true
}
