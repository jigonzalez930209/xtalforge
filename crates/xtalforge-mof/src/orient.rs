//! Deterministic node orientation (Task 6.13.1 subset).
//!
//! PORMAKE fits a node's connection directions onto a slot's directions by
//! RMSD over all rotations. This crate uses a documented deterministic
//! subset:
//!
//! 1. The local reference direction is the normalized mean of the node's
//!    connection-point vectors. When that mean is zero (an octahedral or
//!    trigonal node), the lowest-index connection point supplies the
//!    reference.
//! 2. The target direction is the normalized mean of the incident edge
//!    directions. When that mean is zero, the lowest-index edge supplies the
//!    target.
//! 3. The rotation is the minimal proper rotation between those directions.
//!    The antiparallel case uses a fixed perpendicular axis, so it is
//!    deterministic.
//! 4. Connection points are then matched to edges greedily: the pair with the
//!    largest cosine similarity is assigned first, with ties broken by edge
//!    index and connection-point index. Only positive similarities count.
//!
//! This is not the RMSD-optimal PORMAKE locator and does not explore the
//! symmetry-equivalent fits of a highly symmetric node; the public
//! [`crate::node_oriented`] helper exposes the rotation for tests.

use crate::{Block, Error};

const EPSILON: f64 = 1e-9;

/// One incident edge at a vertex.
#[derive(Debug, Clone, Copy)]
pub(crate) struct Slot {
    /// Index into the topology edge list.
    pub edge: usize,
    /// `0` when the vertex is the first endpoint, `1` when it is the second.
    pub side: usize,
    /// Unit direction from the vertex towards the other endpoint.
    pub direction: [f64; 3],
}

/// Rotation matrix that maps the block's reference connection direction onto
/// `target_direction`.
///
/// The reference is the normalized mean of the block's connection-point
/// vectors; when that mean is zero, the first connection point is used. The
/// result is the minimal proper rotation, deterministic for every input.
///
/// # Errors
///
/// Returns [`Error::InvalidInput`] for a block without connection points, a
/// connection point at the block origin (no direction), or a target direction
/// that is not finite and non-zero.
pub fn node_oriented(block: &Block, target_direction: [f64; 3]) -> Result<[[f64; 3]; 3], Error> {
    let source = reference_direction(block)?;
    let target = normalize(target_direction).ok_or_else(|| {
        Error::InvalidInput("target direction must be finite and non-zero".to_owned())
    })?;
    Ok(rotation_from_to(source, target))
}

fn reference_direction(block: &Block) -> Result<[f64; 3], Error> {
    let first = block
        .connection_points
        .first()
        .ok_or_else(|| {
            Error::InvalidInput(format!("block `{}` has no connection points", block.name))
        })
        .and_then(|point| {
            normalize(*point).ok_or_else(|| {
                Error::InvalidInput(format!(
                    "block `{}` has a connection point at the origin",
                    block.name
                ))
            })
        })?;
    let mut sum = [0.0; 3];
    for point in &block.connection_points {
        if let Some(direction) = normalize(*point) {
            sum = add(sum, direction);
        }
    }
    Ok(normalize(sum).unwrap_or(first))
}

pub(crate) fn assign_connections(connections: &[[f64; 3]], slots: &[Slot]) -> Vec<Option<usize>> {
    let mut candidates = Vec::new();
    for (slot_position, slot) in slots.iter().enumerate() {
        let Some(slot_direction) = normalize(slot.direction) else {
            continue;
        };
        for (connection_position, connection) in connections.iter().enumerate() {
            let Some(connection_direction) = normalize(*connection) else {
                continue;
            };
            let similarity = dot(slot_direction, connection_direction);
            if similarity > 0.0 {
                candidates.push((similarity, slot_position, connection_position));
            }
        }
    }
    candidates.sort_by(|left, right| {
        right
            .0
            .total_cmp(&left.0)
            .then(left.1.cmp(&right.1))
            .then(left.2.cmp(&right.2))
    });
    let mut used_connections = vec![false; connections.len()];
    let mut assignment = vec![None; slots.len()];
    for (_, slot_position, connection_position) in candidates {
        if used_connections
            .get(connection_position)
            .copied()
            .unwrap_or(true)
            || assignment.get(slot_position).is_some_and(Option::is_some)
        {
            continue;
        }
        if let Some(slot) = assignment.get_mut(slot_position) {
            *slot = Some(connection_position);
        }
        if let Some(used) = used_connections.get_mut(connection_position) {
            *used = true;
        }
    }
    assignment
}

pub(crate) fn rotation_from_to(from: [f64; 3], to: [f64; 3]) -> [[f64; 3]; 3] {
    let cosine = dot(from, to);
    if cosine > 1.0 - 1e-12 {
        return identity();
    }
    if cosine < -1.0 + 1e-12 {
        return half_turn(from);
    }
    let axis = cross(from, to);
    let sine_squared = dot(axis, axis);
    let coefficient = (1.0 - cosine) / sine_squared;
    let [x, y, z] = axis;
    [
        [
            cosine + coefficient * x * x,
            -z + coefficient * x * y,
            y + coefficient * x * z,
        ],
        [
            z + coefficient * y * x,
            cosine + coefficient * y * y,
            -x + coefficient * y * z,
        ],
        [
            -y + coefficient * z * x,
            x + coefficient * z * y,
            cosine + coefficient * z * z,
        ],
    ]
}

fn half_turn(from: [f64; 3]) -> [[f64; 3]; 3] {
    let [x, y, z] = from;
    let absolute = [x.abs(), y.abs(), z.abs()];
    let first = absolute.first().copied().unwrap_or(0.0);
    let second = absolute.get(1).copied().unwrap_or(0.0);
    let third = absolute.get(2).copied().unwrap_or(0.0);
    let basis = if first <= second && first <= third {
        [1.0, 0.0, 0.0]
    } else if second <= third {
        [0.0, 1.0, 0.0]
    } else {
        [0.0, 0.0, 1.0]
    };
    let axis = normalize(cross(from, basis)).unwrap_or([0.0, 0.0, 1.0]);
    let [ax, ay, az] = axis;
    [
        [2.0 * ax * ax - 1.0, 2.0 * ax * ay, 2.0 * ax * az],
        [2.0 * ay * ax, 2.0 * ay * ay - 1.0, 2.0 * ay * az],
        [2.0 * az * ax, 2.0 * az * ay, 2.0 * az * az - 1.0],
    ]
}

pub(crate) fn matrix_apply(matrix: [[f64; 3]; 3], vector: [f64; 3]) -> [f64; 3] {
    let [x, y, z] = vector;
    let [first, second, third] = matrix;
    [
        row(first, x, y, z),
        row(second, x, y, z),
        row(third, x, y, z),
    ]
}

fn row(row: [f64; 3], x: f64, y: f64, z: f64) -> f64 {
    let [a, b, c] = row;
    a * x + b * y + c * z
}

pub(crate) const fn identity() -> [[f64; 3]; 3] {
    [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]]
}

pub(crate) fn add(first: [f64; 3], second: [f64; 3]) -> [f64; 3] {
    let [x1, y1, z1] = first;
    let [x2, y2, z2] = second;
    [x1 + x2, y1 + y2, z1 + z2]
}

pub(crate) fn sub(first: [f64; 3], second: [f64; 3]) -> [f64; 3] {
    let [x1, y1, z1] = first;
    let [x2, y2, z2] = second;
    [x1 - x2, y1 - y2, z1 - z2]
}

pub(crate) fn scale(vector: [f64; 3], factor: f64) -> [f64; 3] {
    let [x, y, z] = vector;
    [x * factor, y * factor, z * factor]
}

pub(crate) fn dot(first: [f64; 3], second: [f64; 3]) -> f64 {
    let [x1, y1, z1] = first;
    let [x2, y2, z2] = second;
    x1 * x2 + y1 * y2 + z1 * z2
}

fn cross(first: [f64; 3], second: [f64; 3]) -> [f64; 3] {
    let [x1, y1, z1] = first;
    let [x2, y2, z2] = second;
    [y1 * z2 - z1 * y2, z1 * x2 - x1 * z2, x1 * y2 - y1 * x2]
}

pub(crate) fn norm(vector: [f64; 3]) -> f64 {
    dot(vector, vector).sqrt()
}

pub(crate) fn normalize(vector: [f64; 3]) -> Option<[f64; 3]> {
    let length = norm(vector);
    if length < EPSILON || !length.is_finite() {
        None
    } else {
        Some(scale(vector, 1.0 / length))
    }
}
