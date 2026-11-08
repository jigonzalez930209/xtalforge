//! Wigner-Seitz Brillouin zone (Task 2.10.1).
//!
//! The zone is the intersection of the half-spaces
//! `G . x <= |G|^2 / 2` over the reciprocal lattice vectors `G`, computed by
//! clipping a large cube. Reciprocal vectors include `2 pi`, so distances are
//! in `1/Angstrom` and match the band-structure x axis (Task 2.10.3).

use serde::{Deserialize, Serialize};
use xtalforge_core::{Error, Lattice};

/// A convex polyhedron: vertices, faces and edges.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct BrillouinZone {
    /// Cartesian vertices in `1/Angstrom`.
    pub vertices: Vec<[f64; 3]>,
    /// Face vertex loops, counter-clockwise seen from outside.
    pub faces: Vec<Vec<usize>>,
    /// Unique edges as vertex-index pairs.
    pub edges: Vec<(usize, usize)>,
}

impl BrillouinZone {
    /// Number of faces.
    #[must_use]
    pub fn n_faces(&self) -> usize {
        self.faces.len()
    }

    /// Number of vertices.
    #[must_use]
    pub fn n_vertices(&self) -> usize {
        self.vertices.len()
    }
}

#[derive(Debug, Clone)]
struct Face {
    points: Vec<[f64; 3]>,
    normal: [f64; 3],
}

/// Builds the Wigner-Seitz cell of the reciprocal lattice.
///
/// # Errors
///
/// Returns [`Error::InvalidInput`] when the reciprocal lattice is degenerate or
/// the clipping produces no volume.
pub fn brillouin_zone(lattice: &Lattice) -> Result<BrillouinZone, Error> {
    let reciprocal = lattice.reciprocal(true)?;
    let mut vectors = Vec::new();
    for h in -3..=3 {
        for k in -3..=3 {
            for l in -3..=3 {
                if h == 0 && k == 0 && l == 0 {
                    continue;
                }
                let vector = reciprocal.to_cart([f64::from(h), f64::from(k), f64::from(l)]);
                vectors.push(vector);
            }
        }
    }
    vectors.sort_by(|first, second| {
        norm(*first)
            .partial_cmp(&norm(*second))
            .unwrap_or(std::cmp::Ordering::Equal)
    });
    let shortest = vectors.first().map_or(0.0, |vector| norm(*vector));
    if shortest <= 0.0 {
        return Err(Error::InvalidInput(
            "the reciprocal lattice has no nonzero vector".to_owned(),
        ));
    }
    let radius = 2.0 * shortest;
    let mut faces = cube(radius);
    for vector in vectors {
        let length = norm(vector);
        if length > 2.0 * radius {
            continue;
        }
        faces = clip(faces, vector, length * length / 2.0);
        if faces.is_empty() {
            return Err(Error::InvalidInput(
                "the Brillouin zone clipping produced no volume".to_owned(),
            ));
        }
    }
    assemble(faces)
}

fn cube(radius: f64) -> Vec<Face> {
    let corners = [
        [-radius, -radius, -radius],
        [radius, -radius, -radius],
        [radius, radius, -radius],
        [-radius, radius, -radius],
        [-radius, -radius, radius],
        [radius, -radius, radius],
        [radius, radius, radius],
        [-radius, radius, radius],
    ];
    let faces = [
        (vec![0, 3, 2, 1], [0.0, 0.0, -1.0]),
        (vec![4, 5, 6, 7], [0.0, 0.0, 1.0]),
        (vec![0, 1, 5, 4], [0.0, -1.0, 0.0]),
        (vec![3, 7, 6, 2], [0.0, 1.0, 0.0]),
        (vec![0, 4, 7, 3], [-1.0, 0.0, 0.0]),
        (vec![1, 2, 6, 5], [1.0, 0.0, 0.0]),
    ];
    faces
        .into_iter()
        .map(|(indices, normal)| Face {
            points: indices
                .iter()
                .filter_map(|index| corners.get(*index).copied())
                .collect(),
            normal,
        })
        .collect()
}

fn clip(faces: Vec<Face>, normal: [f64; 3], offset: f64) -> Vec<Face> {
    let mut clipped = Vec::new();
    let mut cap: Vec<[f64; 3]> = Vec::new();
    for face in faces {
        let mut points = Vec::new();
        for (first, second) in face.points.iter().zip(face.points.iter().cycle().skip(1)) {
            let first = *first;
            let second = *second;
            let first_distance = dot(first, normal) - offset;
            let second_distance = dot(second, normal) - offset;
            let first_inside = first_distance <= 1e-10;
            let second_inside = second_distance <= 1e-10;
            if first_inside {
                points.push(first);
            }
            if first_inside != second_inside {
                let denominator = first_distance - second_distance;
                if denominator.abs() > 1e-15 {
                    let t = first_distance / denominator;
                    let point = [
                        first[0] + t * (second[0] - first[0]),
                        first[1] + t * (second[1] - first[1]),
                        first[2] + t * (second[2] - first[2]),
                    ];
                    points.push(point);
                    cap.push(point);
                }
            }
        }
        if points.len() >= 3 {
            clipped.push(Face {
                points: dedupe_loop(points),
                normal: face.normal,
            });
        }
    }
    if cap.len() >= 3
        && let Some(points) = order_cap(cap, normal)
    {
        clipped.push(Face { points, normal });
    }
    clipped
}

fn order_cap(points: Vec<[f64; 3]>, normal: [f64; 3]) -> Option<Vec<[f64; 3]>> {
    let mut unique: Vec<[f64; 3]> = Vec::new();
    for point in points {
        if !unique
            .iter()
            .any(|existing| norm(subtract(*existing, point)) < 1e-8)
        {
            unique.push(point);
        }
    }
    if unique.len() < 3 {
        return None;
    }
    let mut centroid = [0.0; 3];
    for point in &unique {
        centroid = add(centroid, *point);
    }
    let count = unique.len() as f64;
    centroid = scale(centroid, 1.0 / count);
    let reference = if normal[0].abs() < 0.9 {
        [1.0, 0.0, 0.0]
    } else {
        [0.0, 1.0, 0.0]
    };
    let u = normalize(cross(normal, reference))?;
    let v = cross(normal, u);
    unique.sort_by(|first, second| {
        let first_delta = subtract(*first, centroid);
        let second_delta = subtract(*second, centroid);
        let first_angle = dot(first_delta, v).atan2(dot(first_delta, u));
        let second_angle = dot(second_delta, v).atan2(dot(second_delta, u));
        first_angle
            .partial_cmp(&second_angle)
            .unwrap_or(std::cmp::Ordering::Equal)
    });
    Some(unique)
}

fn dedupe_loop(points: Vec<[f64; 3]>) -> Vec<[f64; 3]> {
    let mut unique: Vec<[f64; 3]> = Vec::new();
    for point in points {
        if !unique
            .iter()
            .any(|existing| norm(subtract(*existing, point)) < 1e-8)
        {
            unique.push(point);
        }
    }
    unique
}

fn assemble(faces: Vec<Face>) -> Result<BrillouinZone, Error> {
    let mut vertices: Vec<[f64; 3]> = Vec::new();
    let mut indexed_faces = Vec::new();
    for face in faces {
        let mut indices = Vec::new();
        for point in face.points {
            let index = vertices
                .iter()
                .position(|existing| norm(subtract(*existing, point)) < 1e-8)
                .unwrap_or_else(|| {
                    vertices.push(point);
                    vertices.len() - 1
                });
            indices.push(index);
        }
        if indices.len() >= 3 {
            indexed_faces.push(indices);
        }
    }
    if vertices.len() < 4 || indexed_faces.len() < 4 {
        return Err(Error::InvalidInput(
            "the Brillouin zone has fewer than four faces".to_owned(),
        ));
    }
    let mut edges: Vec<(usize, usize)> = Vec::new();
    for face in &indexed_faces {
        for (first, second) in face.iter().zip(face.iter().cycle().skip(1)) {
            let first = *first;
            let second = *second;
            let edge = if first < second {
                (first, second)
            } else {
                (second, first)
            };
            if !edges.contains(&edge) {
                edges.push(edge);
            }
        }
    }
    Ok(BrillouinZone {
        vertices,
        faces: indexed_faces,
        edges,
    })
}

fn norm(vector: [f64; 3]) -> f64 {
    dot(vector, vector).sqrt()
}

fn subtract(first: [f64; 3], second: [f64; 3]) -> [f64; 3] {
    [
        first[0] - second[0],
        first[1] - second[1],
        first[2] - second[2],
    ]
}

fn add(first: [f64; 3], second: [f64; 3]) -> [f64; 3] {
    [
        first[0] + second[0],
        first[1] + second[1],
        first[2] + second[2],
    ]
}

fn scale(vector: [f64; 3], factor: f64) -> [f64; 3] {
    [vector[0] * factor, vector[1] * factor, vector[2] * factor]
}

fn dot(first: [f64; 3], second: [f64; 3]) -> f64 {
    first[0] * second[0] + first[1] * second[1] + first[2] * second[2]
}

fn cross(first: [f64; 3], second: [f64; 3]) -> [f64; 3] {
    [
        first[1] * second[2] - first[2] * second[1],
        first[2] * second[0] - first[0] * second[2],
        first[0] * second[1] - first[1] * second[0],
    ]
}

fn normalize(vector: [f64; 3]) -> Option<[f64; 3]> {
    let length = norm(vector);
    if length < 1e-12 {
        None
    } else {
        Some(scale(vector, 1.0 / length))
    }
}
