#![allow(clippy::indexing_slicing, clippy::needless_range_loop)] // fixed-size numeric kernels

//! Convex polyhedra as intersections of half-spaces (Task 9.3.1, 9.4.3).
//!
//! The periodic Voronoi cell and the Wulff construction are both the bounded
//! intersection of a set of half-spaces `n . x <= d`. This module builds that
//! polyhedron by clipping a large cube face by face, and reports its volume,
//! area, centroid and per-face geometry. The algorithm is exact for convex
//! input and deterministic: clipping order is the plane order.

use std::f64::consts::PI;

use serde::{Deserialize, Serialize};

use crate::Error;

/// One half-space `normal . x <= offset`.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct HalfSpace {
    /// Outward normal, not necessarily normalized.
    pub normal: [f64; 3],
    /// Plane offset along the normal.
    pub offset: f64,
}

impl HalfSpace {
    /// Creates a half-space.
    ///
    /// # Errors
    ///
    /// Returns [`Error::InvalidInput`] for a non-finite or zero normal and a
    /// non-finite offset.
    pub fn new(normal: [f64; 3], offset: f64) -> Result<Self, Error> {
        if !normal.iter().all(|value| value.is_finite()) {
            return Err(Error::InvalidInput(
                "half-space normal must be finite".to_owned(),
            ));
        }
        let norm = (normal[0] * normal[0] + normal[1] * normal[1] + normal[2] * normal[2]).sqrt();
        if norm < 1e-12 {
            return Err(Error::InvalidInput(
                "half-space normal must not be zero".to_owned(),
            ));
        }
        if !offset.is_finite() {
            return Err(Error::InvalidInput(
                "half-space offset must be finite".to_owned(),
            ));
        }
        Ok(Self {
            normal: [normal[0] / norm, normal[1] / norm, normal[2] / norm],
            offset: offset / norm,
        })
    }

    /// Signed distance of a point from the boundary, positive inside.
    #[must_use]
    pub fn signed_distance(&self, point: [f64; 3]) -> f64 {
        self.offset - dot(self.normal, point)
    }
}

/// A convex polyhedron with ordered faces.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct Polyhedron {
    /// Vertex coordinates.
    pub vertices: Vec<[f64; 3]>,
    /// Faces as ordered vertex indices with an outward orientation.
    pub faces: Vec<Vec<usize>>,
}

impl Polyhedron {
    /// Intersects the half-spaces, starting from a box large enough to contain
    /// every finite intersection.
    ///
    /// # Errors
    ///
    /// Returns [`Error::InvalidInput`] for an empty plane list or a
    /// degenerate intersection.
    pub fn from_half_spaces(planes: &[HalfSpace]) -> Result<Self, Error> {
        if planes.is_empty() {
            return Err(Error::InvalidInput(
                "a polyhedron needs at least one half-space".to_owned(),
            ));
        }
        let extent = planes
            .iter()
            .map(|plane| plane.offset.abs())
            .fold(1.0_f64, f64::max)
            * 4.0
            + 10.0;
        let mut polyhedron = Self::cube(extent);
        for plane in planes {
            polyhedron = polyhedron.clip(plane, 1e-10 * extent);
            if polyhedron.faces.is_empty() {
                return Err(Error::InvalidInput(
                    "the half-spaces have an empty intersection".to_owned(),
                ));
            }
        }
        Ok(polyhedron)
    }

    /// Axis-aligned cube centred at the origin.
    #[must_use]
    pub fn cube(extent: f64) -> Self {
        let e = extent.abs().max(1e-3);
        let vertices = vec![
            [-e, -e, -e],
            [e, -e, -e],
            [e, e, -e],
            [-e, e, -e],
            [-e, -e, e],
            [e, -e, e],
            [e, e, e],
            [-e, e, e],
        ];
        let faces = vec![
            vec![0, 3, 2, 1],
            vec![4, 5, 6, 7],
            vec![0, 1, 5, 4],
            vec![1, 2, 6, 5],
            vec![2, 3, 7, 6],
            vec![3, 0, 4, 7],
        ];
        Self { vertices, faces }
    }

    /// Clips the polyhedron by one half-space.
    #[must_use]
    pub fn clip(&self, plane: &HalfSpace, epsilon: f64) -> Self {
        let distances: Vec<f64> = self
            .vertices
            .iter()
            .map(|vertex| plane.signed_distance(*vertex))
            .collect();
        let mut vertices = Vec::new();
        let mut faces = Vec::new();
        let mut cap_points: Vec<[f64; 3]> = Vec::new();
        for face in &self.faces {
            let mut polygon: Vec<usize> = Vec::new();
            for (position, &current) in face.iter().enumerate() {
                let Some(&next) = face.get((position + 1) % face.len()) else {
                    continue;
                };
                let (Some(&d_current), Some(&d_next)) =
                    (distances.get(current), distances.get(next))
                else {
                    continue;
                };
                let inside_current = d_current >= -epsilon;
                let inside_next = d_next >= -epsilon;
                if inside_current {
                    push_vertex(&mut vertices, &mut polygon, self.vertices[current]);
                }
                if inside_current != inside_next {
                    let denominator = d_current - d_next;
                    if denominator.abs() > 1e-30 {
                        let t = (d_current / denominator).clamp(0.0, 1.0);
                        let Some(&a) = self.vertices.get(current) else {
                            continue;
                        };
                        let Some(&b) = self.vertices.get(next) else {
                            continue;
                        };
                        let point = [
                            a[0] + (b[0] - a[0]) * t,
                            a[1] + (b[1] - a[1]) * t,
                            a[2] + (b[2] - a[2]) * t,
                        ];
                        push_vertex(&mut vertices, &mut polygon, point);
                        cap_points.push(point);
                    }
                }
            }
            if polygon.len() >= 3 {
                faces.push(polygon);
            }
        }
        if let Some(cap) = cap_face(vertices.as_slice(), &cap_points, plane.normal) {
            faces.push(cap);
        }
        Self { vertices, faces }
    }

    /// Unit outward normal of a face, from its vertex order.
    #[must_use]
    pub fn face_normal(&self, face: usize) -> Option<[f64; 3]> {
        let indices = self.faces.get(face)?;
        let first = *self.vertices.get(*indices.first()?)?;
        let mut normal = [0.0; 3];
        for (position, &index) in indices.iter().enumerate() {
            let next = *indices.get((position + 1) % indices.len())?;
            let a = self.vertices.get(index)?;
            let b = self.vertices.get(next)?;
            normal[0] += (a[1] - b[1]) * (a[2] + b[2]);
            normal[1] += (a[2] - b[2]) * (a[0] + b[0]);
            normal[2] += (a[0] - b[0]) * (a[1] + b[1]);
        }
        let norm = norm(normal);
        if norm < 1e-12 {
            return None;
        }
        let outward = dot(normal, subtract(first, self.centroid())) > 0.0;
        let sign = if outward { 1.0 } else { -1.0 };
        Some([
            normal[0] * sign / norm,
            normal[1] * sign / norm,
            normal[2] * sign / norm,
        ])
    }

    /// Area of one face.
    #[must_use]
    pub fn face_area(&self, face: usize) -> f64 {
        let Some(indices) = self.faces.get(face) else {
            return 0.0;
        };
        let Some(&first) = indices.first() else {
            return 0.0;
        };
        let Some(&origin) = self.vertices.get(first) else {
            return 0.0;
        };
        let mut area = 0.0;
        for position in 1..indices.len().saturating_sub(1) {
            let (Some(&i), Some(&j)) = (indices.get(position), indices.get(position + 1)) else {
                continue;
            };
            let (Some(&a), Some(&b)) = (self.vertices.get(i), self.vertices.get(j)) else {
                continue;
            };
            area += 0.5 * norm(cross(subtract(a, origin), subtract(b, origin)));
        }
        area
    }

    /// Total surface area.
    #[must_use]
    pub fn surface_area(&self) -> f64 {
        (0..self.faces.len()).map(|face| self.face_area(face)).sum()
    }

    /// Volume by the signed tetrahedra of each face around the origin.
    #[must_use]
    pub fn volume(&self) -> f64 {
        let mut volume = 0.0;
        for indices in &self.faces {
            let Some(&first) = indices.first() else {
                continue;
            };
            let Some(&origin) = self.vertices.get(first) else {
                continue;
            };
            for position in 1..indices.len().saturating_sub(1) {
                let (Some(&i), Some(&j)) = (indices.get(position), indices.get(position + 1))
                else {
                    continue;
                };
                let (Some(&a), Some(&b)) = (self.vertices.get(i), self.vertices.get(j)) else {
                    continue;
                };
                volume += dot(origin, cross(a, b)) / 6.0;
            }
        }
        volume.abs()
    }

    /// Centroid of the vertices; equal to the volume centroid for a simplex,
    /// and a stable centre for clipping.
    #[must_use]
    pub fn centroid(&self) -> [f64; 3] {
        let mut total = [0.0; 3];
        for vertex in &self.vertices {
            total = add(total, *vertex);
        }
        let count = self.vertices.len().max(1) as f64;
        [total[0] / count, total[1] / count, total[2] / count]
    }

    /// Face centroid.
    #[must_use]
    pub fn face_centroid(&self, face: usize) -> [f64; 3] {
        let Some(indices) = self.faces.get(face) else {
            return [0.0; 3];
        };
        let mut total = [0.0; 3];
        let mut count = 0.0;
        for &index in indices {
            if let Some(&vertex) = self.vertices.get(index) {
                total = add(total, vertex);
                count += 1.0;
            }
        }
        if count == 0.0 {
            return [0.0; 3];
        }
        [total[0] / count, total[1] / count, total[2] / count]
    }

    /// Solid angle of one face as seen from `origin`, in steradians.
    #[must_use]
    pub fn face_solid_angle(&self, face: usize, origin: [f64; 3]) -> f64 {
        let area = self.face_area(face);
        if area <= 0.0 {
            return 0.0;
        }
        let centre = subtract(self.face_centroid(face), origin);
        let distance = norm(centre);
        if distance < 1e-12 {
            return 0.0;
        }
        let Some(normal) = self.face_normal(face) else {
            return 0.0;
        };
        let cosine = dot(normal, centre) / distance;
        if cosine <= 0.0 {
            0.0
        } else {
            area * cosine / (distance * distance)
        }
    }

    /// Triangles of the surface with outward normals.
    #[must_use]
    pub fn triangles(&self) -> Vec<([f64; 3], [[f64; 3]; 3])> {
        let mut triangles = Vec::new();
        for (face, indices) in self.faces.iter().enumerate() {
            let Some(normal) = self.face_normal(face) else {
                continue;
            };
            let Some(&first) = indices.first() else {
                continue;
            };
            let Some(&origin) = self.vertices.get(first) else {
                continue;
            };
            for position in 1..indices.len().saturating_sub(1) {
                let (Some(&i), Some(&j)) = (indices.get(position), indices.get(position + 1))
                else {
                    continue;
                };
                let (Some(&a), Some(&b)) = (self.vertices.get(i), self.vertices.get(j)) else {
                    continue;
                };
                triangles.push((normal, [origin, a, b]));
            }
        }
        triangles
    }
}

/// Spherical area of a full sphere in steradians.
pub const FULL_SPHERE: f64 = 4.0 * PI;

fn cap_face(vertices: &[[f64; 3]], points: &[[f64; 3]], normal: [f64; 3]) -> Option<Vec<usize>> {
    if points.len() < 3 {
        return None;
    }
    let mut total = [0.0; 3];
    for point in points {
        total = add(total, *point);
    }
    let count = points.len() as f64;
    let centre = [total[0] / count, total[1] / count, total[2] / count];
    let (u, v) = plane_basis(normal);
    let mut ordered: Vec<(f64, [f64; 3])> = points
        .iter()
        .map(|point| {
            let delta = subtract(*point, centre);
            (dot(delta, v).atan2(dot(delta, u)), *point)
        })
        .collect();
    ordered.sort_by(|first, second| {
        first
            .0
            .partial_cmp(&second.0)
            .unwrap_or(std::cmp::Ordering::Equal)
    });
    let mut face = Vec::with_capacity(ordered.len());
    for (_, point) in ordered {
        let mut index = None;
        for (candidate, vertex) in vertices.iter().enumerate() {
            if norm(subtract(*vertex, point)) < 1e-8 {
                index = Some(candidate);
                break;
            }
        }
        if let Some(index) = index {
            face.push(index);
        }
    }
    if face.len() < 3 {
        return None;
    }
    let area_normal = polygon_normal(vertices, &face);
    if dot(area_normal, normal) < 0.0 {
        face.reverse();
    }
    Some(face)
}

fn polygon_normal(vertices: &[[f64; 3]], face: &[usize]) -> [f64; 3] {
    let mut normal = [0.0; 3];
    for (position, &index) in face.iter().enumerate() {
        let next = face[(position + 1) % face.len()];
        let (Some(a), Some(b)) = (vertices.get(index), vertices.get(next)) else {
            continue;
        };
        normal[0] += (a[1] - b[1]) * (a[2] + b[2]);
        normal[1] += (a[2] - b[2]) * (a[0] + b[0]);
        normal[2] += (a[0] - b[0]) * (a[1] + b[1]);
    }
    normal
}

fn plane_basis(normal: [f64; 3]) -> ([f64; 3], [f64; 3]) {
    let reference = if normal[0].abs() < 0.9 {
        [1.0, 0.0, 0.0]
    } else {
        [0.0, 1.0, 0.0]
    };
    let u = normalize(cross(normal, reference));
    let v = normalize(cross(normal, u));
    (u, v)
}

fn push_vertex(vertices: &mut Vec<[f64; 3]>, polygon: &mut Vec<usize>, point: [f64; 3]) {
    for (index, vertex) in vertices.iter().enumerate() {
        if norm(subtract(*vertex, point)) < 1e-9 {
            polygon.push(index);
            return;
        }
    }
    let index = vertices.len();
    vertices.push(point);
    polygon.push(index);
}

fn add(first: [f64; 3], second: [f64; 3]) -> [f64; 3] {
    [
        first[0] + second[0],
        first[1] + second[1],
        first[2] + second[2],
    ]
}

fn subtract(first: [f64; 3], second: [f64; 3]) -> [f64; 3] {
    [
        first[0] - second[0],
        first[1] - second[1],
        first[2] - second[2],
    ]
}

fn cross(first: [f64; 3], second: [f64; 3]) -> [f64; 3] {
    [
        first[1] * second[2] - first[2] * second[1],
        first[2] * second[0] - first[0] * second[2],
        first[0] * second[1] - first[1] * second[0],
    ]
}

fn dot(first: [f64; 3], second: [f64; 3]) -> f64 {
    first[0] * second[0] + first[1] * second[1] + first[2] * second[2]
}

fn norm(vector: [f64; 3]) -> f64 {
    dot(vector, vector).sqrt()
}

fn normalize(vector: [f64; 3]) -> [f64; 3] {
    let length = norm(vector);
    if length < 1e-12 {
        [0.0; 3]
    } else {
        [vector[0] / length, vector[1] / length, vector[2] / length]
    }
}
