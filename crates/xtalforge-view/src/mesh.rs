//! Printable triangle mesh and ASCII STL (Task 3.15.1).

use serde::{Deserialize, Serialize};

use crate::scene::SceneData;

/// One triangle with its outward normal.
#[derive(Debug, Clone, Copy, PartialEq, Serialize, Deserialize)]
pub struct Triangle {
    /// Outward unit normal.
    pub normal: [f64; 3],
    /// Three vertices.
    pub vertices: [[f64; 3]; 3],
}

/// Builds a printable mesh: each atom is an octahedron and each bond a
/// rectangular prism. The result is watertight enough for a slicer when the
/// radii overlap.
#[must_use]
pub fn build_mesh(scene: &SceneData) -> Vec<Triangle> {
    let mut triangles = Vec::new();
    for atom in &scene.atoms {
        if atom.hidden || atom.radius <= 0.0 {
            continue;
        }
        let r = atom.radius;
        let [x, y, z] = atom.position;
        let vertices = [
            [x + r, y, z],
            [x - r, y, z],
            [x, y + r, z],
            [x, y - r, z],
            [x, y, z + r],
            [x, y, z - r],
        ];
        for (a, b, c) in [
            (0, 2, 4),
            (2, 1, 4),
            (1, 3, 4),
            (3, 0, 4),
            (2, 0, 5),
            (1, 2, 5),
            (3, 1, 5),
            (0, 3, 5),
        ] {
            if let (Some(first), Some(second), Some(third)) =
                (vertices.get(a), vertices.get(b), vertices.get(c))
            {
                push_triangle(&mut triangles, *first, *second, *third);
            }
        }
    }
    for bond in &scene.bonds {
        let axis = subtract(bond.end, bond.start);
        let length = norm(axis);
        if length < 1e-9 {
            continue;
        }
        let direction = scale(axis, 1.0 / length);
        let reference = if direction[0].abs() < 0.9 {
            [1.0, 0.0, 0.0]
        } else {
            [0.0, 1.0, 0.0]
        };
        let side = normalize(cross(direction, reference));
        let up = cross(direction, side);
        let r = bond.radius;
        let start = bond.start;
        let end = bond.end;
        let corners = [
            add(start, add(scale(side, r), scale(up, r))),
            add(start, add(scale(side, -r), scale(up, r))),
            add(start, add(scale(side, -r), scale(up, -r))),
            add(start, add(scale(side, r), scale(up, -r))),
            add(end, add(scale(side, r), scale(up, r))),
            add(end, add(scale(side, -r), scale(up, r))),
            add(end, add(scale(side, -r), scale(up, -r))),
            add(end, add(scale(side, r), scale(up, -r))),
        ];
        for (a, b, c) in [
            (0, 1, 2),
            (0, 2, 3),
            (4, 6, 5),
            (4, 7, 6),
            (0, 4, 5),
            (0, 5, 1),
            (1, 5, 6),
            (1, 6, 2),
            (2, 6, 7),
            (2, 7, 3),
            (3, 7, 4),
            (3, 4, 0),
        ] {
            if let (Some(first), Some(second), Some(third)) =
                (corners.get(a), corners.get(b), corners.get(c))
            {
                push_triangle(&mut triangles, *first, *second, *third);
            }
        }
    }
    triangles
}

fn push_triangle(
    triangles: &mut Vec<Triangle>,
    first: [f64; 3],
    second: [f64; 3],
    third: [f64; 3],
) {
    let normal = normalize(cross(subtract(second, first), subtract(third, first)));
    triangles.push(Triangle {
        normal,
        vertices: [first, second, third],
    });
}

/// Writes an ASCII STL.
#[must_use]
pub fn write_stl(triangles: &[Triangle], name: &str) -> String {
    let mut out = format!("solid {name}\n");
    for triangle in triangles {
        out.push_str(&format!(
            "  facet normal {:.6} {:.6} {:.6}\n    outer loop\n",
            triangle.normal[0], triangle.normal[1], triangle.normal[2]
        ));
        for vertex in triangle.vertices {
            out.push_str(&format!(
                "      vertex {:.6} {:.6} {:.6}\n",
                vertex[0], vertex[1], vertex[2]
            ));
        }
        out.push_str("    endloop\n  endfacet\n");
    }
    out.push_str(&format!("endsolid {name}\n"));
    out
}

pub(crate) fn subtract(first: [f64; 3], second: [f64; 3]) -> [f64; 3] {
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

pub(crate) fn cross(first: [f64; 3], second: [f64; 3]) -> [f64; 3] {
    [
        first[1] * second[2] - first[2] * second[1],
        first[2] * second[0] - first[0] * second[2],
        first[0] * second[1] - first[1] * second[0],
    ]
}

pub(crate) fn norm(vector: [f64; 3]) -> f64 {
    (vector[0] * vector[0] + vector[1] * vector[1] + vector[2] * vector[2]).sqrt()
}

pub(crate) fn normalize(vector: [f64; 3]) -> [f64; 3] {
    let length = norm(vector);
    if length < 1e-12 {
        [0.0, 0.0, 0.0]
    } else {
        scale(vector, 1.0 / length)
    }
}
