//! Periodic marching tetrahedra for scalar grids (Task 5.4).

#![allow(clippy::indexing_slicing)] // fixed grid indexing, bounds checked by loops

use serde::{Deserialize, Serialize};
use xtalforge_core::Error;
use xtalforge_io::grid::Grid3;

/// A triangle mesh in Cartesian coordinates.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct TriangleMesh {
    /// Vertices in Angstrom.
    pub vertices: Vec<[f64; 3]>,
    /// Triangles as vertex index triples.
    pub triangles: Vec<[usize; 3]>,
}

impl TriangleMesh {
    /// Number of triangles.
    #[must_use]
    pub fn len(&self) -> usize {
        self.triangles.len()
    }

    /// Whether the mesh is empty.
    #[must_use]
    pub fn is_empty(&self) -> bool {
        self.triangles.is_empty()
    }
}

/// Extracts an isosurface with marching tetrahedra.
///
/// The grid is split into six tetrahedra per cell, which avoids the ambiguous
/// cases of marching cubes. Values equal to the isovalue produce no triangle.
///
/// # Errors
///
/// Returns [`Error::InvalidInput`] for a non-finite isovalue or a grid with
/// fewer than two samples along an axis.
pub fn isosurface(grid: &Grid3, isovalue: f64) -> Result<TriangleMesh, Error> {
    if !isovalue.is_finite() {
        return Err(Error::InvalidInput("isovalue must be finite".to_owned()));
    }
    if grid.shape.iter().any(|value| *value < 2) {
        return Err(Error::InvalidInput(
            "an isosurface needs at least two samples per axis".to_owned(),
        ));
    }
    let mut mesh = TriangleMesh {
        vertices: Vec::new(),
        triangles: Vec::new(),
    };
    for z in 0..grid.shape[2] - 1 {
        for y in 0..grid.shape[1] - 1 {
            for x in 0..grid.shape[0] - 1 {
                let corners = [
                    (x, y, z),
                    (x + 1, y, z),
                    (x + 1, y + 1, z),
                    (x, y + 1, z),
                    (x, y, z + 1),
                    (x + 1, y, z + 1),
                    (x + 1, y + 1, z + 1),
                    (x, y + 1, z + 1),
                ];
                // Kuhn decomposition about the 0-7 diagonal, in the corner
                // order of the reference marching-tetrahedra table.
                let tetrahedra = [
                    [0, 4, 7, 6],
                    [0, 4, 5, 6],
                    [0, 3, 7, 6],
                    [0, 3, 2, 6],
                    [0, 1, 5, 6],
                    [0, 1, 2, 6],
                ];
                for tetrahedron in tetrahedra {
                    let points = tetrahedron.map(|index| corners[index]);
                    tetrahedron_triangles(grid, points, isovalue, &mut mesh)?;
                }
            }
        }
    }
    Ok(mesh)
}

fn tetrahedron_triangles(
    grid: &Grid3,
    corners: [(usize, usize, usize); 4],
    isovalue: f64,
    mesh: &mut TriangleMesh,
) -> Result<(), Error> {
    let mut points = [[0.0; 3]; 4];
    let mut values = [0.0; 4];
    for index in 0..4 {
        let corner = corners[index];
        values[index] = grid
            .get(corner.0, corner.1, corner.2)
            .ok_or_else(|| Error::InvalidInput("grid corner is out of range".to_owned()))?;
        points[index] = grid.position(corner.0, corner.1, corner.2);
    }
    let mut pattern = 0u8;
    for (index, value) in values.iter().enumerate() {
        if *value >= isovalue {
            pattern |= 1 << index;
        }
    }
    let triangles: &[[usize; 3]] = match pattern {
        0b0001 | 0b1110 => &[[0, 1, 2]],
        0b0010 | 0b1101 => &[[0, 3, 4]],
        0b0100 | 0b1011 => &[[1, 3, 5]],
        0b1000 | 0b0111 => &[[2, 4, 5]],
        0b0011 | 0b1100 => &[[1, 2, 4], [1, 4, 3]],
        0b0101 | 0b1010 => &[[0, 2, 5], [0, 5, 3]],
        0b1001 | 0b0110 => &[[0, 1, 5], [0, 5, 4]],
        _ => &[],
    };
    for triangle in triangles {
        let mut crossing = [[0.0; 3]; 3];
        for (slot, edge) in triangle.iter().enumerate() {
            let (u, v) = tet_edge(*edge);
            crossing[slot] = crossing_point(points[u], points[v], values[u], values[v], isovalue);
        }
        add_triangle(mesh, crossing[0], crossing[1], crossing[2]);
    }
    Ok(())
}

/// The six edges of a tetrahedron, in the reference table's order.
const TET_EDGES: [(usize, usize); 6] = [(0, 1), (0, 2), (0, 3), (1, 2), (1, 3), (2, 3)];

fn tet_edge(index: usize) -> (usize, usize) {
    TET_EDGES.get(index).copied().unwrap_or((0, 1))
}

fn crossing_point(
    start: [f64; 3],
    end: [f64; 3],
    start_value: f64,
    end_value: f64,
    isovalue: f64,
) -> [f64; 3] {
    let span = end_value - start_value;
    let t = if span.abs() < 1e-12 {
        0.5
    } else {
        (isovalue - start_value) / span
    };
    [
        start[0] + t * (end[0] - start[0]),
        start[1] + t * (end[1] - start[1]),
        start[2] + t * (end[2] - start[2]),
    ]
}

fn add_triangle(mesh: &mut TriangleMesh, a: [f64; 3], b: [f64; 3], c: [f64; 3]) {
    let base = mesh.vertices.len();
    mesh.vertices.extend_from_slice(&[a, b, c]);
    mesh.triangles.push([base, base + 1, base + 2]);
}
