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
                let tetrahedra = [
                    [0, 1, 2, 6],
                    [0, 2, 3, 6],
                    [0, 3, 7, 6],
                    [0, 7, 4, 6],
                    [0, 4, 5, 6],
                    [0, 5, 1, 6],
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
    let mut inside = Vec::new();
    let mut outside = Vec::new();
    for corner in corners {
        let value = grid
            .get(corner.0, corner.1, corner.2)
            .ok_or_else(|| Error::InvalidInput("grid corner is out of range".to_owned()))?;
        if value >= isovalue {
            inside.push((corner, value));
        } else {
            outside.push((corner, value));
        }
    }
    if inside.is_empty() || outside.is_empty() {
        return Ok(());
    }
    let mut points = Vec::new();
    for (a, va) in &inside {
        for (b, vb) in &outside {
            let denominator = va - vb;
            if denominator.abs() < 1e-15 {
                continue;
            }
            let t = (isovalue - vb) / denominator;
            let pa = grid.position(a.0, a.1, a.2);
            let pb = grid.position(b.0, b.1, b.2);
            points.push([
                pa[0] + t * (pb[0] - pa[0]),
                pa[1] + t * (pb[1] - pa[1]),
                pa[2] + t * (pb[2] - pa[2]),
            ]);
        }
    }
    if points.len() == 3 {
        add_triangle(mesh, points[0], points[1], points[2]);
    } else if points.len() == 4 {
        add_triangle(mesh, points[0], points[1], points[2]);
        add_triangle(mesh, points[0], points[2], points[3]);
    }
    Ok(())
}

fn add_triangle(mesh: &mut TriangleMesh, a: [f64; 3], b: [f64; 3], c: [f64; 3]) {
    let base = mesh.vertices.len();
    mesh.vertices.extend_from_slice(&[a, b, c]);
    mesh.triangles.push([base, base + 1, base + 2]);
}
