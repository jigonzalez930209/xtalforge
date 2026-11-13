//! Indexed mesh and vertex normals (Task 14.4.4).
//!
//! [`crate::mesh::build_mesh`] emits independent triangles with one flat
//! normal each, as STL expects. OBJ and glTF share vertices instead: a vertex
//! at the same position in several triangles gets one averaged normal. Welding
//! uses the exact bit pattern of the coordinates, so two positions that differ
//! in the last bit stay separate.

use std::collections::BTreeMap;

use serde::{Deserialize, Serialize};

use crate::mesh::{Triangle, cross, normalize, subtract};

/// One welded mesh vertex with its averaged unit normal.
#[derive(Debug, Clone, Copy, PartialEq, Serialize, Deserialize)]
pub struct MeshVertex {
    /// Cartesian position in Angstrom.
    pub position: [f64; 3],
    /// Unit normal averaged over the incident triangles.
    ///
    /// It is the zero vector when the incident triangles are degenerate.
    pub normal: [f64; 3],
}

/// Triangle mesh with shared vertices, ready for OBJ and glTF.
#[derive(Debug, Clone, Default, PartialEq, Serialize, Deserialize)]
pub struct IndexedMesh {
    /// Unique vertices in first-seen order.
    pub vertices: Vec<MeshVertex>,
    /// Triangle corner indices into `vertices`.
    pub indices: Vec<[usize; 3]>,
}

impl IndexedMesh {
    /// Number of triangles.
    #[must_use]
    pub fn triangle_count(&self) -> usize {
        self.indices.len()
    }

    /// Whether the mesh has no triangles.
    #[must_use]
    pub fn is_empty(&self) -> bool {
        self.indices.is_empty()
    }
}

/// Welds the triangle corners by position and averages area-weighted normals.
///
/// The result is deterministic: vertices keep the order in which they first
/// appear and normals are summed in triangle order. Each normal is the
/// normalized sum of the incident triangle cross products, which weights a
/// triangle by twice its area; a vertex with only degenerate incident
/// triangles gets the zero vector.
#[must_use]
pub fn build_indexed_mesh(triangles: &[Triangle]) -> IndexedMesh {
    let mut vertices: Vec<MeshVertex> = Vec::new();
    let mut index_of: BTreeMap<[u64; 3], usize> = BTreeMap::new();
    let mut indices: Vec<[usize; 3]> = Vec::with_capacity(triangles.len());
    let mut accumulated: Vec<[f64; 3]> = Vec::new();
    for triangle in triangles {
        let [first, second, third] = triangle.vertices;
        let weighted = cross(subtract(second, first), subtract(third, first));
        let mut corners = [0usize; 3];
        for (slot, position) in corners.iter_mut().zip(triangle.vertices.iter()) {
            let key = position_key(*position);
            let index = match index_of.get(&key) {
                Some(index) => *index,
                None => {
                    let index = vertices.len();
                    index_of.insert(key, index);
                    vertices.push(MeshVertex {
                        position: *position,
                        normal: [0.0; 3],
                    });
                    accumulated.push([0.0; 3]);
                    index
                }
            };
            *slot = index;
        }
        for index in corners {
            if let Some(sum) = accumulated.get_mut(index) {
                sum[0] += weighted[0];
                sum[1] += weighted[1];
                sum[2] += weighted[2];
            }
        }
        indices.push(corners);
    }
    for (vertex, sum) in vertices.iter_mut().zip(accumulated.iter()) {
        vertex.normal = normalize(*sum);
    }
    IndexedMesh { vertices, indices }
}

/// Returns one normal per triangle corner, averaged over shared positions.
///
/// The layout matches the input: entry `t[c]` is the normal of corner `c` of
/// triangle `t`. Corners that share an exact position share the same averaged
/// normal, which makes the result suitable for smooth shading. A corner whose
/// incident triangles are all degenerate keeps the zero vector.
#[must_use]
pub fn vertex_normals(triangles: &[Triangle]) -> Vec<[[f64; 3]; 3]> {
    let mesh = build_indexed_mesh(triangles);
    mesh.indices
        .iter()
        .map(|corners| {
            let mut normals = [[0.0; 3]; 3];
            for (slot, index) in normals.iter_mut().zip(corners.iter()) {
                if let Some(vertex) = mesh.vertices.get(*index) {
                    *slot = vertex.normal;
                }
            }
            normals
        })
        .collect()
}

/// Exact coordinate key, with `-0.0` folded into `0.0`.
fn position_key(position: [f64; 3]) -> [u64; 3] {
    position.map(|value| {
        if value == 0.0 {
            0.0f64.to_bits()
        } else {
            value.to_bits()
        }
    })
}
