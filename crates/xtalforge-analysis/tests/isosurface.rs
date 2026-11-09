//! Integration tests for marching tetrahedra (Task 5.4).

#![allow(
    clippy::unwrap_used,
    clippy::expect_used,
    clippy::panic,
    clippy::indexing_slicing
)]

use xtalforge_analysis::isosurface;
use xtalforge_io::grid::Grid3;

#[test]
fn sphere_isosurface_has_vertices_and_triangles() {
    let n = 20;
    let step = 0.1;
    let mut data = Vec::with_capacity(n * n * n);
    for z in 0..n {
        for y in 0..n {
            for x in 0..n {
                let position = [x as f64 * step, y as f64 * step, z as f64 * step];
                let radius = ((position[0] - 1.0).powi(2)
                    + (position[1] - 1.0).powi(2)
                    + (position[2] - 1.0).powi(2))
                .sqrt();
                data.push(radius);
            }
        }
    }
    let grid = Grid3::new(
        [0.0; 3],
        [[step, 0.0, 0.0], [0.0, step, 0.0], [0.0, 0.0, step]],
        [n, n, n],
        data,
        false,
    )
    .unwrap();
    let mesh = isosurface(&grid, 0.47).unwrap();
    assert!(!mesh.is_empty());
    assert_eq!(mesh.triangles.len(), mesh.vertices.len() / 3);
    assert!(
        mesh.vertices
            .iter()
            .all(|vertex| vertex.iter().all(|v| v.is_finite()))
    );
}

#[test]
fn invalid_grids_and_values_are_rejected() {
    let grid = Grid3::new(
        [0.0; 3],
        [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]],
        [2, 2, 2],
        vec![0.0; 8],
        false,
    )
    .unwrap();
    assert!(isosurface(&grid, f64::NAN).is_err());
    let flat = Grid3::new(
        [0.0; 3],
        [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]],
        [1, 2, 2],
        vec![0.0; 4],
        false,
    )
    .unwrap();
    assert!(isosurface(&flat, 0.5).is_err());
}

#[test]
fn marching_tetrahedra_produces_a_watertight_sphere() {
    use std::collections::BTreeMap;

    let n = 24;
    let step = 0.1;
    let mut data = Vec::with_capacity(n * n * n);
    for z in 0..n {
        for y in 0..n {
            for x in 0..n {
                let position = [x as f64 * step, y as f64 * step, z as f64 * step];
                data.push(
                    ((position[0] - 1.05).powi(2)
                        + (position[1] - 1.05).powi(2)
                        + (position[2] - 1.05).powi(2))
                    .sqrt(),
                );
            }
        }
    }
    let grid = Grid3::new(
        [0.0; 3],
        [[step, 0.0, 0.0], [0.0, step, 0.0], [0.0, 0.0, step]],
        [n, n, n],
        data,
        false,
    )
    .unwrap();
    let mesh = isosurface(&grid, 0.47).unwrap();
    let mut map = BTreeMap::new();
    let mut vertices = Vec::new();
    for point in &mesh.vertices {
        let key = [point[0].to_bits(), point[1].to_bits(), point[2].to_bits()];
        map.entry(key).or_insert_with(|| {
            vertices.push(*point);
            vertices.len() - 1
        });
    }
    let mut remap = BTreeMap::new();
    for (position, key) in map.keys().enumerate() {
        remap.insert(*key, position);
    }
    let _ = &vertices;
    let mut edges: BTreeMap<(usize, usize), usize> = BTreeMap::new();
    for triangle in &mesh.triangles {
        let corners = triangle.map(|corner| {
            let point = mesh.vertices[corner];
            let key = [point[0].to_bits(), point[1].to_bits(), point[2].to_bits()];
            *remap.get(&key).unwrap()
        });
        for side in 0..3 {
            let a = corners[side];
            let b = corners[(side + 1) % 3];
            *edges.entry((a.min(b), a.max(b))).or_default() += 1;
        }
    }
    let mut histogram = BTreeMap::new();
    for count in edges.values() {
        *histogram.entry(*count).or_insert(0usize) += 1;
    }
    println!(
        "histogram {histogram:?} vertices {} triangles {}",
        vertices.len(),
        mesh.triangles.len()
    );
    assert!(edges.values().all(|count| *count == 2));
    assert_eq!(
        vertices.len() as i64 - edges.len() as i64 + mesh.triangles.len() as i64,
        2
    );
}
