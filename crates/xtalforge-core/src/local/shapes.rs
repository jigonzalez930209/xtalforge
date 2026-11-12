#![allow(clippy::indexing_slicing)] // fixed-size numeric kernels

//! Continuous shape measures against ideal polyhedra (Task 9.3.2).

use serde::{Deserialize, Serialize};

/// One continuous shape measure result.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct ShapeMeasure {
    /// Ideal shape name.
    pub shape: String,
    /// Continuous shape measure; zero is the ideal polyhedron.
    pub cshm: f64,
}

/// Continuous shape measures for a set of neighbour vectors around the
/// origin, compared with the ideal polyhedra of the same size.
#[must_use]
pub fn shape_measures(points: &[[f64; 3]]) -> Vec<ShapeMeasure> {
    let count = points.len();
    if count < 2 {
        return Vec::new();
    }
    let mut results = Vec::new();
    for (name, ideal) in ideal_shapes(count) {
        if let Some(cshm) = continuous_shape_measure(points, &ideal) {
            results.push(ShapeMeasure {
                shape: name.to_owned(),
                cshm,
            });
        }
    }
    results.sort_by(|first, second| {
        first
            .cshm
            .partial_cmp(&second.cshm)
            .unwrap_or(std::cmp::Ordering::Equal)
            .then_with(|| first.shape.cmp(&second.shape))
    });
    results
}

/// Continuous shape measure of a point set against one ideal polyhedron.
///
/// The value is `100 * sum |q - p|^2 / sum |p|^2` after the best rotation and
/// vertex assignment, which is zero for the ideal shape.
fn continuous_shape_measure(points: &[[f64; 3]], ideal: &[[f64; 3]]) -> Option<f64> {
    if points.len() != ideal.len() || points.len() < 2 {
        return None;
    }
    let actual = normalise(points)?;
    let target = normalise(ideal)?;
    let mut assignment: Vec<usize> = (0..target.len()).collect();
    let mut best = f64::INFINITY;
    for _ in 0..4 {
        let aligned: Vec<[f64; 3]> = assignment
            .iter()
            .map(|&index| *actual.get(index).unwrap_or(&[0.0; 3]))
            .collect();
        let alignment = crate::compare::kabsch_alignment(&aligned, &target).ok()?;
        let mut cost = vec![vec![0.0; target.len()]; actual.len()];
        for (row, point) in actual.iter().enumerate() {
            let moved = apply(alignment.rotation, *point);
            for (column, ideal_point) in target.iter().enumerate() {
                let delta = [
                    moved[0] - ideal_point[0],
                    moved[1] - ideal_point[1],
                    moved[2] - ideal_point[2],
                ];
                if let Some(entry) = cost.get_mut(row).and_then(|r| r.get_mut(column)) {
                    *entry = delta[0] * delta[0] + delta[1] * delta[1] + delta[2] * delta[2];
                }
            }
        }
        assignment = hungarian(&cost);
        let total: f64 = assignment
            .iter()
            .enumerate()
            .filter_map(|(row, &column)| cost.get(row)?.get(column).copied())
            .sum();
        if total >= best - 1e-12 {
            best = best.min(total);
            break;
        }
        best = total;
    }
    let denominator: f64 = target.iter().map(|point| squared(*point)).sum();
    if denominator < 1e-18 {
        return None;
    }
    Some(100.0 * best / denominator)
}

fn normalise(points: &[[f64; 3]]) -> Option<Vec<[f64; 3]>> {
    let mut centroid = [0.0; 3];
    for point in points {
        centroid[0] += point[0];
        centroid[1] += point[1];
        centroid[2] += point[2];
    }
    let count = points.len() as f64;
    centroid = [
        centroid[0] / count,
        centroid[1] / count,
        centroid[2] / count,
    ];
    let mut centred: Vec<[f64; 3]> = points
        .iter()
        .map(|point| {
            [
                point[0] - centroid[0],
                point[1] - centroid[1],
                point[2] - centroid[2],
            ]
        })
        .collect();
    let rms = (centred.iter().map(|point| squared(*point)).sum::<f64>() / count).sqrt();
    if rms < 1e-12 {
        return None;
    }
    for point in &mut centred {
        *point = [point[0] / rms, point[1] / rms, point[2] / rms];
    }
    Some(centred)
}

fn ideal_shapes(count: usize) -> Vec<(&'static str, Vec<[f64; 3]>)> {
    let mut shapes = Vec::new();
    match count {
        2 => shapes.push(("linear", vec![[1.0, 0.0, 0.0], [-1.0, 0.0, 0.0]])),
        3 => shapes.push((
            "trigonal planar",
            vec![
                [1.0, 0.0, 0.0],
                [-0.5, 0.866_025_403_784, 0.0],
                [-0.5, -0.866_025_403_784, 0.0],
            ],
        )),
        4 => {
            shapes.push((
                "tetrahedral",
                vec![
                    [1.0, 1.0, 1.0],
                    [1.0, -1.0, -1.0],
                    [-1.0, 1.0, -1.0],
                    [-1.0, -1.0, 1.0],
                ],
            ));
            shapes.push((
                "square planar",
                vec![
                    [1.0, 0.0, 0.0],
                    [0.0, 1.0, 0.0],
                    [-1.0, 0.0, 0.0],
                    [0.0, -1.0, 0.0],
                ],
            ));
        }
        5 => {
            shapes.push((
                "trigonal bipyramid",
                vec![
                    [0.0, 0.0, 1.0],
                    [0.0, 0.0, -1.0],
                    [1.0, 0.0, 0.0],
                    [-0.5, 0.866_025_403_784, 0.0],
                    [-0.5, -0.866_025_403_784, 0.0],
                ],
            ));
            shapes.push((
                "square pyramid",
                vec![
                    [0.0, 0.0, 1.0],
                    [1.0, 0.0, 0.0],
                    [0.0, 1.0, 0.0],
                    [-1.0, 0.0, 0.0],
                    [0.0, -1.0, 0.0],
                ],
            ));
        }
        6 => shapes.push((
            "octahedral",
            vec![
                [1.0, 0.0, 0.0],
                [-1.0, 0.0, 0.0],
                [0.0, 1.0, 0.0],
                [0.0, -1.0, 0.0],
                [0.0, 0.0, 1.0],
                [0.0, 0.0, -1.0],
            ],
        )),
        7 => shapes.push((
            "pentagonal bipyramid",
            vec![
                [0.0, 0.0, 1.0],
                [0.0, 0.0, -1.0],
                [1.0, 0.0, 0.0],
                [0.309_016_994_375, 0.951_056_516_295, 0.0],
                [-0.809_016_994_375, 0.587_785_252_292, 0.0],
                [-0.809_016_994_375, -0.587_785_252_292, 0.0],
                [0.309_016_994_375, -0.951_056_516_295, 0.0],
            ],
        )),
        8 => {
            let mut cube = Vec::new();
            for x in [-1.0, 1.0] {
                for y in [-1.0, 1.0] {
                    for z in [-1.0, 1.0] {
                        cube.push([x, y, z]);
                    }
                }
            }
            shapes.push(("cube", cube));
            shapes.push((
                "square antiprism",
                vec![
                    [1.0, 0.0, 0.5],
                    [0.0, 1.0, 0.5],
                    [-1.0, 0.0, 0.5],
                    [0.0, -1.0, 0.5],
                    [
                        std::f64::consts::FRAC_1_SQRT_2,
                        std::f64::consts::FRAC_1_SQRT_2,
                        -0.5,
                    ],
                    [
                        -std::f64::consts::FRAC_1_SQRT_2,
                        std::f64::consts::FRAC_1_SQRT_2,
                        -0.5,
                    ],
                    [
                        -std::f64::consts::FRAC_1_SQRT_2,
                        -std::f64::consts::FRAC_1_SQRT_2,
                        -0.5,
                    ],
                    [
                        std::f64::consts::FRAC_1_SQRT_2,
                        -std::f64::consts::FRAC_1_SQRT_2,
                        -0.5,
                    ],
                ],
            ));
        }
        12 => {
            let golden = (1.0 + 5.0_f64.sqrt()) / 2.0;
            let mut icosahedron = Vec::new();
            for first in [-1.0, 1.0] {
                for second in [-golden, golden] {
                    icosahedron.push([0.0, first, second]);
                    icosahedron.push([first, second, 0.0]);
                    icosahedron.push([second, 0.0, first]);
                }
            }
            shapes.push(("icosahedral", icosahedron));
        }
        _ => {}
    }
    shapes
}

/// Hungarian assignment for a square cost matrix (minimization).
fn hungarian(cost: &[Vec<f64>]) -> Vec<usize> {
    let n = cost.len();
    if n == 0 {
        return Vec::new();
    }
    let mut u = vec![0.0; n + 1];
    let mut v = vec![0.0; n + 1];
    let mut p = vec![0usize; n + 1];
    let mut way = vec![0usize; n + 1];
    for i in 1..=n {
        p[0] = i;
        let mut j0 = 0usize;
        let mut minimum = vec![f64::INFINITY; n + 1];
        let mut used = vec![false; n + 1];
        loop {
            used[j0] = true;
            let i0 = p[j0];
            let mut delta = f64::INFINITY;
            let mut j1 = 0usize;
            for j in 1..=n {
                if used[j] {
                    continue;
                }
                let Some(&entry) = cost.get(i0 - 1).and_then(|row| row.get(j - 1)) else {
                    continue;
                };
                let current = entry - u[i0] - v[j];
                if current < minimum[j] {
                    minimum[j] = current;
                    way[j] = j0;
                }
                if minimum[j] < delta {
                    delta = minimum[j];
                    j1 = j;
                }
            }
            for j in 0..=n {
                if used[j] {
                    u[p[j]] += delta;
                    v[j] -= delta;
                } else {
                    minimum[j] -= delta;
                }
            }
            j0 = j1;
            if p[j0] == 0 {
                break;
            }
        }
        loop {
            let j1 = way[j0];
            p[j0] = p[j1];
            j0 = j1;
            if j0 == 0 {
                break;
            }
        }
    }
    let mut assignment = vec![0usize; n];
    for j in 1..=n {
        if let Some(&row) = p.get(j)
            && row > 0
        {
            assignment[row - 1] = j - 1;
        }
    }
    assignment
}

fn apply(matrix: [[f64; 3]; 3], point: [f64; 3]) -> [f64; 3] {
    [
        matrix[0][0] * point[0] + matrix[0][1] * point[1] + matrix[0][2] * point[2],
        matrix[1][0] * point[0] + matrix[1][1] * point[1] + matrix[1][2] * point[2],
        matrix[2][0] * point[0] + matrix[2][1] * point[1] + matrix[2][2] * point[2],
    ]
}

fn squared(point: [f64; 3]) -> f64 {
    point[0] * point[0] + point[1] * point[1] + point[2] * point[2]
}
