//! Tangential relaxation and projection of a remeshed sheet (Task 6.8.2).

#![allow(clippy::indexing_slicing)] // bounded local vertex arrays

use xtalforge_core::Error;

use super::field::LevelSet;
use super::mesh::Mesh;
use super::net::{add, dot, norm, sub, transform};

/// Moves every vertex part of the way to its neighbours, within the sheet,
/// and puts it back on the level set with a capped Newton step.
pub(crate) fn relax(
    mesh: &Mesh,
    field: &dyn LevelSet,
    level: f64,
    damping: f64,
) -> Result<Mesh, Error> {
    let mut adjacency: Vec<Vec<(usize, [i32; 3])>> = vec![Vec::new(); mesh.frac.len()];
    for (key, _, _) in mesh.edges() {
        if let Some(list) = adjacency.get_mut(key.lo) {
            list.push((key.hi, key.rel));
        }
        if let Some(list) = adjacency.get_mut(key.hi) {
            list.push((key.lo, key.rel.map(|value| -value)));
        }
    }
    let matrix = mesh.lattice.matrix();
    let inverse = mesh.lattice.inverse();
    let mut moved = Vec::with_capacity(mesh.frac.len());
    for (index, center) in mesh.frac.iter().enumerate() {
        let neighbours = adjacency.get(index).map_or(&[][..], Vec::as_slice);
        if neighbours.is_empty() {
            moved.push(*center);
            continue;
        }
        let mut average = [0.0; 3];
        for (other, image) in neighbours {
            average = add(average, add(mesh.frac[*other], image.map(f64::from)));
        }
        let count = neighbours.len() as f64;
        let average = average.map(|value| value / count);
        let move_vector = transform(sub(average, *center), matrix);
        let (_, gradient) = field.value_gradient(*center);
        let normal = gradient.map(|value| value / norm(gradient).max(1e-12));
        let tangent = sub(
            move_vector,
            normal.map(|value| value * dot(move_vector, normal)),
        );
        let tangent = transform(tangent, inverse);
        let mut point = add(*center, tangent.map(|value| value * damping));
        for _ in 0..4 {
            let (field_value, gradient) = field.value_gradient(point);
            let denominator = dot(gradient, gradient).max(1e-12);
            let mut step =
                gradient.map(|component| component * (field_value - level) / denominator);
            let length = norm(step);
            if length > 0.5 {
                step = step.map(|value| value * 0.5 / length);
            }
            let step = transform(step, inverse);
            point = sub(point, step);
        }
        moved.push(point);
    }
    let home: Vec<[i32; 3]> = moved
        .iter()
        .map(|point| point.map(f64::floor) as [f64; 3])
        .map(|home| [home[0] as i32, home[1] as i32, home[2] as i32])
        .collect();
    let frac: Vec<[f64; 3]> = moved
        .iter()
        .zip(home.iter())
        .map(|(point, home)| {
            [
                point[0] - home[0] as f64,
                point[1] - home[1] as f64,
                point[2] - home[2] as f64,
            ]
        })
        .collect();
    let mut shift = mesh.shift.clone();
    for (triangle, corners) in mesh.tri.iter().enumerate() {
        for column in 0..3 {
            let vertex = corners[column];
            shift[triangle][column] = [
                shift[triangle][column][0] + home[vertex][0],
                shift[triangle][column][1] + home[vertex][1],
                shift[triangle][column][2] + home[vertex][2],
            ];
        }
    }
    Ok(Mesh {
        lattice: mesh.lattice,
        frac,
        tri: mesh.tri.clone(),
        shift,
    })
}
