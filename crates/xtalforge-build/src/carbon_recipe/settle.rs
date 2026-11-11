//! Spring relaxation of a periodic three-coordinate carbon graph (Task 6.8.5).
//!
//! The dual of a remeshed sheet puts one carbon at each patch centre, so the
//! raw bond lengths follow the mesh and two carbons can sit far closer than a
//! C-C bond. This pass keeps the topology and moves the carbons with a
//! gradient descent on three springs: every bond towards 1.42 Angstrom, every
//! pair of bonds sharing a carbon towards the 2.46 Angstrom 1-3 distance of a
//! 120 degree angle, and every other pair closer than 2.2 Angstrom apart.
//! The cell is fixed and the result is deterministic.

#![allow(clippy::indexing_slicing)] // bounded graph arrays

use std::collections::BTreeSet;

use super::graph::{CarbonGraph, GraphBond};
use super::net::{add, norm, sub, transform};

/// Graphene C-C bond, Angstrom.
const BOND: f64 = 1.42;
/// 1-3 distance of a 120 degree angle between two graphene bonds, Angstrom.
const SECOND: f64 = 2.459;
/// Non-bonded pairs closer than this repel, Angstrom.
const CONTACT: f64 = 2.2;
/// Pairs inside this radius are kept in the contact list, Angstrom.
const LIST_RADIUS: f64 = 4.0;
/// Gradient-descent iterations.
const ITERATIONS: usize = 400;
/// Iterations between two rebuilds of the contact list.
const LIST_EVERY: usize = 50;
/// Largest displacement of one carbon in one step, Angstrom.
const MAX_STEP: f64 = 0.1;
/// Gradient step.
const RATE: f64 = 0.25;

/// A pair term `(first, second, image of second, target, stiffness)`.
type Term = (usize, usize, [i32; 3], f64, f64);

/// Relaxes the carbon positions and returns the graph with wrapped centres.
pub(crate) fn settle(graph: &CarbonGraph) -> CarbonGraph {
    let n = graph.len();
    if n == 0 {
        return graph.clone();
    }
    let matrix = graph.lattice.matrix();
    let inverse = graph.lattice.inverse();
    let mut springs: Vec<Term> = graph
        .bonds
        .iter()
        .map(|bond| (bond.a, bond.b, bond.image, BOND, 1.0))
        .collect();
    let adjacency = graph.adjacency();
    let mut excluded: BTreeSet<(usize, usize, [i32; 3])> = BTreeSet::new();
    for bond in &graph.bonds {
        excluded.insert(key(bond.a, bond.b, bond.image));
    }
    for arms in &adjacency {
        for (index, &(first, first_image)) in arms.iter().enumerate() {
            for &(second, second_image) in arms.iter().skip(index + 1) {
                let image = [0, 1, 2].map(|axis| second_image[axis] - first_image[axis]);
                if first == second && image == [0; 3] {
                    continue;
                }
                if excluded.insert(key(first, second, image)) {
                    springs.push((first, second, image, SECOND, 0.4));
                }
            }
        }
    }
    let mut frac = graph.centers.clone();
    let mut contacts: Vec<(usize, usize, [i32; 3])> = Vec::new();
    for step in 0..ITERATIONS {
        if step % LIST_EVERY == 0 {
            contacts = contact_list(&frac, matrix, &excluded);
        }
        let mut forces = vec![[0.0; 3]; n];
        let mut pull =
            |a: usize, b: usize, image: [i32; 3], target: f64, k: f64, only_close: bool| {
                let other = add(frac[b], image.map(f64::from));
                let delta = transform(sub(other, frac[a]), matrix);
                let distance = norm(delta);
                if distance < 1e-9 || (only_close && distance >= target) {
                    return;
                }
                let scale = k * (distance - target) / distance;
                let force = delta.map(|value| value * scale);
                forces[a] = add(forces[a], force);
                forces[b] = sub(forces[b], force);
            };
        for &(a, b, image, target, k) in &springs {
            pull(a, b, image, target, k, false);
        }
        for &(a, b, image) in &contacts {
            pull(a, b, image, CONTACT, 0.5, true);
        }
        for (point, force) in frac.iter_mut().zip(&forces) {
            let length = norm(*force) * RATE;
            let scale = if length > MAX_STEP {
                MAX_STEP / length
            } else {
                1.0
            };
            let moved = force.map(|value| value * RATE * scale);
            *point = add(*point, transform(moved, inverse));
        }
    }
    let homes: Vec<[i32; 3]> = frac
        .iter()
        .map(|point| point.map(|value| value.floor() as i32))
        .collect();
    let centers = frac
        .iter()
        .zip(&homes)
        .map(|(point, home)| [0, 1, 2].map(|axis| point[axis] - f64::from(home[axis])))
        .collect();
    let bonds = graph
        .bonds
        .iter()
        .map(|bond| GraphBond {
            a: bond.a,
            b: bond.b,
            image: [0, 1, 2]
                .map(|axis| bond.image[axis] + homes[bond.b][axis] - homes[bond.a][axis]),
        })
        .collect();
    CarbonGraph {
        lattice: graph.lattice,
        centers,
        normals: graph.normals.clone(),
        bonds,
    }
}

/// Canonical key of a pair and the image of its second member.
fn key(a: usize, b: usize, image: [i32; 3]) -> (usize, usize, [i32; 3]) {
    if a < b || (a == b && image >= [0; 3]) {
        (a, b, image)
    } else {
        (b, a, image.map(|value| -value))
    }
}

/// Minimum-image pairs inside [`LIST_RADIUS`] that are not bonded or 1-3.
fn contact_list(
    frac: &[[f64; 3]],
    matrix: [[f64; 3]; 3],
    excluded: &BTreeSet<(usize, usize, [i32; 3])>,
) -> Vec<(usize, usize, [i32; 3])> {
    let mut out = Vec::new();
    for a in 0..frac.len() {
        for b in (a + 1)..frac.len() {
            let raw = sub(frac[b], frac[a]);
            let image = raw.map(|value| -(value.round()) as i32);
            let delta = transform(add(raw, image.map(f64::from)), matrix);
            if norm(delta) < LIST_RADIUS && !excluded.contains(&key(a, b, image)) {
                out.push((a, b, image));
            }
        }
    }
    out
}
