//! Carbon graph from the periodic triangle mesh (Task 6.8.2).
//!
//! The mesh is clustered into patches of roughly one graphene cell
//! (`2.46^2 sqrt(3)/2` square Angstrom); one carbon sits at each patch
//! centre and two patches that share a mesh edge are bonded. This is the
//! dual of the meshed sheet, with the clustering standing in for the
//! reference's even remesh.

#![allow(clippy::indexing_slicing)] // bounded graph arrays

use std::collections::{BTreeMap, VecDeque};

use xtalforge_core::{Error, Lattice};

use super::mesh::Mesh;
use super::net::{norm, sub, transform};

/// One bond of the carbon graph, with the lattice image of `b`.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub(crate) struct GraphBond {
    /// First carbon.
    pub(crate) a: usize,
    /// Second carbon.
    pub(crate) b: usize,
    /// Lattice image of `b`.
    pub(crate) image: [i32; 3],
}

/// A periodic three-coordinate carbon graph.
#[derive(Debug, Clone)]
pub(crate) struct CarbonGraph {
    /// Cell.
    pub(crate) lattice: Lattice,
    /// Wrapped fractional centres.
    pub(crate) centers: Vec<[f64; 3]>,
    /// Outward unit normal per carbon.
    pub(crate) normals: Vec<[f64; 3]>,
    /// Bonds.
    pub(crate) bonds: Vec<GraphBond>,
}

impl CarbonGraph {
    /// Number of carbons.
    pub(crate) fn len(&self) -> usize {
        self.centers.len()
    }

    /// Degree of every carbon.
    pub(crate) fn degree(&self) -> Vec<usize> {
        let mut degree = vec![0usize; self.centers.len()];
        for bond in &self.bonds {
            if let Some(value) = degree.get_mut(bond.a) {
                *value += 1;
            }
            if let Some(value) = degree.get_mut(bond.b) {
                *value += 1;
            }
        }
        degree
    }

    /// Connected component label per carbon.
    pub(crate) fn components(&self) -> Vec<usize> {
        let mut parent: Vec<usize> = (0..self.centers.len()).collect();
        for bond in &self.bonds {
            let first = find(&mut parent, bond.a);
            let second = find(&mut parent, bond.b);
            if let Some(slot) = parent.get_mut(first) {
                *slot = second;
            }
        }
        (0..self.centers.len())
            .map(|index| find(&mut parent, index))
            .collect()
    }

    /// The biggest connected component, alone.
    pub(crate) fn largest_piece(&self) -> Self {
        let labels = self.components();
        let mut counts: BTreeMap<usize, usize> = BTreeMap::new();
        for label in &labels {
            *counts.entry(*label).or_default() += 1;
        }
        let largest = counts
            .into_iter()
            .max_by_key(|(label, count)| (*count, usize::MAX - *label))
            .map_or(0, |(label, _)| label);
        let keep: Vec<bool> = labels.iter().map(|label| *label == largest).collect();
        self.subset(&keep)
    }

    /// Keeps the carbons `keep` marks.
    pub(crate) fn subset(&self, keep: &[bool]) -> Self {
        let mut map = vec![usize::MAX; self.centers.len()];
        let mut centers = Vec::new();
        let mut normals = Vec::new();
        for (index, selected) in keep.iter().enumerate() {
            if !*selected {
                continue;
            }
            if let Some(slot) = map.get_mut(index) {
                *slot = centers.len();
            }
            centers.push(self.centers[index]);
            normals.push(self.normals[index]);
        }
        let mut bonds = Vec::new();
        for bond in &self.bonds {
            let Some(a) = map.get(bond.a).copied() else {
                continue;
            };
            let Some(b) = map.get(bond.b).copied() else {
                continue;
            };
            if a == usize::MAX || b == usize::MAX {
                continue;
            }
            bonds.push(GraphBond {
                a,
                b,
                image: bond.image,
            });
        }
        Self {
            lattice: self.lattice,
            centers,
            normals,
            bonds,
        }
    }

    /// Drops carbons with fewer than two neighbours, repeatedly.
    pub(crate) fn prune(&self) -> Self {
        let mut graph = self.clone();
        loop {
            let degree = graph.degree();
            let keep: Vec<bool> = degree.iter().map(|value| *value >= 2).collect();
            if keep.iter().all(|value| *value) {
                return graph;
            }
            graph = graph.subset(&keep);
        }
    }

    /// Number of independent directions the graph runs through the cell.
    pub(crate) fn periodicity(&self) -> usize {
        let adjacency = self.adjacency();
        let mut where_is: BTreeMap<usize, [i32; 3]> = BTreeMap::new();
        let mut loops: Vec<[i32; 3]> = Vec::new();
        for start in 0..self.centers.len() {
            if where_is.contains_key(&start) {
                continue;
            }
            where_is.insert(start, [0; 3]);
            let mut stack = vec![start];
            while let Some(current) = stack.pop() {
                let current_image = where_is.get(&current).copied().unwrap_or([0; 3]);
                let Some(neighbours) = adjacency.get(current) else {
                    continue;
                };
                for (next, image) in neighbours {
                    let at = [
                        current_image[0] + image[0],
                        current_image[1] + image[1],
                        current_image[2] + image[2],
                    ];
                    if let Some(existing) = where_is.get(next) {
                        if *existing != at {
                            loops.push([
                                at[0] - existing[0],
                                at[1] - existing[1],
                                at[2] - existing[2],
                            ]);
                            if loops.len() > 64 && rank3(&loops) == 3 {
                                return 3;
                            }
                        }
                    } else {
                        where_is.insert(*next, at);
                        stack.push(*next);
                    }
                }
            }
        }
        rank3(&loops)
    }

    /// Directed adjacency with bond images.
    pub(crate) fn adjacency(&self) -> Vec<Vec<(usize, [i32; 3])>> {
        let mut out = vec![Vec::new(); self.centers.len()];
        for bond in &self.bonds {
            if let Some(list) = out.get_mut(bond.a) {
                list.push((bond.b, bond.image));
            }
            if let Some(list) = out.get_mut(bond.b) {
                list.push((bond.a, bond.image.map(|value| -value)));
            }
        }
        out
    }
}

/// One carbon at every triangle centre, bonded across every shared edge.
pub(crate) fn dualise(mesh: &Mesh) -> Result<CarbonGraph, Error> {
    let triangles = mesh.triangles();
    if triangles < 4 {
        return Err(Error::InvalidInput(
            "the carbon mesh has too few triangles to dualise".to_owned(),
        ));
    }
    let matrix = mesh.lattice.matrix();
    let mut centers = Vec::with_capacity(triangles);
    let mut unwrapped = Vec::with_capacity(triangles);
    let mut normals = Vec::with_capacity(triangles);
    for index in 0..triangles {
        let Some(corners) = mesh.corners(index) else {
            continue;
        };
        let center = [
            (corners[0][0] + corners[1][0] + corners[2][0]) / 3.0,
            (corners[0][1] + corners[1][1] + corners[2][1]) / 3.0,
            (corners[0][2] + corners[1][2] + corners[2][2]) / 3.0,
        ];
        let normal = cross(sub(corners[1], corners[0]), sub(corners[2], corners[0]));
        let length = norm(normal).max(1e-12);
        let local = mesh.lattice.to_frac(center);
        unwrapped.push(local);
        centers.push(local.map(|value| value.rem_euclid(1.0)));
        normals.push(normal.map(|component| component / length));
    }
    let _ = matrix;
    let mut bonds: BTreeMap<(usize, usize), GraphBond> = BTreeMap::new();
    for (pair, _edge) in mesh.neighbours() {
        if pair[0] == pair[1] {
            continue;
        }
        let Some(first) = unwrapped.get(pair[0]).copied() else {
            continue;
        };
        let Some(second) = unwrapped.get(pair[1]).copied() else {
            continue;
        };
        let wrapped_first = centers.get(pair[0]).copied().unwrap_or(first);
        let wrapped_second = centers.get(pair[1]).copied().unwrap_or(second);
        // The centre of a triangle that crosses a cell face is an average
        // in several images, so choose the bond image by the shortest
        // Cartesian contact rather than by that average.
        let matrix = mesh.lattice.matrix();
        let mut best = [0i32; 3];
        let mut best_length = f64::INFINITY;
        for image_x in -2..=2i32 {
            for image_y in -2..=2i32 {
                for image_z in -2..=2i32 {
                    let image = [image_x, image_y, image_z];
                    let end = [
                        wrapped_second[0] + f64::from(image[0]),
                        wrapped_second[1] + f64::from(image[1]),
                        wrapped_second[2] + f64::from(image[2]),
                    ];
                    let length = norm(transform(sub(end, wrapped_first), matrix));
                    if length < best_length {
                        best_length = length;
                        best = image;
                    }
                }
            }
        }
        let image = best;
        let (a, b, image) = if pair[0] < pair[1] {
            (pair[0], pair[1], image)
        } else {
            (pair[1], pair[0], image.map(|value| -value))
        };
        bonds.entry((a, b)).or_insert(GraphBond { a, b, image });
    }
    let bonds: Vec<GraphBond> = bonds.into_values().collect();
    let graph = CarbonGraph {
        lattice: mesh.lattice,
        centers,
        normals,
        bonds,
    };
    Ok(graph)
}

/// Cuts the graph into ribbons along the net edges.
///
/// Each edge gets a shortest path over the triangle graph from the anchor
/// of one net vertex to the anchor of the other; the paths share anchors,
/// so their union is the net's own network. The kept carbon is the set of
/// graph nodes within a number of steps of that network, widened until the
/// largest piece holds `target` carbons and percolates in all directions.
pub(crate) fn cut(
    graph: &CarbonGraph,
    segments: &[(usize, usize, [f64; 3], [f64; 3])],
    vertices: &[[f64; 3]],
    target: usize,
) -> Result<CarbonGraph, Error> {
    if target == 0 {
        return Err(Error::InvalidInput(
            "the carbon target must be positive".to_owned(),
        ));
    }
    if target >= graph.len() {
        return Ok(graph.clone());
    }
    let matrix = graph.lattice.matrix();
    let anchor = |point: [f64; 3]| -> usize {
        let mut best = (f64::INFINITY, 0usize);
        for (index, center) in graph.centers.iter().enumerate() {
            let cart = transform(*center, matrix);
            let mut distance = norm(sub(cart, point));
            for shift_x in -1..=1i32 {
                for shift_y in -1..=1i32 {
                    for shift_z in -1..=1i32 {
                        let image = transform(
                            [f64::from(shift_x), f64::from(shift_y), f64::from(shift_z)],
                            matrix,
                        );
                        distance = distance.min(norm(sub(
                            cart,
                            [
                                point[0] + image[0],
                                point[1] + image[1],
                                point[2] + image[2],
                            ],
                        )));
                    }
                }
            }
            if distance < best.0 {
                best = (distance, index);
            }
        }
        best.1
    };
    let anchors: Vec<usize> = vertices.iter().map(|vertex| anchor(*vertex)).collect();
    let adjacency = graph.adjacency();
    let mut on_path = vec![false; graph.len()];
    for (first, second, _start, _end) in segments {
        let Some(start) = anchors.get(*first).copied() else {
            continue;
        };
        let Some(goal) = anchors.get(*second).copied() else {
            continue;
        };
        if let Some(path) = shortest_path(&adjacency, start, goal) {
            for node in path {
                if let Some(slot) = on_path.get_mut(node) {
                    *slot = true;
                }
            }
        }
    }
    if !on_path.iter().any(|value| *value) {
        return Ok(graph.prune().largest_piece());
    }
    let mut distance = vec![usize::MAX; graph.len()];
    let mut queue = VecDeque::new();
    for (index, value) in on_path.iter().enumerate() {
        if *value {
            if let Some(slot) = distance.get_mut(index) {
                *slot = 0;
            }
            queue.push_back(index);
        }
    }
    while let Some(current) = queue.pop_front() {
        let next = distance.get(current).copied().unwrap_or(0) + 1;
        for (other, _image) in adjacency.get(current).into_iter().flatten() {
            if let Some(slot) = distance.get_mut(*other)
                && *slot == usize::MAX
            {
                *slot = next;
                queue.push_back(*other);
            }
        }
    }
    let mut order: Vec<usize> = (0..graph.len()).collect();
    order.sort_by_key(|index| distance.get(*index).copied().unwrap_or(usize::MAX));
    let path_count = on_path.iter().filter(|value| **value).count();
    let mut keep = on_path.clone();
    if path_count > target {
        let cut = graph.prune().subset(&keep).largest_piece();
        if cut.periodicity() == 3 {
            return Ok(cut);
        }
    }
    let mut cursor = 0usize;
    let mut best = graph.prune().subset(&keep).largest_piece();
    while best.len() < target && cursor < order.len() {
        let Some(index) = order.get(cursor).copied() else {
            break;
        };
        cursor += 1;
        if on_path.get(index).copied().unwrap_or(false) {
            continue;
        }
        if let Some(slot) = keep.get_mut(index) {
            *slot = true;
        }
        if cursor.is_multiple_of(8) || best.len() < target {
            best = graph.prune().subset(&keep).largest_piece();
        }
    }
    best = graph.prune().subset(&keep).largest_piece();
    if best.periodicity() == 3 {
        return Ok(best);
    }
    // The path network percolates by construction; if a poor anchor lost
    // it, widen to every node rather than return a fragment.
    Ok(graph.prune().largest_piece())
}

fn shortest_path(
    adjacency: &[Vec<(usize, [i32; 3])>],
    start: usize,
    goal: usize,
) -> Option<Vec<usize>> {
    if start == goal {
        return Some(vec![start]);
    }
    let mut parent: Vec<usize> = vec![usize::MAX; adjacency.len()];
    let mut seen = vec![false; adjacency.len()];
    let mut queue = VecDeque::from([start]);
    seen[start] = true;
    while let Some(current) = queue.pop_front() {
        for (next, _image) in adjacency.get(current).into_iter().flatten() {
            if seen.get(*next).copied().unwrap_or(true) {
                continue;
            }
            if let Some(slot) = seen.get_mut(*next) {
                *slot = true;
            }
            if let Some(slot) = parent.get_mut(*next) {
                *slot = current;
            }
            if *next == goal {
                let mut path = vec![goal];
                let mut walk = goal;
                while let Some(previous) = parent.get(walk).copied() {
                    if previous == usize::MAX {
                        break;
                    }
                    path.push(previous);
                    walk = previous;
                }
                path.reverse();
                return Some(path);
            }
            queue.push_back(*next);
        }
    }
    None
}

fn find(parent: &mut [usize], index: usize) -> usize {
    let mut current = index;
    while parent.get(current).copied().unwrap_or(current) != current {
        current = parent.get(current).copied().unwrap_or(current);
    }
    let root = current;
    let mut walk = index;
    while parent.get(walk).copied().unwrap_or(walk) != walk {
        let next = parent.get(walk).copied().unwrap_or(walk);
        if let Some(slot) = parent.get_mut(walk) {
            *slot = root;
        }
        walk = next;
    }
    root
}

fn rank3(vectors: &[[i32; 3]]) -> usize {
    let mut basis: Vec<[f64; 3]> = Vec::new();
    for vector in vectors {
        let mut value = vector.map(f64::from);
        for existing in &basis {
            let dot = value[0] * existing[0] + value[1] * existing[1] + value[2] * existing[2];
            value[0] -= dot * existing[0];
            value[1] -= dot * existing[1];
            value[2] -= dot * existing[2];
        }
        let length = (value[0] * value[0] + value[1] * value[1] + value[2] * value[2]).sqrt();
        if length > 1e-9 {
            basis.push(value.map(|component| component / length));
            if basis.len() == 3 {
                return 3;
            }
        }
    }
    basis.len()
}
fn cross(first: [f64; 3], second: [f64; 3]) -> [f64; 3] {
    [
        first[1] * second[2] - first[2] * second[1],
        first[2] * second[0] - first[0] * second[2],
        first[0] * second[1] - first[1] * second[0],
    ]
}
