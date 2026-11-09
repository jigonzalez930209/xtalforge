//! Periodic triangle mesh and its topology (Task 6.8.2).
//!
//! Vertices are wrapped into the cell; every triangle corner carries the
//! lattice image it uses, so an edge can cross a cell face and the sheet
//! stays closed. The marching tetrahedra output is welded by fractional
//! position before this type sees it.

#![allow(clippy::indexing_slicing)] // bounded mesh arrays

use std::collections::HashMap;

use xtalforge_core::{Error, Lattice};

use super::net::{add, norm, sub, transform};

/// One undirected mesh edge with the translation from `lo` to `hi`.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub(crate) struct EdgeKey {
    pub(crate) lo: usize,
    pub(crate) hi: usize,
    pub(crate) rel: [i32; 3],
}

/// A closed periodic triangle mesh.
#[derive(Debug, Clone)]
pub(crate) struct Mesh {
    pub(crate) lattice: Lattice,
    pub(crate) frac: Vec<[f64; 3]>,
    pub(crate) tri: Vec<[usize; 3]>,
    pub(crate) shift: Vec<[[i32; 3]; 3]>,
}

impl Mesh {
    /// Number of triangles.
    pub(crate) fn triangles(&self) -> usize {
        self.tri.len()
    }

    /// Cartesian corners of one triangle.
    pub(crate) fn corners(&self, index: usize) -> Option<[[f64; 3]; 3]> {
        let triangle = self.tri.get(index)?;
        let shift = self.shift.get(index)?;
        let matrix = self.lattice.matrix();
        let mut out = [[0.0; 3]; 3];
        for column in 0..3 {
            let vertex = self.frac.get(triangle[column])?;
            out[column] = add(
                self.lattice.to_cart(*vertex),
                transform(shift[column].map(f64::from), matrix),
            );
        }
        Some(out)
    }

    /// Total surface area in square Angstrom.
    pub(crate) fn area(&self) -> f64 {
        (0..self.tri.len())
            .filter_map(|index| self.corners(index))
            .map(|corners| {
                0.5 * norm(cross(
                    sub(corners[1], corners[0]),
                    sub(corners[2], corners[0]),
                ))
            })
            .sum()
    }

    /// Every side as an undirected edge, with the translation from the
    /// smaller vertex to the larger one.
    pub(crate) fn edges(&self) -> Vec<(EdgeKey, usize, usize)> {
        let mut out = Vec::with_capacity(self.tri.len() * 3);
        for (triangle, corners) in self.tri.iter().enumerate() {
            let shift = self.shift.get(triangle).copied().unwrap_or_default();
            for side in 0..3 {
                let next = (side + 1) % 3;
                let mut lo = corners[side];
                let mut hi = corners[next];
                let mut rel = [
                    shift[next][0] - shift[side][0],
                    shift[next][1] - shift[side][1],
                    shift[next][2] - shift[side][2],
                ];
                if lo > hi || (lo == hi && lexicographic_negative(rel)) {
                    std::mem::swap(&mut lo, &mut hi);
                    rel = rel.map(|value| -value);
                }
                out.push((EdgeKey { lo, hi, rel }, triangle, side));
            }
        }
        out
    }

    /// `V - E + F` on the welded mesh.
    pub(crate) fn euler(&self) -> i64 {
        let mut used = vec![false; self.frac.len()];
        for triangle in &self.tri {
            for vertex in triangle {
                if let Some(slot) = used.get_mut(*vertex) {
                    *slot = true;
                }
            }
        }
        let vertices = used.iter().filter(|value| **value).count() as i64;
        let mut keys: Vec<EdgeKey> = self.edges().into_iter().map(|(key, _, _)| key).collect();
        keys.sort_by_key(|key| (key.lo, key.hi, key.rel));
        keys.dedup();
        vertices - keys.len() as i64 + self.tri.len() as i64
    }

    /// Triangle pairs that share an edge, with the shared edge index.
    pub(crate) fn neighbours(&self) -> Vec<([usize; 2], EdgeKey)> {
        let mut by_edge: HashMap<EdgeKey, Vec<usize>> = HashMap::new();
        for (key, triangle, _side) in self.edges() {
            by_edge.entry(key).or_default().push(triangle);
        }
        let mut out = Vec::new();
        for (key, mut triangles) in by_edge {
            triangles.sort_unstable();
            triangles.dedup();
            for pair in triangles.windows(2) {
                if let [first, second] = pair {
                    out.push(([*first, *second], key));
                }
            }
        }
        out
    }

    /// Connected component label per triangle.
    pub(crate) fn components(&self) -> Vec<usize> {
        let mut parent: Vec<usize> = (0..self.tri.len()).collect();
        for (pair, _edge) in self.neighbours() {
            let first = find(&mut parent, pair[0]);
            let second = find(&mut parent, pair[1]);
            if let Some(slot) = parent.get_mut(first) {
                *slot = second;
            }
        }
        (0..self.tri.len())
            .map(|index| find(&mut parent, index))
            .collect()
    }
}

/// Welds the marching-tetrahedra output by fractional position.
pub(crate) fn weld(
    points: &[[f64; 3]],
    triangles: &[[usize; 3]],
    lattice: &Lattice,
    shape: [usize; 3],
) -> Result<Mesh, Error> {
    let scale = [
        (shape[0].max(1) * 1000) as i64,
        (shape[1].max(1) * 1000) as i64,
        (shape[2].max(1) * 1000) as i64,
    ];
    let mut lookup: HashMap<[i64; 3], usize> = HashMap::new();
    let mut frac = Vec::new();
    let mut index = vec![0usize; points.len()];
    let mut shifts = vec![[0i32; 3]; points.len()];
    for (position, point) in points.iter().enumerate() {
        let local = lattice.to_frac(*point);
        let mut q = [0i64; 3];
        let mut wrapped = [0i64; 3];
        for axis in 0..3 {
            q[axis] = (local[axis] * scale[axis] as f64).round() as i64;
            wrapped[axis] = q[axis].rem_euclid(scale[axis]);
            shifts[position][axis] = ((q[axis] - wrapped[axis]) / scale[axis]) as i32;
        }
        let slot = *lookup.entry(wrapped).or_insert_with(|| {
            frac.push([
                wrapped[0] as f64 / scale[0] as f64,
                wrapped[1] as f64 / scale[1] as f64,
                wrapped[2] as f64 / scale[2] as f64,
            ]);
            frac.len() - 1
        });
        index[position] = slot;
    }
    let mut kept_tri = Vec::with_capacity(triangles.len());
    let mut kept_shift = Vec::with_capacity(triangles.len());
    for triangle in triangles {
        let corners = triangle.map(|corner| index.get(corner).copied().unwrap_or(0));
        let corner_shift = triangle.map(|corner| shifts.get(corner).copied().unwrap_or_default());
        let same = (0..3).any(|side| {
            let next = (side + 1) % 3;
            corners[side] == corners[next] && corner_shift[side] == corner_shift[next]
        });
        if same {
            continue;
        }
        kept_tri.push(corners);
        kept_shift.push(corner_shift);
    }
    if kept_tri.is_empty() {
        return Err(Error::InvalidInput(
            "the carbon surface produced no triangle".to_owned(),
        ));
    }
    Ok(Mesh {
        lattice: *lattice,
        frac,
        tri: kept_tri,
        shift: kept_shift,
    })
}

fn lexicographic_negative(value: [i32; 3]) -> bool {
    for component in value {
        if component != 0 {
            return component < 0;
        }
    }
    false
}

fn cross(first: [f64; 3], second: [f64; 3]) -> [f64; 3] {
    [
        first[1] * second[2] - first[2] * second[1],
        first[2] * second[0] - first[0] * second[2],
        first[0] * second[1] - first[1] * second[0],
    ]
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
