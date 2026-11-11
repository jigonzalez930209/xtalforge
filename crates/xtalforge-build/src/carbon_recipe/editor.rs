//! Isotropic remeshing of a periodic triangle mesh (Task 6.8.2).
//!
//! Split edges longer than the target, collapse edges shorter than it,
//! and flip towards valence six. Every operation keeps the sheet closed
//! and oriented, so the Euler characteristic cannot change.

#![allow(clippy::indexing_slicing)] // bounded local mesh arrays

use std::collections::{BTreeMap, BTreeSet};

use xtalforge_core::{Error, Lattice};

use super::mesh::Mesh;
use super::net::{add, norm, sub};

/// An even triangulation's edge; the dual's bonds are then 1.42 A.
pub(crate) const TARGET_EDGE: f64 = 3.55;

/// A mesh editor that preserves the closed oriented manifold.
pub(crate) struct Editor {
    pub(super) lattice: Lattice,
    pub(super) pos: Vec<[f64; 3]>,
    pub(super) tris: BTreeMap<usize, [usize; 3]>,
    pub(super) shifts: BTreeMap<usize, [[i32; 3]; 3]>,
    pub(super) at: Vec<BTreeSet<usize>>,
    next: usize,
}

impl Editor {
    pub(crate) fn new(mesh: &Mesh) -> Self {
        let mut editor = Self {
            lattice: mesh.lattice,
            pos: mesh.frac.clone(),
            tris: BTreeMap::new(),
            shifts: BTreeMap::new(),
            at: vec![BTreeSet::new(); mesh.frac.len()],
            next: 0,
        };
        for (index, triangle) in mesh.tri.iter().enumerate() {
            editor.add(*triangle, mesh.shift[index]);
        }
        editor
    }

    pub(super) fn add(&mut self, corners: [usize; 3], shift: [[i32; 3]; 3]) -> usize {
        let id = self.next;
        self.next += 1;
        self.tris.insert(id, corners);
        self.shifts.insert(id, shift);
        for vertex in corners {
            if let Some(set) = self.at.get_mut(vertex) {
                set.insert(id);
            }
        }
        id
    }

    pub(super) fn remove(&mut self, id: usize) {
        if let Some(corners) = self.tris.remove(&id) {
            for vertex in corners {
                if let Some(set) = self.at.get_mut(vertex) {
                    set.remove(&id);
                }
            }
        }
        self.shifts.remove(&id);
    }

    pub(crate) fn mesh(&self) -> Mesh {
        let mut used = BTreeSet::new();
        for corners in self.tris.values() {
            used.extend(corners.iter().copied());
        }
        let mut map = BTreeMap::new();
        let mut frac = Vec::new();
        for vertex in used {
            map.insert(vertex, frac.len());
            frac.push(self.pos[vertex]);
        }
        let mut tri = Vec::with_capacity(self.tris.len());
        let mut shift = Vec::with_capacity(self.tris.len());
        for (id, corners) in &self.tris {
            let (Some(a), Some(b), Some(c)) = (
                map.get(&corners[0]).copied(),
                map.get(&corners[1]).copied(),
                map.get(&corners[2]).copied(),
            ) else {
                continue;
            };
            tri.push([a, b, c]);
            shift.push(self.shifts.get(id).copied().unwrap_or_default());
        }
        Mesh {
            lattice: self.lattice,
            frac,
            tri,
            shift,
        }
    }

    pub(super) fn valence(&self, vertex: usize) -> usize {
        self.at.get(vertex).map_or(0, BTreeSet::len)
    }

    pub(super) fn corner(&self, triangle: usize, vertex: usize) -> Option<usize> {
        self.tris
            .get(&triangle)?
            .iter()
            .position(|corner| *corner == vertex)
    }

    fn sides(&self) -> Vec<(f64, usize, usize, [i32; 3])> {
        let mut keys: BTreeSet<(usize, usize, [i32; 3])> = BTreeSet::new();
        for (triangle, corners) in &self.tris {
            let Some(shift) = self.shifts.get(triangle) else {
                continue;
            };
            for side in 0..3 {
                let mut a = corners[side];
                let mut b = corners[(side + 1) % 3];
                let mut rel = sub3(shift[(side + 1) % 3], shift[side]);
                if a > b || (a == b && lexicographic_negative(rel)) {
                    std::mem::swap(&mut a, &mut b);
                    rel = rel.map(|value| -value);
                }
                keys.insert((a, b, rel));
            }
        }
        let matrix = self.lattice.matrix();
        keys.into_iter()
            .map(|(a, b, rel)| {
                let start = self.pos[a];
                let end = add(self.pos[b], rel.map(f64::from));
                (norm(transform(sub(end, start), matrix)), a, b, rel)
            })
            .collect()
    }

    fn split_long(&mut self, longest: f64) -> usize {
        let mut done = 0;
        let mut sides = self.sides();
        sides.sort_by(|first, second| second.0.total_cmp(&first.0));
        for (length, a, b, rel) in sides {
            if length <= longest {
                break;
            }
            if let Some((triangle, side)) = self.find(a, b, rel)
                && self.split(triangle, side).is_some()
            {
                done += 1;
            }
        }
        done
    }

    fn collapse_short(&mut self, shortest: f64, longest: f64) -> usize {
        let mut done = 0;
        let mut sides = self.sides();
        sides.sort_by(|first, second| first.0.total_cmp(&second.0));
        for (length, a, b, rel) in sides {
            if length >= shortest {
                break;
            }
            if self.at.get(a).is_none_or(BTreeSet::is_empty)
                || self.at.get(b).is_none_or(BTreeSet::is_empty)
            {
                continue;
            }
            if let Some((triangle, side)) = self.find(a, b, rel)
                && self.collapse(triangle, side, longest)
            {
                done += 1;
            }
        }
        done
    }

    fn flip_to_six(&mut self) -> usize {
        let mut done = 0;
        for (_length, a, b, rel) in self.sides() {
            let Some((triangle, side)) = self.find(a, b, rel) else {
                continue;
            };
            let Some((other, other_side)) = self.opposite(triangle, side) else {
                continue;
            };
            let corners = self.tris.get(&triangle).copied().unwrap_or_default();
            let other_corners = self.tris.get(&other).copied().unwrap_or_default();
            let c = corners[(side + 2) % 3];
            let d = other_corners[(other_side + 2) % 3];
            let ends = [
                self.valence(a),
                self.valence(b),
                self.valence(c),
                self.valence(d),
            ];
            let before: i64 = ends.iter().map(|value| (*value as i64 - 6).pow(2)).sum();
            let after: i64 = [
                ends[0].saturating_sub(1),
                ends[1].saturating_sub(1),
                ends[2] + 1,
                ends[3] + 1,
            ]
            .iter()
            .map(|value| (*value as i64 - 6).pow(2))
            .sum();
            if after >= before {
                continue;
            }
            if self.flip(triangle, side) {
                done += 1;
            }
        }
        done
    }
}

/// Remeshes the sheet towards an even `TARGET_EDGE` triangulation.
pub(crate) fn remesh(
    mesh: &Mesh,
    field: &dyn super::field::LevelSet,
    level: f64,
    iterations: usize,
    target: f64,
) -> Result<Mesh, Error> {
    if !target.is_finite() || target <= 0.0 {
        return Err(Error::InvalidInput(
            "the carbon target edge must be positive and finite".to_owned(),
        ));
    }
    let aim = 0.968 * target;
    let longest = 1.3 * aim;
    let shortest = 1.05 * aim;
    let ceiling = 4.0 * mesh.area() / (3.0_f64.sqrt() / 4.0 * target * target);
    let mut current = mesh.clone();
    for _ in 0..iterations {
        let mut editor = Editor::new(&current);
        editor.split_long(longest);
        editor.collapse_short(shortest, longest);
        editor.flip_to_six();
        current = editor.mesh();
        current = super::relax::relax(&current, field, level, 0.0)?;

        let _ = ceiling;
    }
    Ok(current)
}

pub(super) fn normal(
    first: [f64; 3],
    second: [f64; 3],
    third: [f64; 3],
    matrix: [[f64; 3]; 3],
) -> [f64; 3] {
    let a = transform(first, matrix);
    let b = transform(second, matrix);
    let c = transform(third, matrix);
    let value = cross(sub(b, a), sub(c, a));
    let length = norm(value).max(1e-12);
    value.map(|component| component / length)
}

pub(super) fn sub3(first: [i32; 3], second: [i32; 3]) -> [i32; 3] {
    [
        first[0] - second[0],
        first[1] - second[1],
        first[2] - second[2],
    ]
}

pub(super) fn add3(first: [i32; 3], second: [i32; 3]) -> [i32; 3] {
    [
        first[0] + second[0],
        first[1] + second[1],
        first[2] + second[2],
    ]
}

fn lexicographic_negative(value: [i32; 3]) -> bool {
    for component in value {
        if component != 0 {
            return component < 0;
        }
    }
    false
}

fn transform(vector: [f64; 3], matrix: [[f64; 3]; 3]) -> [f64; 3] {
    [
        vector[0] * matrix[0][0] + vector[1] * matrix[1][0] + vector[2] * matrix[2][0],
        vector[0] * matrix[0][1] + vector[1] * matrix[1][1] + vector[2] * matrix[2][1],
        vector[0] * matrix[0][2] + vector[1] * matrix[1][2] + vector[2] * matrix[2][2],
    ]
}

pub(super) fn cross(first: [f64; 3], second: [f64; 3]) -> [f64; 3] {
    [
        first[1] * second[2] - first[2] * second[1],
        first[2] * second[0] - first[0] * second[2],
        first[0] * second[1] - first[1] * second[0],
    ]
}
