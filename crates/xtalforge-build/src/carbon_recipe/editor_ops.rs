//! Split, collapse and flip operations of the carbon remesher (Task 6.8.2).

#![allow(clippy::indexing_slicing)] // bounded local mesh arrays

use std::collections::{BTreeMap, BTreeSet};

use super::editor::{Editor, add3, cross, normal, sub3};
use super::net::{add, dot, norm, sub, transform};

impl Editor {
    pub(super) fn neighbours(&self, vertex: usize) -> BTreeSet<(usize, [i32; 3])> {
        let mut out = BTreeSet::new();
        let Some(triangles) = self.at.get(vertex) else {
            return out;
        };
        for triangle in triangles {
            let Some(corner) = self.corner(*triangle, vertex) else {
                continue;
            };
            let Some(corners) = self.tris.get(triangle) else {
                continue;
            };
            let Some(shift) = self.shifts.get(triangle) else {
                continue;
            };
            let base = shift[corner];
            for step in [1usize, 2] {
                let image = sub3(shift[(corner + step) % 3], base);
                out.insert((corners[(corner + step) % 3], image));
            }
        }
        out
    }

    pub(super) fn opposite(&self, triangle: usize, side: usize) -> Option<(usize, usize)> {
        let corners = self.tris.get(&triangle)?;
        let shift = self.shifts.get(&triangle)?;
        let a = corners[side];
        let b = corners[(side + 1) % 3];
        let rel = sub3(shift[(side + 1) % 3], shift[side]);
        for other in self.at.get(a)?.intersection(self.at.get(b)?) {
            if *other == triangle {
                continue;
            }
            let other_corners = self.tris.get(other)?;
            let other_shift = self.shifts.get(other)?;
            for m in 0..3 {
                if other_corners[m] == b
                    && other_corners[(m + 1) % 3] == a
                    && sub3(other_shift[(m + 1) % 3], other_shift[m]) == [-rel[0], -rel[1], -rel[2]]
                {
                    return Some((*other, m));
                }
            }
        }
        None
    }

    pub(super) fn find(&self, a: usize, b: usize, rel: [i32; 3]) -> Option<(usize, usize)> {
        for triangle in self.at.get(a)? {
            let corners = self.tris.get(triangle)?;
            let shift = self.shifts.get(triangle)?;
            for side in 0..3 {
                if corners[side] == a
                    && corners[(side + 1) % 3] == b
                    && sub3(shift[(side + 1) % 3], shift[side]) == rel
                {
                    return Some((*triangle, side));
                }
            }
        }
        None
    }

    pub(super) fn split(&mut self, triangle: usize, side: usize) -> Option<usize> {
        let (other, other_side) = self.opposite(triangle, side)?;
        let corners = self.tris.get(&triangle).copied()?;
        let shift = self.shifts.get(&triangle).copied()?;
        let other_corners = self.tris.get(&other).copied()?;
        let other_shift = self.shifts.get(&other).copied()?;
        let a = corners[side];
        let b = corners[(side + 1) % 3];
        let c = corners[(side + 2) % 3];
        let d = other_corners[(other_side + 2) % 3];
        let sa = shift[side];
        let sb = shift[(side + 1) % 3];
        let sc = shift[(side + 2) % 3];
        let frame = sub3(shift[side], other_shift[(other_side + 1) % 3]);
        let sd = add3(other_shift[(other_side + 2) % 3], frame);
        let middle = [
            0.5 * (self.pos[a][0] + sa[0] as f64 + self.pos[b][0] + sb[0] as f64),
            0.5 * (self.pos[a][1] + sa[1] as f64 + self.pos[b][1] + sb[1] as f64),
            0.5 * (self.pos[a][2] + sa[2] as f64 + self.pos[b][2] + sb[2] as f64),
        ];
        let home = middle.map(f64::floor);
        let vertex = self.pos.len();
        self.pos.push([
            middle[0] - home[0],
            middle[1] - home[1],
            middle[2] - home[2],
        ]);
        self.at.push(BTreeSet::new());
        let sm = home.map(|value| value as i32);
        self.remove(triangle);
        self.remove(other);
        self.add([a, vertex, c], [sa, sm, sc]);
        self.add([vertex, b, c], [sm, sb, sc]);
        self.add([b, vertex, d], [sb, sm, sd]);
        self.add([vertex, a, d], [sm, sa, sd]);
        Some(vertex)
    }

    pub(super) fn collapse(&mut self, triangle: usize, side: usize, longest: f64) -> bool {
        let Some((other, other_side)) = self.opposite(triangle, side) else {
            return false;
        };
        let corners = self.tris.get(&triangle).copied().unwrap_or_default();
        let shift = self.shifts.get(&triangle).copied().unwrap_or_default();
        let other_corners = self.tris.get(&other).copied().unwrap_or_default();
        let other_shift = self.shifts.get(&other).copied().unwrap_or_default();
        let a = corners[side];
        let b = corners[(side + 1) % 3];
        if a == b {
            return false;
        }
        let c = corners[(side + 2) % 3];
        let d = other_corners[(other_side + 2) % 3];
        if self.valence(c) <= 3 || self.valence(d) <= 3 {
            return false;
        }
        let sa = shift[side];
        let sb = shift[(side + 1) % 3];
        let frame = sub3(shift[side], other_shift[(other_side + 1) % 3]);
        let sd = add3(other_shift[(other_side + 2) % 3], frame);
        let sc = shift[(side + 2) % 3];
        let delta = sub3(sb, sa);
        let mut around_a = self.neighbours(a);
        around_a.remove(&(b, delta));
        let around_b: BTreeSet<(usize, [i32; 3])> = self
            .neighbours(b)
            .into_iter()
            .map(|(vertex, image)| (vertex, add3(image, delta)))
            .collect();
        let mut around_b = around_b;
        around_b.remove(&(a, [0, 0, 0]));
        let common: BTreeSet<(usize, [i32; 3])> =
            around_a.intersection(&around_b).copied().collect();
        let expected: BTreeSet<(usize, [i32; 3])> = [(c, sc), (d, sd)]
            .into_iter()
            .map(|(vertex, image)| (vertex, sub3(image, sa)))
            .collect();
        if common != expected {
            return false;
        }
        let pa = add(self.pos[a], sa.map(f64::from));
        let pb = add(self.pos[b], sb.map(f64::from));
        let middle = [
            0.5 * (pa[0] + pb[0]),
            0.5 * (pa[1] + pb[1]),
            0.5 * (pa[2] + pb[2]),
        ];
        let home = middle.map(f64::floor);
        let merged = [
            middle[0] - home[0],
            middle[1] - home[1],
            middle[2] - home[2],
        ];
        // Every triangle round either end, with the image its corner needs
        // after the move.
        let mut moved: BTreeMap<usize, (usize, [i32; 3])> = BTreeMap::new();
        for (vertex, start) in [(a, pa), (b, pb)] {
            let Some(triangles) = self.at.get(vertex) else {
                continue;
            };
            for other_triangle in triangles.clone() {
                if other_triangle == triangle || other_triangle == other {
                    continue;
                }
                let Some(corner) = self.corner(other_triangle, vertex) else {
                    continue;
                };
                let Some(other_shift) = self.shifts.get(&other_triangle).copied() else {
                    continue;
                };
                let here = add(
                    add(self.pos[vertex], other_shift[corner].map(f64::from)),
                    sub(middle, start),
                );
                let image = [
                    (here[0] - merged[0]).round() as i32,
                    (here[1] - merged[1]).round() as i32,
                    (here[2] - merged[2]).round() as i32,
                ];
                moved.insert(other_triangle, (corner, image));
            }
        }
        let matrix = self.lattice.matrix();
        for (other_triangle, (corner, image)) in &moved {
            let Some(other_corners) = self.tris.get(other_triangle).copied() else {
                continue;
            };
            let Some(other_shift) = self.shifts.get(other_triangle).copied() else {
                continue;
            };
            let mut before = [[0.0; 3]; 3];
            for corner_index in 0..3 {
                let vertex = other_corners[corner_index];
                before[corner_index] =
                    add(self.pos[vertex], other_shift[corner_index].map(f64::from));
            }
            let mut after = before;
            after[*corner] = add(merged, image.map(f64::from));
            let old = cross(
                sub(transform(before[1], matrix), transform(before[0], matrix)),
                sub(transform(before[2], matrix), transform(before[0], matrix)),
            );
            let new = cross(
                sub(transform(after[1], matrix), transform(after[0], matrix)),
                sub(transform(after[2], matrix), transform(after[0], matrix)),
            );
            if norm(old) < 1e-18 || norm(new) < 1e-18 {
                return false;
            }
            if dot(old, new) / (norm(old) * norm(new)) <= 0.2 {
                return false;
            }
            for side in 0..3 {
                let start = after[side];
                let end = after[(side + 1) % 3];
                if norm(transform(sub(end, start), matrix)) > longest {
                    return false;
                }
            }
            before = after;
            let _ = before;
        }
        self.remove(triangle);
        self.remove(other);
        if let Some(slot) = self.pos.get_mut(a) {
            *slot = merged;
        }
        for (other_triangle, (corner, image)) in moved {
            let Some(other_corners) = self.tris.get(&other_triangle).copied() else {
                continue;
            };
            let Some(mut other_shift) = self.shifts.get(&other_triangle).copied() else {
                continue;
            };
            let mut corners = other_corners;
            corners[corner] = a;
            other_shift[corner] = image;
            self.remove(other_triangle);
            self.add(corners, other_shift);
        }
        true
    }

    pub(super) fn flip(&mut self, triangle: usize, side: usize) -> bool {
        let Some((other, other_side)) = self.opposite(triangle, side) else {
            return false;
        };
        let corners = self.tris.get(&triangle).copied().unwrap_or_default();
        let shift = self.shifts.get(&triangle).copied().unwrap_or_default();
        let other_corners = self.tris.get(&other).copied().unwrap_or_default();
        let other_shift = self.shifts.get(&other).copied().unwrap_or_default();
        let a = corners[side];
        let b = corners[(side + 1) % 3];
        let c = corners[(side + 2) % 3];
        let d = other_corners[(other_side + 2) % 3];
        if self.valence(a) <= 3 || self.valence(b) <= 3 {
            return false;
        }
        let sa = shift[side];
        let sb = shift[(side + 1) % 3];
        let sc = shift[(side + 2) % 3];
        let frame = sub3(shift[side], other_shift[(other_side + 1) % 3]);
        let sd = add3(other_shift[(other_side + 2) % 3], frame);
        if (c == d && sc == sd) || self.neighbours(c).contains(&(d, sub3(sd, sc))) {
            return false;
        }
        let pa = add(self.pos[a], sa.map(f64::from));
        let pb = add(self.pos[b], sb.map(f64::from));
        let pc = add(self.pos[c], sc.map(f64::from));
        let pd = add(self.pos[d], sd.map(f64::from));
        let matrix = self.lattice.matrix();
        let before = add(normal(pa, pb, pc, matrix), normal(pb, pa, pd, matrix));
        if dot(normal(pa, pd, pc, matrix), before) <= 0.2
            || dot(normal(pd, pb, pc, matrix), before) <= 0.2
        {
            return false;
        }
        self.remove(triangle);
        self.remove(other);
        self.add([a, d, c], [sa, sd, sc]);
        self.add([d, b, c], [sd, sb, sc]);
        true
    }
}
