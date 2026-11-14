//! Periodic connected components of the probe-accessible void (Task 17.1.2).
//!
//! The samples whose surface distance is at least the probe radius are joined
//! to their six periodic neighbours with a union-find. A component percolates
//! along `a`, `b` or `c` when it contains a sample on both opposite cell faces
//! of that axis. This is the same condition Zeo++ uses for its channel
//! dimensionality, applied to the voxelised field instead of a Voronoi
//! decomposition: a channel narrower than one grid step can be missed, so the
//! result converges as the spacing shrinks.
//!
//! Component sizes are sample counts and are returned largest first; ties are
//! broken by the flat index of the component root, so the order is stable.

use std::collections::BTreeMap;

use serde::{Deserialize, Serialize};
use xtalforge_core::Error;

use super::DistanceGrid;

/// The probe-accessible void split into periodic connected components.
#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
pub struct Percolation {
    /// Whether the accessible void percolates along `a`, `b` and `c`.
    pub directions: [bool; 3],
    /// Number of samples in each accessible component, largest first.
    pub component_sizes: Vec<usize>,
}

impl Percolation {
    /// Number of percolating directions, from 0 to 3.
    #[must_use]
    pub fn dimension(&self) -> usize {
        self.directions.iter().filter(|value| **value).count()
    }

    /// Whether the void percolates along at least one direction.
    #[must_use]
    pub fn percolates(&self) -> bool {
        self.dimension() > 0
    }
}

impl DistanceGrid {
    /// Connected components of the samples where a probe of `probe_radius`
    /// fits, and the directions they percolate through.
    ///
    /// A negative probe radius is treated as zero. All loops run in a fixed
    /// order, so the component sizes are deterministic.
    ///
    /// # Errors
    ///
    /// Returns [`Error::InvalidInput`] when the probe radius is not finite.
    pub fn percolation(&self, probe_radius: f64) -> Result<Percolation, Error> {
        let threshold = validated_probe(probe_radius)?;
        let [n_x, n_y, n_z] = self.grid.shape;
        let mut components = UnionFind::new(self.len(), self.grid.shape);
        for z in 0..n_z {
            for y in 0..n_y {
                for x in 0..n_x {
                    let Some(flat) = self.grid.index(x, y, z) else {
                        continue;
                    };
                    if self
                        .grid
                        .data
                        .get(flat)
                        .copied()
                        .unwrap_or(f64::NEG_INFINITY)
                        >= threshold
                    {
                        components.activate(flat, x, y, z);
                    }
                }
            }
        }
        Ok(components.finish())
    }

    /// Largest probe radius at which a component spans the cell in `axis`.
    ///
    /// `axis` is `0`, `1` or `2` for `a`, `b` or `c`; `None` accepts any axis.
    /// The sweep starts with the largest sample value and activates samples in
    /// decreasing order, so the first threshold at which a component touches
    /// both opposite faces is the answer. Returns zero when nothing percolates
    /// at a non-negative probe radius, and infinity for a grid without atoms.
    ///
    /// # Errors
    ///
    /// Returns [`Error::InvalidInput`] when `axis` is not `0`, `1`, `2` or
    /// `None`.
    pub fn percolation_radius(&self, axis: Option<usize>) -> Result<f64, Error> {
        if let Some(axis) = axis
            && axis >= 3
        {
            return Err(Error::InvalidInput(format!(
                "lattice axis must be 0, 1 or 2, got {axis}"
            )));
        }
        let [n_x, n_y, n_z] = self.grid.shape;
        let mut values: Vec<(f64, usize)> = self
            .grid
            .data
            .iter()
            .copied()
            .enumerate()
            .map(|(index, value)| (value, index))
            .collect();
        values.sort_by(|first, second| {
            second
                .0
                .total_cmp(&first.0)
                .then_with(|| first.1.cmp(&second.1))
        });
        let mut components = UnionFind::new(self.len(), self.grid.shape);
        let mut position = 0;
        let mut newly_active: Vec<usize> = Vec::new();
        while let Some(&(threshold, _)) = values.get(position) {
            newly_active.clear();
            while let Some(&(value, flat)) = values.get(position) {
                if value != threshold {
                    break;
                }
                let (x, y, z) = coordinates(flat, [n_x, n_y, n_z]);
                if !components.is_active(flat) {
                    components.activate(flat, x, y, z);
                    newly_active.push(flat);
                }
                position += 1;
            }
            for flat in &newly_active {
                let root = components.find(*flat);
                if components.spans(root, axis) {
                    return Ok(threshold.max(0.0));
                }
            }
        }
        Ok(0.0)
    }

    /// Pore-limiting diameter in Angstrom: twice
    /// [`percolation_radius`](DistanceGrid::percolation_radius) over any axis.
    ///
    /// Returns `0` when no component percolates at a non-negative probe radius
    /// and infinity for a grid without atoms.
    ///
    /// # Errors
    ///
    /// Returns the errors of [`DistanceGrid::percolation_radius`].
    pub fn pore_limiting_diameter(&self) -> Result<f64, Error> {
        Ok(2.0 * self.percolation_radius(None)?)
    }
}

fn validated_probe(probe_radius: f64) -> Result<f64, Error> {
    if !probe_radius.is_finite() {
        return Err(Error::InvalidInput(format!(
            "probe radius must be finite, got {probe_radius}"
        )));
    }
    Ok(probe_radius.max(0.0))
}

/// Six boundary bits per component: `(0, 1)` for the `a` low/high faces,
/// `(2, 3)` for `b`, `(4, 5)` for `c`.
fn boundary_bits(x: usize, y: usize, z: usize, shape: [usize; 3]) -> u8 {
    let [n_x, n_y, n_z] = shape;
    let mut bits = 0u8;
    if x == 0 {
        bits |= 1 << 0;
    }
    if x + 1 == n_x {
        bits |= 1 << 1;
    }
    if y == 0 {
        bits |= 1 << 2;
    }
    if y + 1 == n_y {
        bits |= 1 << 3;
    }
    if z == 0 {
        bits |= 1 << 4;
    }
    if z + 1 == n_z {
        bits |= 1 << 5;
    }
    bits
}

fn spans(faces: u8, axis: Option<usize>) -> bool {
    (0..3).any(|position| {
        axis.is_none_or(|axis| axis == position)
            && faces & (1 << (2 * position)) != 0
            && faces & (1 << (2 * position + 1)) != 0
    })
}

fn coordinates(flat: usize, shape: [usize; 3]) -> (usize, usize, usize) {
    let [n_x, n_y, _] = shape;
    if n_x == 0 {
        return (0, 0, 0);
    }
    let x = flat % n_x;
    let y = (flat / n_x) % n_y;
    let z = flat / n_x.saturating_mul(n_y).max(1);
    (x, y, z)
}

fn step_forward(value: usize, size: usize) -> usize {
    value
        .checked_add(1)
        .filter(|next| *next < size)
        .unwrap_or(0)
}

fn step_backward(value: usize, size: usize) -> usize {
    value
        .checked_sub(1)
        .unwrap_or_else(|| size.saturating_sub(1))
}

struct UnionFind {
    parent: Vec<usize>,
    rank: Vec<u8>,
    faces: Vec<u8>,
    active: Vec<bool>,
    shape: [usize; 3],
}

impl UnionFind {
    fn new(len: usize, shape: [usize; 3]) -> Self {
        Self {
            parent: (0..len).collect(),
            rank: vec![0; len],
            faces: vec![0; len],
            active: vec![false; len],
            shape,
        }
    }

    fn is_active(&self, flat: usize) -> bool {
        self.active.get(flat).copied().unwrap_or(false)
    }

    fn activate(&mut self, flat: usize, x: usize, y: usize, z: usize) {
        if self.is_active(flat) {
            return;
        }
        if let Some(slot) = self.active.get_mut(flat) {
            *slot = true;
        }
        let bits = boundary_bits(x, y, z, self.shape);
        if let Some(slot) = self.faces.get_mut(flat) {
            *slot = bits;
        }
        let [n_x, n_y, n_z] = self.shape;
        let neighbours = [
            (step_backward(x, n_x), y, z),
            (step_forward(x, n_x), y, z),
            (x, step_backward(y, n_y), z),
            (x, step_forward(y, n_y), z),
            (x, y, step_backward(z, n_z)),
            (x, y, step_forward(z, n_z)),
        ];
        for (n_x, n_y, n_z) in neighbours {
            let Some(neighbour) = self.index(n_x, n_y, n_z) else {
                continue;
            };
            if self.is_active(neighbour) {
                union(
                    &mut self.parent,
                    &mut self.rank,
                    &mut self.faces,
                    flat,
                    neighbour,
                );
            }
        }
    }

    fn index(&self, x: usize, y: usize, z: usize) -> Option<usize> {
        let [n_x, n_y, n_z] = self.shape;
        if x >= n_x || y >= n_y || z >= n_z {
            return None;
        }
        (z * n_y + y).checked_mul(n_x)?.checked_add(x)
    }

    fn find(&mut self, node: usize) -> usize {
        let mut root = node;
        while let Some(&parent) = self.parent.get(root) {
            if parent == root {
                break;
            }
            root = parent;
        }
        let mut current = node;
        while let Some(&parent) = self.parent.get(current) {
            if parent == root {
                break;
            }
            if let Some(slot) = self.parent.get_mut(current) {
                *slot = root;
            }
            current = parent;
        }
        root
    }

    fn spans(&self, root: usize, axis: Option<usize>) -> bool {
        spans(self.faces.get(root).copied().unwrap_or(0), axis)
    }

    fn finish(mut self) -> Percolation {
        let mut counts: BTreeMap<usize, usize> = BTreeMap::new();
        for flat in 0..self.parent.len() {
            if !self.is_active(flat) {
                continue;
            }
            let root = self.find(flat);
            *counts.entry(root).or_insert(0) += 1;
        }
        let mut directions = [false; 3];
        for root in counts.keys() {
            let faces = self.faces.get(*root).copied().unwrap_or(0);
            for (axis, found) in directions.iter_mut().enumerate() {
                if spans(faces, Some(axis)) {
                    *found = true;
                }
            }
        }
        let mut ordered: Vec<(usize, usize)> = counts.into_iter().collect();
        ordered
            .sort_by(|first, second| second.1.cmp(&first.1).then_with(|| first.0.cmp(&second.0)));
        Percolation {
            directions,
            component_sizes: ordered.into_iter().map(|(_, size)| size).collect(),
        }
    }
}

fn union(
    parent: &mut [usize],
    rank: &mut [u8],
    faces: &mut [u8],
    first: usize,
    second: usize,
) -> Option<usize> {
    let find_root = |node: usize| {
        let mut root = node;
        while let Some(&parent) = parent.get(root) {
            if parent == root {
                break;
            }
            root = parent;
        }
        root
    };
    let first_root = find_root(first);
    let second_root = find_root(second);
    if first_root == second_root {
        return None;
    }
    let (low, high) = if rank.get(first_root).copied().unwrap_or(0)
        < rank.get(second_root).copied().unwrap_or(0)
    {
        (second_root, first_root)
    } else {
        (first_root, second_root)
    };
    if let Some(slot) = parent.get_mut(high) {
        *slot = low;
    }
    if rank.get(low).copied().unwrap_or(0) == rank.get(high).copied().unwrap_or(0)
        && let Some(slot) = rank.get_mut(low)
    {
        *slot = slot.saturating_add(1);
    }
    let merged = faces.get(low).copied().unwrap_or(0) | faces.get(high).copied().unwrap_or(0);
    if let Some(slot) = faces.get_mut(low) {
        *slot = merged;
    }
    Some(low)
}
