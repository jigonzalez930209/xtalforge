//! Fractional-space binning shared by the P1, neighbor and supercell paths.
//!
//! A bin count per axis comes from the perpendicular width of the cell, so a
//! difference of at most `size` Angstrom along that axis always lands in the
//! same or an adjacent bin, including periodic wraps.

use std::collections::HashMap;

use crate::{Error, Lattice};

/// Perpendicular width of every cell axis, in Angstrom.
pub(crate) fn perpendicular_widths(lattice: &Lattice) -> Result<[f64; 3], Error> {
    let reciprocal = lattice.reciprocal(false)?;
    Ok(reciprocal
        .matrix()
        .map(|row| 1.0 / row.iter().map(|value| value * value).sum::<f64>().sqrt()))
}

/// Bin count per axis for a search size in Angstrom.
pub(crate) fn bin_counts(widths: [f64; 3], size: f64) -> [usize; 3] {
    widths.map(|width| {
        let bins = (width / size).floor();
        if bins >= 1.0 { bins as usize } else { 1 }
    })
}

/// Bin index of a fractional coordinate on one axis.
pub(crate) fn bin_index(value: f64, count: usize) -> usize {
    let scaled = value.rem_euclid(1.0) * count as f64;
    let index = scaled.floor() as usize;
    index.min(count.saturating_sub(1))
}

/// Neighbouring bin offsets to visit; duplicates are removed for one or two
/// bins, where `-1` and `+1` wrap to the same bin.
pub(crate) fn offsets(count: usize) -> &'static [isize] {
    match count {
        0 | 1 => &[0],
        2 => &[-1, 0],
        _ => &[-1, 0, 1],
    }
}

/// Wraps a bin index with an offset into `0..count`.
pub(crate) fn wrap(index: usize, offset: isize, count: usize) -> usize {
    let count = count as isize;
    (index as isize + offset).rem_euclid(count) as usize
}

/// Fractional grid of atom indices, used for duplicate detection.
pub(crate) struct FractionalGrid<'a> {
    counts: [usize; 3],
    bins: HashMap<(&'a str, [usize; 3]), Vec<usize>>,
}

impl<'a> FractionalGrid<'a> {
    /// Creates an empty grid for the given cell and tolerance.
    pub(crate) fn new(widths: [f64; 3], tolerance: f64) -> Self {
        let counts = if tolerance > 0.0 {
            bin_counts(widths, tolerance)
        } else {
            [1, 1, 1]
        };
        Self {
            counts,
            bins: HashMap::new(),
        }
    }

    /// Visits every stored index that can be within the tolerance of `frac`.
    ///
    /// The visitor returns `false` to stop the walk.
    pub(crate) fn for_each_candidate(
        &self,
        element: &str,
        frac: [f64; 3],
        mut visit: impl FnMut(usize) -> bool,
    ) {
        let [count_x, count_y, count_z] = self.counts;
        let [frac_x, frac_y, frac_z] = frac;
        let base = [
            bin_index(frac_x, count_x),
            bin_index(frac_y, count_y),
            bin_index(frac_z, count_z),
        ];
        let [base_x, base_y, base_z] = base;
        for &offset_x in offsets(count_x) {
            let target_x = wrap(base_x, offset_x, count_x);
            for &offset_y in offsets(count_y) {
                let target_y = wrap(base_y, offset_y, count_y);
                for &offset_z in offsets(count_z) {
                    let target_z = wrap(base_z, offset_z, count_z);
                    let Some(bucket) = self.bins.get(&(element, [target_x, target_y, target_z]))
                    else {
                        continue;
                    };
                    for &index in bucket {
                        if !visit(index) {
                            return;
                        }
                    }
                }
            }
        }
    }

    /// Inserts an atom index at its fractional position.
    pub(crate) fn insert(&mut self, element: &'a str, frac: [f64; 3], index: usize) {
        let [count_x, count_y, count_z] = self.counts;
        let [frac_x, frac_y, frac_z] = frac;
        let key = (
            element,
            [
                bin_index(frac_x, count_x),
                bin_index(frac_y, count_y),
                bin_index(frac_z, count_z),
            ],
        );
        self.bins.entry(key).or_default().push(index);
    }
}
