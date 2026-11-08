//! Periodic-image helpers shared by QEq, EQeq and Ewald.
//!
//! The image range is derived from the perpendicular width of the cell in
//! each direction, the same rule Crystal-Builder uses: an atom sees every
//! image within the cutoff, however skewed the cell is.

use xtalforge_core::Lattice;

use crate::math::{cross, norm};

/// Perpendicular widths of the cell along `a`, `b` and `c`, in Angstrom.
pub(crate) fn widths(lattice: &Lattice) -> [f64; 3] {
    let matrix = lattice.matrix();
    let volume = lattice.volume();
    let first = norm(cross(matrix[1], matrix[2]));
    let second = norm(cross(matrix[2], matrix[0]));
    let third = norm(cross(matrix[0], matrix[1]));
    [
        volume / first.max(1e-12),
        volume / second.max(1e-12),
        volume / third.max(1e-12),
    ]
}

/// Integer image range needed to reach `cutoff` in every direction.
pub(crate) fn image_range(lattice: &Lattice, cutoff: f64) -> [i32; 3] {
    let widths = widths(lattice);
    [
        (cutoff / widths[0] + 0.5).ceil().max(1.0) as i32,
        (cutoff / widths[1] + 0.5).ceil().max(1.0) as i32,
        (cutoff / widths[2] + 0.5).ceil().max(1.0) as i32,
    ]
}
