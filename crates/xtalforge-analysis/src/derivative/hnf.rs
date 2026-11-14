//! Hermite-normal-form supercell helper and coordinate expansion.

use std::collections::BTreeSet;

use xtalforge_core::{Error, Lattice};

/// Largest HNF index accepted by [`hnf_supercells`].
pub const MAX_HNF_INDEX: u32 = 128;

/// Enumerates the HNF supercells of one index in deterministic order.
///
/// The matrices are upper triangular,
///
/// ```text
/// [[a, b, c],
///  [0, d, e],
///  [0, 0, f]]
/// ```
///
/// with `a d f = index`, `0 <= b < d`, `0 <= c < f` and `0 <= e < f`. The
/// order is `(a, d, f, b, c, e)`.
///
/// # Errors
///
/// Returns [`Error::InvalidInput`] for a zero index or one above
/// [`MAX_HNF_INDEX`].
pub fn hnf_supercells(index: u32) -> Result<Vec<[[i32; 3]; 3]>, Error> {
    if index == 0 || index > MAX_HNF_INDEX {
        return Err(Error::InvalidInput(format!(
            "HNF index must be between 1 and {MAX_HNF_INDEX}, got {index}"
        )));
    }
    let mut matrices = Vec::new();
    for a in 1..=index {
        if !index.is_multiple_of(a) {
            continue;
        }
        let rest = index / a;
        for d in 1..=rest {
            if !rest.is_multiple_of(d) {
                continue;
            }
            let f = rest / d;
            for b in 0..d {
                for c in 0..f {
                    for e in 0..f {
                        matrices.push([
                            [a as i32, b as i32, c as i32],
                            [0, d as i32, e as i32],
                            [0, 0, f as i32],
                        ]);
                    }
                }
            }
        }
    }
    matrices.sort_by_key(|matrix| {
        let [[a, b, c], [_, d, e], [_, _, f]] = *matrix;
        (a, d, f, b, c, e)
    });
    Ok(matrices)
}

/// Expands one fractional position into its supercell images.
pub(super) fn supercell_positions(
    original: [f64; 3],
    matrix: [[i32; 3]; 3],
) -> Result<Vec<[f64; 3]>, Error> {
    let integer = matrix.map(|row| row.map(f64::from));
    let inverse = Lattice::new(integer)?.inverse();
    let [[a, _, _], [_, d, _], [_, _, f]] = matrix;
    let [x, y, z] = original;
    let mut positions = Vec::with_capacity((a * d * f) as usize);
    let mut seen = BTreeSet::new();
    for first in 0..a {
        for second in 0..d {
            for third in 0..f {
                let vector = [
                    x + f64::from(first),
                    y + f64::from(second),
                    z + f64::from(third),
                ];
                let frac = Lattice::wrap(multiply_row(vector, inverse));
                let key = position_key(frac);
                if seen.insert(key) {
                    positions.push(frac);
                }
            }
        }
    }
    Ok(positions)
}

fn position_key(frac: [f64; 3]) -> [i64; 3] {
    frac.map(|value| (value * 1_000_000.0).round() as i64)
}

fn multiply_row(vector: [f64; 3], matrix: [[f64; 3]; 3]) -> [f64; 3] {
    let [x, y, z] = vector;
    let [[a, b, c], [d, e, f], [g, h, i]] = matrix;
    [
        x * a + y * d + z * g,
        x * b + y * e + z * h,
        x * c + y * f + z * i,
    ]
}
