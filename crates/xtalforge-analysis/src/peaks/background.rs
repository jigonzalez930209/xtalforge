//! SNIP background envelope for measured patterns (Task 7.3.1).
//!
//! The classical statistics-sensitive non-linear iterative peak-clipping
//! algorithm (Ryan et al., 1988) works on the log-log transform
//! `t = ln(ln(sqrt(I + 1) + 1) + 1)`, clips every point against the average of
//! its neighbours with a window that shrinks from `window` to one, and
//! inverts the transform. Peaks are narrowed at every pass while a smooth
//! baseline is preserved, which is what separates a weak line from a sloping
//! amorphous background.

use xtalforge_io::xy::PatternPoint;

/// Background envelope of a pattern, in intensity units.
///
/// A zero window (or an empty pattern) returns the input intensities.
#[must_use]
pub fn background_envelope(points: &[PatternPoint], window: usize) -> Vec<f64> {
    let count = points.len();
    if count == 0 || window == 0 {
        return points.iter().map(|point| point.intensity).collect();
    }
    let mut clipped: Vec<f64> = points
        .iter()
        .map(|point| {
            let value = point.intensity.max(0.0);
            ((value.sqrt() + 1.0).ln() + 1.0).ln()
        })
        .collect();
    for step in (1..=window.min(count / 2)).rev() {
        let previous = clipped.clone();
        for index in 0..count {
            let left = index
                .checked_sub(step)
                .and_then(|slot| previous.get(slot))
                .copied()
                .unwrap_or_else(|| previous.get(index).copied().unwrap_or(0.0));
            let right = previous
                .get(index + step)
                .copied()
                .or_else(|| previous.get(index).copied())
                .unwrap_or(0.0);
            let average = 0.5 * (left + right);
            if let Some(cell) = clipped.get_mut(index)
                && average < *cell
            {
                *cell = average;
            }
        }
    }
    clipped
        .into_iter()
        .map(|value| {
            let inner = value.exp() - 1.0;
            let root = inner.exp() - 1.0;
            (root * root - 1.0).max(0.0)
        })
        .collect()
}
