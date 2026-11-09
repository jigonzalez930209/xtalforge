//! Small dense linear algebra for the Pawley fit (Task 7.6.2).

#![allow(clippy::indexing_slicing)] // fixed-size local systems

/// Solves a square system by Gauss-Jordan elimination with partial pivoting.
pub(crate) fn solve(matrix: &[Vec<f64>], rhs: &[f64]) -> Option<Vec<f64>> {
    let n = matrix.len();
    let mut augmented: Vec<Vec<f64>> = matrix
        .iter()
        .zip(rhs.iter())
        .map(|(row, value)| {
            let mut line = row.clone();
            line.push(*value);
            line
        })
        .collect();
    eliminate(&mut augmented, n)?;
    Some(augmented.iter().map(|row| row[n]).collect())
}

/// Inverts a square matrix, or returns `None` when it is singular.
pub(crate) fn invert(matrix: &[Vec<f64>]) -> Option<Vec<Vec<f64>>> {
    let n = matrix.len();
    let mut augmented: Vec<Vec<f64>> = matrix
        .iter()
        .enumerate()
        .map(|(row, values)| {
            let mut line = values.clone();
            line.extend((0..n).map(|column| if column == row { 1.0 } else { 0.0 }));
            line
        })
        .collect();
    eliminate(&mut augmented, n)?;
    Some(
        augmented
            .iter()
            .map(|row| row[n..].to_vec())
            .collect::<Vec<_>>(),
    )
}

fn eliminate(augmented: &mut [Vec<f64>], n: usize) -> Option<()> {
    for column in 0..n {
        let pivot = (column..n).max_by(|a, b| {
            augmented
                .get(*a)
                .and_then(|row| row.get(column))
                .unwrap_or(&0.0)
                .abs()
                .partial_cmp(
                    &augmented
                        .get(*b)
                        .and_then(|row| row.get(column))
                        .unwrap_or(&0.0)
                        .abs(),
                )
                .unwrap_or(std::cmp::Ordering::Equal)
        })?;
        let pivot_value = *augmented.get(pivot)?.get(column)?;
        if pivot_value.abs() < 1e-14 {
            return None;
        }
        augmented.swap(column, pivot);
        let divisor = *augmented.get(column)?.get(column)?;
        for value in augmented.get_mut(column)?.iter_mut() {
            *value /= divisor;
        }
        let pivot_row = augmented.get(column)?.clone();
        for (row, line) in augmented.iter_mut().enumerate() {
            if row == column {
                continue;
            }
            let factor = *line.get(column)?;
            for (slot, pivot) in line.iter_mut().zip(pivot_row.iter()) {
                *slot -= factor * pivot;
            }
        }
    }
    Some(())
}
