//! Small vector, matrix and linear-solver helpers shared by this crate.
//!
//! The crate deliberately avoids a hard dependency on a linear-algebra crate:
//! the systems it solves are at most `(n + 1) x (n + 1)` for charge
//! equilibration and `4 x 4` for an equation of state, and the geometry code
//! only needs 3-vectors and 3x3 matrices.  Everything here is deterministic
//! and free of unchecked indexing.

/// Component-wise sum.
pub(crate) fn add(first: [f64; 3], second: [f64; 3]) -> [f64; 3] {
    [
        first[0] + second[0],
        first[1] + second[1],
        first[2] + second[2],
    ]
}

/// Component-wise difference.
pub(crate) fn sub(first: [f64; 3], second: [f64; 3]) -> [f64; 3] {
    [
        first[0] - second[0],
        first[1] - second[1],
        first[2] - second[2],
    ]
}

/// Scalar multiple.
pub(crate) fn scale(vector: [f64; 3], factor: f64) -> [f64; 3] {
    [vector[0] * factor, vector[1] * factor, vector[2] * factor]
}

/// Dot product.
pub(crate) fn dot(first: [f64; 3], second: [f64; 3]) -> f64 {
    first[0] * second[0] + first[1] * second[1] + first[2] * second[2]
}

/// Cross product.
pub(crate) fn cross(first: [f64; 3], second: [f64; 3]) -> [f64; 3] {
    [
        first[1] * second[2] - first[2] * second[1],
        first[2] * second[0] - first[0] * second[2],
        first[0] * second[1] - first[1] * second[0],
    ]
}

/// Euclidean norm.
pub(crate) fn norm(vector: [f64; 3]) -> f64 {
    dot(vector, vector).sqrt()
}

/// Unit vector, or `None` for a vector shorter than `1e-14`.
pub(crate) fn unit(vector: [f64; 3]) -> Option<[f64; 3]> {
    let length = norm(vector);
    if length < 1e-14 {
        None
    } else {
        Some(scale(vector, 1.0 / length))
    }
}

/// Identity matrix.
pub(crate) fn identity3() -> [[f64; 3]; 3] {
    [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]]
}

/// Matrix product.
pub(crate) fn mat3_mul(first: [[f64; 3]; 3], second: [[f64; 3]; 3]) -> [[f64; 3]; 3] {
    let mut result = [[0.0; 3]; 3];
    for (row, output) in result.iter_mut().enumerate() {
        let Some(first_row) = first.get(row) else {
            continue;
        };
        for (column, slot) in output.iter_mut().enumerate() {
            let mut value = 0.0;
            for (index, first_value) in first_row.iter().enumerate() {
                let second_value = second
                    .get(index)
                    .and_then(|values| values.get(column))
                    .copied()
                    .unwrap_or(0.0);
                value += first_value * second_value;
            }
            *slot = value;
        }
    }
    result
}

/// Matrix transpose.
pub(crate) fn mat3_transpose(matrix: [[f64; 3]; 3]) -> [[f64; 3]; 3] {
    [
        [matrix[0][0], matrix[1][0], matrix[2][0]],
        [matrix[0][1], matrix[1][1], matrix[2][1]],
        [matrix[0][2], matrix[1][2], matrix[2][2]],
    ]
}

/// Row vector times matrix, `v @ m`.
pub(crate) fn vec_mat(vector: [f64; 3], matrix: [[f64; 3]; 3]) -> [f64; 3] {
    [
        vector[0] * matrix[0][0] + vector[1] * matrix[1][0] + vector[2] * matrix[2][0],
        vector[0] * matrix[0][1] + vector[1] * matrix[1][1] + vector[2] * matrix[2][1],
        vector[0] * matrix[0][2] + vector[1] * matrix[1][2] + vector[2] * matrix[2][2],
    ]
}

/// Component-wise matrix sum.
pub(crate) fn mat3_add(first: [[f64; 3]; 3], second: [[f64; 3]; 3]) -> [[f64; 3]; 3] {
    let mut result = [[0.0; 3]; 3];
    for (row, output) in result.iter_mut().enumerate() {
        for (column, slot) in output.iter_mut().enumerate() {
            let a = first
                .get(row)
                .and_then(|values| values.get(column))
                .copied()
                .unwrap_or(0.0);
            let b = second
                .get(row)
                .and_then(|values| values.get(column))
                .copied()
                .unwrap_or(0.0);
            *slot = a + b;
        }
    }
    result
}

/// Scalar multiple of a matrix.
pub(crate) fn mat3_scale(matrix: [[f64; 3]; 3], factor: f64) -> [[f64; 3]; 3] {
    let mut result = [[0.0; 3]; 3];
    for (row, output) in result.iter_mut().enumerate() {
        for (column, slot) in output.iter_mut().enumerate() {
            let value = matrix
                .get(row)
                .and_then(|values| values.get(column))
                .copied()
                .unwrap_or(0.0);
            *slot = value * factor;
        }
    }
    result
}

/// Solves `matrix * x = right` by Gaussian elimination with partial pivoting.
///
/// `matrix` and `right` are modified in place.  Returns `None` when the
/// system is not square or is singular to working precision; callers turn
/// that into their crate error.
pub(crate) fn solve_linear(matrix: &mut [Vec<f64>], right: &mut [f64]) -> Option<Vec<f64>> {
    let n = right.len();
    if matrix.len() != n {
        return None;
    }
    for column in 0..n {
        let mut pivot = column;
        let mut best = matrix.get(column)?.get(column)?.abs();
        for row in (column + 1)..n {
            let value = matrix.get(row)?.get(column)?.abs();
            if value > best {
                best = value;
                pivot = row;
            }
        }
        if best < 1e-13 {
            return None;
        }
        if pivot != column {
            matrix.swap(column, pivot);
            right.swap(column, pivot);
        }
        let pivot_value = *matrix.get(column)?.get(column)?;
        if pivot_value == 0.0 {
            return None;
        }
        let pivot_row = matrix.get(column)?.clone();
        for row in (column + 1)..n {
            let factor = *matrix.get(row)?.get(column)? / pivot_value;
            if factor == 0.0 {
                continue;
            }
            let values = matrix.get_mut(row)?;
            for index in (column + 1)..n {
                let pivot_component = *pivot_row.get(index)?;
                let slot = values.get_mut(index)?;
                *slot -= factor * pivot_component;
            }
            if let Some(slot) = values.get_mut(column) {
                *slot = 0.0;
            }
            let pivot_right = *right.get(column)?;
            if let Some(value) = right.get_mut(row) {
                *value -= factor * pivot_right;
            }
        }
    }

    let mut solution = vec![0.0; n];
    for row in (0..n).rev() {
        let values = matrix.get(row)?;
        let mut total = *right.get(row)?;
        for index in (row + 1)..n {
            let coefficient = *values.get(index)?;
            total -= coefficient * solution.get(index)?;
        }
        let diagonal = *values.get(row)?;
        if diagonal.abs() < 1e-13 {
            return None;
        }
        if let Some(slot) = solution.get_mut(row) {
            *slot = total / diagonal;
        }
    }
    Some(solution)
}
