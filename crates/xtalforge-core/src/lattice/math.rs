//! Small fixed-size matrix helpers for the lattice module.

pub(super) fn norm(vector: [f64; 3]) -> f64 {
    vector.iter().map(|value| value * value).sum::<f64>().sqrt()
}

pub(super) fn subtract(first: [f64; 3], second: [f64; 3]) -> [f64; 3] {
    let [a, b, c] = first;
    let [d, e, f] = second;
    [a - d, b - e, c - f]
}

fn add_row(first: [f64; 3], second: [f64; 3]) -> [f64; 3] {
    let [a, b, c] = first;
    let [d, e, f] = second;
    [a + d, b + e, c + f]
}

pub(super) fn add(first: [[f64; 3]; 3], second: [[f64; 3]; 3]) -> [[f64; 3]; 3] {
    let [a, b, c] = first;
    let [d, e, f] = second;
    [add_row(a, d), add_row(b, e), add_row(c, f)]
}

fn dot(first: [f64; 3], second: [f64; 3]) -> f64 {
    let [a, b, c] = first;
    let [d, e, f] = second;
    a * d + b * e + c * f
}

pub(super) fn multiply(first: [[f64; 3]; 3], second: [[f64; 3]; 3]) -> [[f64; 3]; 3] {
    let [a, b, c] = first;
    let [[d0, d1, d2], [e0, e1, e2], [f0, f1, f2]] = second;
    let column0 = [d0, e0, f0];
    let column1 = [d1, e1, f1];
    let column2 = [d2, e2, f2];
    [
        [dot(a, column0), dot(a, column1), dot(a, column2)],
        [dot(b, column0), dot(b, column1), dot(b, column2)],
        [dot(c, column0), dot(c, column1), dot(c, column2)],
    ]
}

pub(super) fn multiply_row(vector: [f64; 3], matrix: [[f64; 3]; 3]) -> [f64; 3] {
    let [x, y, z] = vector;
    let [[a, b, c], [d, e, f], [g, h, i]] = matrix;
    [
        x * a + y * d + z * g,
        x * b + y * e + z * h,
        x * c + y * f + z * i,
    ]
}

pub(super) fn transpose(matrix: [[f64; 3]; 3]) -> [[f64; 3]; 3] {
    let [[a, b, c], [d, e, f], [g, h, i]] = matrix;
    [[a, d, g], [b, e, h], [c, f, i]]
}

pub(super) fn determinant(matrix: [[f64; 3]; 3]) -> f64 {
    let [[a, b, c], [d, e, f], [g, h, i]] = matrix;
    a * (e * i - f * h) - b * (d * i - f * g) + c * (d * h - e * g)
}

pub(super) fn inverse(matrix: [[f64; 3]; 3]) -> [[f64; 3]; 3] {
    let [[a, b, c], [d, e, f], [g, h, i]] = matrix;
    let det = determinant(matrix);
    if det.abs() < super::MIN_VOLUME {
        return [[f64::NAN; 3]; 3];
    }
    [
        [
            (e * i - f * h) / det,
            (c * h - b * i) / det,
            (b * f - c * e) / det,
        ],
        [
            (f * g - d * i) / det,
            (a * i - c * g) / det,
            (c * d - a * f) / det,
        ],
        [
            (d * h - e * g) / det,
            (b * g - a * h) / det,
            (a * e - b * d) / det,
        ],
    ]
}

pub(super) fn angle_between(first: [f64; 3], second: [f64; 3], scale: f64) -> f64 {
    let cosine = dot(first, second) / scale;
    cosine.clamp(-1.0, 1.0).acos().to_degrees()
}
